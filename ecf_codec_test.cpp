/**********************************************************************************************************************
 * Copyright (c) Prophesee S.A.                                                                                       *
 *                                                                                                                    *
 * Licensed under the Apache License, Version 2.0 (the "License");                                                    *
 * you may not use this file except in compliance with the License.                                                   *
 * You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0                                 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed   *
 * on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                      *
 * See the License for the specific language governing permissions and limitations under the License.                 *
 **********************************************************************************************************************/

#include <array>
#include <utility>
#include <random>
#include <gtest/gtest.h>

#include "ecf_codec.h"

using namespace ECF;

namespace {
template<typename... T>
auto make_array(T &&...t) -> std::array<std::common_type_t<T...>, sizeof...(t)> {
    return {std::forward<T>(t)...};
}

template<typename T>
std::vector<T> &operator+=(std::vector<T> &v1, const std::vector<T> &v2) {
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

template<typename T, size_t N>
std::vector<T> &operator+=(std::vector<T> &v, const std::array<T, N> &a) {
    v.insert(v.end(), a.begin(), a.end());
    return v;
}

template<typename T>
std::vector<T> &operator+=(std::vector<T> &v, const T &t) {
    v.push_back(t);
    return v;
}
} // namespace

#define EXPLODE16(v) std::uint8_t(v & 0xff), std::uint8_t((v >> 8) & 0xff)
#define EXPLODE32(v)                                                                       \
    std::uint8_t(v & 0xff), std::uint8_t((v >> 8) & 0xff), std::uint8_t((v >> 16) & 0xff), \
        std::uint8_t((v >> 24) & 0xff)
#define EXPLODE64(v)                                                                                    \
    std::uint8_t(v & 0xff), std::uint8_t((v >> 8) & 0xff), std::uint8_t((v >> 16) & 0xff),              \
        std::uint8_t((v >> 24) & 0xff), std::uint8_t((v >> 32) & 0xff), std::uint8_t((v >> 40) & 0xff), \
        std::uint8_t((v >> 48) & 0xff), std::uint8_t((v >> 56) & 0xff)

#define HEADER(n, pyxt, pxt) EXPLODE32((((n) << 2) | (((pyxt) ? 1 : 0) << 1) | ((pxt) ? 1 : 0)))
#define ABSOLUTE_T(t) EXPLODE64(t)
#define RLE_T(t, c) std::uint8_t(((t) << 4) | c)
#define DIFF_T(dt) std::uint8_t(0b11110000 | ((dt)&0b1111))
#define RLE_Y(y, c) EXPLODE16((((y) << 5) | c))
#define RLE_P(p, c) std::uint8_t(((p) << 7) | c)

#define PACK_Y_X_AND_P(e) (((e.y) << 12) | ((e.x) << 1) | (e.p))
#define PACK_YS_XS_AND_PS_4(e0, e1, e2, e3)                                                 \
    EXPLODE32((((PACK_Y_X_AND_P(e0)) << 8) | ((PACK_Y_X_AND_P(e1)) & 0xff))),               \
        EXPLODE32(((((PACK_Y_X_AND_P(e1)) >> 8) << 16) | ((PACK_Y_X_AND_P(e2)) & 0xffff))), \
        EXPLODE32(((((PACK_Y_X_AND_P(e2)) >> 16) << 24) | ((PACK_Y_X_AND_P(e3)) & 0xffffff)))
#define PACK_YS_XS_AND_PS_3(e0, e1, e2) PACK_YS_XS_AND_PS_4(e0, e1, e2, EMPTY_EVENT)
#define PACK_YS_XS_AND_PS_2(e0, e1) PACK_YS_XS_AND_PS_4(e0, e1, EMPTY_EVENT, EMPTY_EVENT)
#define PACK_YS_XS_AND_PS_1(e0) PACK_YS_XS_AND_PS_4(e0, EMPTY_EVENT, EMPTY_EVENT, EMPTY_EVENT)

#define PACK_X_AND_P(e) (((e.x) << 1) | (e.p))
#define PACK_XS_AND_PS_4(e0, e1, e2, e3)                                                   \
    EXPLODE16((((PACK_X_AND_P(e0)) << 4) | ((PACK_X_AND_P(e1)) & 0b1111))),                \
        EXPLODE16(((((PACK_X_AND_P(e1)) >> 4) << 8) | ((PACK_X_AND_P(e2)) & 0b11111111))), \
        EXPLODE16(((((PACK_X_AND_P(e2)) >> 8) << 12) | ((PACK_X_AND_P(e3)) & 0b111111111111)))
#define PACK_XS_AND_PS_3(e0, e1, e2) PACK_XS_AND_PS_4(e0, e1, e2, EMPTY_EVENT)
#define PACK_XS_AND_PS_2(e0, e1) PACK_XS_AND_PS_4(e0, e1, EMPTY_EVENT, EMPTY_EVENT)
#define PACK_XS_AND_PS_1(e0) PACK_XS_AND_PS_4(e0, EMPTY_EVENT, EMPTY_EVENT, EMPTY_EVENT)

#define MASK_X(x, m) EXPLODE16((((x) << 5) | (m)))

#define EMPTY_EVENT (EventCD{(unsigned short)0, (unsigned short)0, (short)0, (long long)0})

TEST(Encoder, encode_too_many_events) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    std::vector<EventCD> events(100000);
    const size_t num_events = events.size();

    std::vector<std::uint8_t> buffer(1000000);
    ASSERT_THROW(encoder(events.data(), events.data() + num_events, buffer.data()), std::runtime_error);
}

TEST(Decoder, decode_too_many_events) {
    Decoder decoder;

    std::vector<EventCD> expected_events(100000);
    const size_t num_expected_events = expected_events.size();
    auto buffer                      = make_array(
        // 4 bytes for header
        HEADER(num_expected_events, false, false),
        // 8 bytes for absolute t
        ABSOLUTE_T(expected_events[0].t),
        // 1 byte for 1 timestamp RLE encoded
        RLE_T(expected_events[0].t, 1),
        // 3*4 bytes for x and p
        PACK_YS_XS_AND_PS_1(expected_events[0]));

    std::vector<EventCD> events(100000);
    ASSERT_THROW(decoder(buffer.data(), buffer.data() + buffer.size(), events.data()), std::runtime_error);
}

TEST(Encoder, encode_one_empty_event) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    std::vector<EventCD> events{EMPTY_EVENT};
    const size_t num_events = events.size();
    auto expected_buffer    = make_array(
        // 4 bytes for header
        HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps),
        // 8 bytes for absolute t
        ABSOLUTE_T(events[0].t),
        // 1 byte for 1 timestamp RLE encoded
        RLE_T(events[0].t, 1),
        // 3*4 bytes for x and p
        PACK_YS_XS_AND_PS_1(events[0]));

    std::vector<std::uint8_t> buffer(1000000);
    ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
    for (size_t i = 0; i < expected_buffer.size(); ++i) {
        EXPECT_EQ(expected_buffer[i], buffer[i]);
    }
}

TEST(Decoder, decode_one_empty_event) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    auto expected_events             = make_array(EMPTY_EVENT);
    const size_t num_expected_events = expected_events.size();
    auto buffer                      = make_array(
        // 4 bytes for header
        HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps),
        // 8 bytes for absolute t
        ABSOLUTE_T(expected_events[0].t),
        // 1 byte for 1 timestamp RLE encoded
        RLE_T(expected_events[0].t, 1),
        // 3*4 bytes for x and p
        PACK_YS_XS_AND_PS_1(expected_events[0]));

    EXPECT_EQ(num_expected_events * sizeof(EventCD),
              decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

    std::vector<EventCD> events(100000);
    ASSERT_EQ(expected_events.size() * sizeof(EventCD),
              decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
    for (size_t i = 0; i < expected_events.size(); ++i) {
        EXPECT_EQ(expected_events[i].x, events[i].x);
        EXPECT_EQ(expected_events[i].y, events[i].y);
        EXPECT_EQ(expected_events[i].p, events[i].p);
        EXPECT_EQ(expected_events[i].t, events[i].t);
    }
}

TEST(Encoder, encode_one_event) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    auto events             = make_array(EventCD{1, 2, 0, 4});
    const size_t num_events = events.size();
    auto expected_buffer    = make_array(
        // 4 bytes for header
        HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps),
        // 8 bytes for absolute t
        ABSOLUTE_T(events[0].t),
        // 1 byte for 1 timestamp RLE encoded
        RLE_T(0, 1),
        // 3*4 bytes for x and p
        PACK_YS_XS_AND_PS_1(events[0]));

    std::vector<std::uint8_t> buffer(1000000);
    ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
    for (size_t i = 0; i < expected_buffer.size(); ++i) {
        EXPECT_EQ(expected_buffer[i], buffer[i]);
    }
}

TEST(Decoder, decode_one_event) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    auto expected_events             = make_array(EventCD{1, 2, 0, 4});
    const size_t num_expected_events = expected_events.size();
    auto buffer                      = make_array(
        // 4 bytes for header
        HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps),
        // 8 bytes for absolute t
        ABSOLUTE_T(expected_events[0].t),
        // 1 byte for 1 timestamp RLE encoded
        RLE_T(0, 1),
        // 3*4 bytes for x and p
        PACK_YS_XS_AND_PS_1(expected_events[0]));

    EXPECT_EQ(num_expected_events * sizeof(EventCD),
              decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

    std::vector<EventCD> events(100000);
    ASSERT_EQ(expected_events.size() * sizeof(EventCD),
              decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
    for (size_t i = 0; i < expected_events.size(); ++i) {
        EXPECT_EQ(expected_events[i].x, events[i].x);
        EXPECT_EQ(expected_events[i].y, events[i].y);
        EXPECT_EQ(expected_events[i].p, events[i].p);
        EXPECT_EQ(expected_events[i].t, events[i].t);
    }
}

TEST(Encoder, encode_up_to_4_empty_events) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EMPTY_EVENT;
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        expected_buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            expected_buffer += make_array(PACK_XS_AND_PS_2(events[0], events[1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_XS_AND_PS_3(events[0], events[1], events[2]));
            break;
        case 4:
            expected_buffer += make_array(PACK_XS_AND_PS_4(events[0], events[1], events[2], events[3]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_empty_events) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EMPTY_EVENT;
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            buffer += make_array(PACK_XS_AND_PS_2(expected_events[0], expected_events[1]));
            break;
        case 3:
            buffer += make_array(PACK_XS_AND_PS_3(expected_events[0], expected_events[1], expected_events[2]));
            break;
        case 4:
            buffer += make_array(
                PACK_XS_AND_PS_4(expected_events[0], expected_events[1], expected_events[2], expected_events[3]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_up_to_4_events_same_t_rle_same_y_rle_xp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EventCD{(unsigned short)(i * 10), 23, 1, 10000};
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        expected_buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            expected_buffer += make_array(PACK_XS_AND_PS_2(events[0], events[1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_XS_AND_PS_3(events[0], events[1], events[2]));
            break;
        case 4:
            expected_buffer += make_array(PACK_XS_AND_PS_4(events[0], events[1], events[2], events[3]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_events_same_t_rle_same_y_rle_xp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EventCD{(unsigned short)(i * 10), 23, 1, 10000};
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            buffer += make_array(PACK_XS_AND_PS_2(expected_events[0], expected_events[1]));
            break;
        case 3:
            buffer += make_array(PACK_XS_AND_PS_3(expected_events[0], expected_events[1], expected_events[2]));
            break;
        case 4:
            buffer += make_array(
                PACK_XS_AND_PS_4(expected_events[0], expected_events[1], expected_events[2], expected_events[3]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_up_to_4_events_diff_t_rle_same_y_rle_xp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EventCD{(unsigned short)(i * 10), 23, 1, (long long)i};
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for each timestamp RLE encoded
        expected_buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            expected_buffer += RLE_T(events[i].t - events[i - 1].t, 1);
        }
        // 2 bytes for 4 equal y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            expected_buffer += make_array(PACK_XS_AND_PS_2(events[0], events[1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_XS_AND_PS_3(events[0], events[1], events[2]));
            break;
        case 4:
            expected_buffer += make_array(PACK_XS_AND_PS_4(events[0], events[1], events[2], events[3]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_events_diff_t_rle_same_y_rle_xp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EventCD{(unsigned short)(i * 10), 23, 1, (long long)i};
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for each timestamp RLE encoded
        buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            buffer += RLE_T(expected_events[i].t - expected_events[i - 1].t, 1);
        }
        // 2 bytes for 4 equal y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            buffer += make_array(PACK_XS_AND_PS_2(expected_events[0], expected_events[1]));
            break;
        case 3:
            buffer += make_array(PACK_XS_AND_PS_3(expected_events[0], expected_events[1], expected_events[2]));
            break;
        case 4:
            buffer += make_array(
                PACK_XS_AND_PS_4(expected_events[0], expected_events[1], expected_events[2], expected_events[3]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_up_to_4_events_big_diff_t_rle_same_y_rle_xp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    const long long diff_t = 5217;
    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EventCD{(unsigned short)(i * 10), 23, 1, (long long)((i + 1) * (i + 2 * diff_t) / 2)};
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for each timestamp RLE encoded
        expected_buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            // multiple bytes for diff to previous timestamp
            expected_buffer += DIFF_T(diff_t + i);
            expected_buffer += DIFF_T((diff_t + i) >> 4);
            expected_buffer += DIFF_T((diff_t + i) >> 8);
            expected_buffer += DIFF_T((diff_t + i) >> 12);
            expected_buffer += RLE_T(0, 1);
        }
        // 2 bytes for 4 equal y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            expected_buffer += make_array(PACK_XS_AND_PS_2(events[0], events[1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_XS_AND_PS_3(events[0], events[1], events[2]));
            break;
        case 4:
            expected_buffer += make_array(PACK_XS_AND_PS_4(events[0], events[1], events[2], events[3]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_events_big_diff_t_rle_same_y_rle_xp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    const long long diff_t = 5217;
    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EventCD{(unsigned short)(i * 10), 23, 1, (long long)((i + 1) * (i + 2 * diff_t) / 2)};
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for each timestamp RLE encoded
        buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            // multiple bytes for diff to previous timestamp
            buffer += DIFF_T(diff_t + i);
            buffer += DIFF_T((diff_t + i) >> 4);
            buffer += DIFF_T((diff_t + i) >> 8);
            buffer += DIFF_T((diff_t + i) >> 12);
            buffer += RLE_T(0, 1);
        }
        // 2 bytes for 4 equal y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, j));
        // 2*2 bytes for x and p
        switch (j) {
        case 2:
            buffer += make_array(PACK_XS_AND_PS_2(expected_events[0], expected_events[1]));
            break;
        case 3:
            buffer += make_array(PACK_XS_AND_PS_3(expected_events[0], expected_events[1], expected_events[2]));
            break;
        case 4:
            buffer += make_array(
                PACK_XS_AND_PS_4(expected_events[0], expected_events[1], expected_events[2], expected_events[3]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_up_to_4_events_diff_t_rle_xyp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EventCD{(unsigned short)(29 + i), (unsigned short)i, 0, (long long)i};
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for each timestamp RLE encoded
        expected_buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            expected_buffer += RLE_T(events[i].t - events[i - 1].t, 1);
        }
        // 3*4 bytes for x, y and p
        switch (j) {
        case 2:
            expected_buffer += make_array(PACK_YS_XS_AND_PS_2(events[0], events[1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_YS_XS_AND_PS_3(events[0], events[1], events[2]));
            break;
        case 4:
            expected_buffer += make_array(PACK_YS_XS_AND_PS_4(events[0], events[1], events[2], events[3]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_events_diff_t_rle_xyp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EventCD{(unsigned short)(29 + i), (unsigned short)i, 0, (long long)i};
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for each timestamp RLE encoded
        buffer += RLE_T(0, 1);
        for (size_t i = 1; i < j; ++i) {
            buffer += RLE_T(expected_events[i].t - expected_events[i - 1].t, 1);
        }
        // 3*4 bytes for x, y and p
        switch (j) {
        case 2:
            buffer += make_array(PACK_YS_XS_AND_PS_2(expected_events[0], expected_events[1]));
            break;
        case 3:
            buffer += make_array(PACK_YS_XS_AND_PS_3(expected_events[0], expected_events[1], expected_events[2]));
            break;
        case 4:
            buffer += make_array(
                PACK_YS_XS_AND_PS_4(expected_events[0], expected_events[1], expected_events[2], expected_events[3]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_up_to_6_events_same_t_rle_same_y_rle_x_masked) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = false;

    for (size_t j = 2; j <= 6; ++j) {
        std::vector<EventCD> events;
        for (size_t i = 0; i < j; ++i) {
            events += EventCD{(unsigned short)(19 + i), (unsigned short)23, 0, (long long)13875};
        }

        const size_t num_events = events.size();
        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        expected_buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, j));
        // 2 byte for x masked
        expected_buffer += make_array(MASK_X(events[0].x, (1 << 5) - (1 << (6 - j))));
        // 1 byte for 1 polarity RLE encoded
        expected_buffer += RLE_P(events[0].p, j);

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_up_to_4_events_same_t_rle_same_y_rle_x_masked) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = false;

    for (size_t j = 2; j <= 4; ++j) {
        std::vector<EventCD> expected_events;
        for (size_t i = 0; i < j; ++i) {
            expected_events += EventCD{(unsigned short)(19 + i), (unsigned short)23, 0, (long long)13875};
        }

        const size_t num_expected_events = expected_events.size();
        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for 1 timestamp RLE encoded
        buffer += RLE_T(0, j);
        // 2 bytes for 4 equal y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, j));
        // 2 byte for x masked
        buffer += make_array(MASK_X(expected_events[0].x, (1 << 5) - (1 << (6 - j))));
        // 1 byte for 1 polarity RLE encoded
        buffer += RLE_P(expected_events[0].p, j);

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_more_than_15_events_same_t_rle_xyp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279), dy(0, 719);
    mt_rand.seed(42);

    auto num_events_list = {16, 25, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_events : num_events_list) {
        std::vector<EventCD> events;
        for (int i = 0; i < num_events; ++i) {
            events += EventCD{(unsigned short)(dx(mt_rand)), (unsigned short)(dy(mt_rand)), 1, 10000};
        }

        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for max count timestamps RLE encoded
        expected_buffer += RLE_T(0, 0b1111);
        // 2 other byte for counts
        expected_buffer += make_array(EXPLODE16(num_events));
        for (int i = 0; (i + 4) <= num_events; i += 4) {
            // 3*4 bytes for 4 y, x and p
            expected_buffer += make_array(PACK_YS_XS_AND_PS_4(events[i], events[i + 1], events[i + 2], events[i + 3]));
        }
        // 3*4 bytes for last 4 y, x and p
        switch (num_events % 4) {
        case 1:
            expected_buffer += make_array(PACK_YS_XS_AND_PS_1(events[4 * (num_events / 4)]));
            break;
        case 2:
            expected_buffer +=
                make_array(PACK_YS_XS_AND_PS_2(events[4 * (num_events / 4)], events[4 * (num_events / 4) + 1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_YS_XS_AND_PS_3(
                events[4 * (num_events / 4)], events[4 * (num_events / 4) + 1], events[4 * (num_events / 4) + 2]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_more_than_15_events_same_t_rle_xyp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = true;
    bool pack_xs_and_ps    = false;

    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279), dy(0, 719);
    mt_rand.seed(42);

    auto num_expected_events_list = {16, 25, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_expected_events : num_expected_events_list) {
        std::vector<EventCD> expected_events;
        for (int i = 0; i < num_expected_events; ++i) {
            expected_events += EventCD{(unsigned short)(dx(mt_rand)), (unsigned short)(dy(mt_rand)), 1, 10000};
        }

        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for max count timestamps RLE encoded
        buffer += RLE_T(0, 0b1111);
        // 2 other byte for counts
        buffer += make_array(EXPLODE16(num_expected_events));
        for (int i = 0; (i + 4) <= num_expected_events; i += 4) {
            // 3*4 bytes for 4 y, x and p
            buffer += make_array(PACK_YS_XS_AND_PS_4(expected_events[i], expected_events[i + 1], expected_events[i + 2],
                                                     expected_events[i + 3]));
        }
        // 3*4 bytes for last 4 y, x and p
        switch (num_expected_events % 4) {
        case 1:
            buffer += make_array(PACK_YS_XS_AND_PS_1(expected_events[4 * (num_expected_events / 4)]));
            break;
        case 2:
            buffer += make_array(PACK_YS_XS_AND_PS_2(expected_events[4 * (num_expected_events / 4)],
                                                     expected_events[4 * (num_expected_events / 4) + 1]));
            break;
        case 3:
            buffer += make_array(PACK_YS_XS_AND_PS_3(expected_events[4 * (num_expected_events / 4)],
                                                     expected_events[4 * (num_expected_events / 4) + 1],
                                                     expected_events[4 * (num_expected_events / 4) + 2]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_more_than_31_events_same_t_rle_same_y_rle_xp_packed) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279);
    mt_rand.seed(42);

    auto num_events_list = {32, 43, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_events : num_events_list) {
        std::vector<EventCD> events;
        for (int i = 0; i < num_events; ++i) {
            events += EventCD{(unsigned short)(dx(mt_rand)), (unsigned short)(873), 1, 10000};
        }

        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for max count timestamps RLE encoded
        expected_buffer += RLE_T(0, 0b1111);
        // 2 other byte for counts
        expected_buffer += make_array(EXPLODE16(num_events));
        // 2 bytes for max count y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, 0b11111));
        // 2 other bytes for counts
        expected_buffer += make_array(EXPLODE16(num_events));
        for (int i = 0; (i + 4) <= num_events; i += 4) {
            // 2*4 bytes for 4 x and p
            expected_buffer += make_array(PACK_XS_AND_PS_4(events[i], events[i + 1], events[i + 2], events[i + 3]));
        }
        // 2*4 bytes for last 4 x and p
        switch (num_events % 4) {
        case 1:
            expected_buffer += make_array(PACK_XS_AND_PS_1(events[4 * (num_events / 4)]));
            break;
        case 2:
            expected_buffer +=
                make_array(PACK_XS_AND_PS_2(events[4 * (num_events / 4)], events[4 * (num_events / 4) + 1]));
            break;
        case 3:
            expected_buffer += make_array(PACK_XS_AND_PS_3(
                events[4 * (num_events / 4)], events[4 * (num_events / 4) + 1], events[4 * (num_events / 4) + 2]));
            break;
        }

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_more_than_31_events_same_t_rle_same_y_rle_xp_packed) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = true;

    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279);
    mt_rand.seed(42);

    auto num_expected_events_list = {32, 43, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_expected_events : num_expected_events_list) {
        std::vector<EventCD> expected_events;
        for (int i = 0; i < num_expected_events; ++i) {
            expected_events += EventCD{(unsigned short)(dx(mt_rand)), (unsigned short)(873), 1, 10000};
        }

        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for max count timestamps RLE encoded
        buffer += RLE_T(0, 0b1111);
        // 2 other byte for counts
        buffer += make_array(EXPLODE16(num_expected_events));
        // 2 bytes for max count y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, 0b11111));
        // 2 other bytes for counts
        buffer += make_array(EXPLODE16(num_expected_events));
        for (int i = 0; (i + 4) <= num_expected_events; i += 4) {
            // 2*4 bytes for 4 x and p
            buffer += make_array(PACK_XS_AND_PS_4(expected_events[i], expected_events[i + 1], expected_events[i + 2],
                                                  expected_events[i + 3]));
        }
        // 2*4 bytes for last 4 x and p
        switch (num_expected_events % 4) {
        case 1:
            buffer += make_array(PACK_XS_AND_PS_1(expected_events[4 * (num_expected_events / 4)]));
            break;
        case 2:
            buffer += make_array(PACK_XS_AND_PS_2(expected_events[4 * (num_expected_events / 4)],
                                                  expected_events[4 * (num_expected_events / 4) + 1]));
            break;
        case 3:
            buffer += make_array(PACK_XS_AND_PS_3(expected_events[4 * (num_expected_events / 4)],
                                                  expected_events[4 * (num_expected_events / 4) + 1],
                                                  expected_events[4 * (num_expected_events / 4) + 2]));
            break;
        }

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(Encoder, encode_more_than_127_events_same_t_rle_same_y_rle_x_masked) {
    Encoder encoder;
    EXPECT_EQ(4 + 8 + 65535 * (4 + 2 + 2 + 1), encoder.getCompressedSize());

    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = false;

    auto num_events_list = {128, 133, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_events : num_events_list) {
        std::vector<EventCD> events;
        for (int i = 0; i < num_events; ++i) {
            events += EventCD{(unsigned short)(i % 1278), (unsigned short)(873), 1, 10000};
        }

        std::vector<std::uint8_t> expected_buffer;
        // 4 bytes for header
        expected_buffer += make_array(HEADER(num_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        expected_buffer += make_array(ABSOLUTE_T(events[0].t));
        // 1 byte for max count timestamps RLE encoded
        expected_buffer += RLE_T(0, 0b1111);
        // 2 other byte for count
        expected_buffer += make_array(EXPLODE16(num_events));
        // 2 bytes for max count y RLE encoded
        expected_buffer += make_array(RLE_Y(events[0].y, 0b11111));
        // 2 other bytes for count
        expected_buffer += make_array(EXPLODE16(num_events));
        for (int i = 0; (i + 6) <= num_events; i += 6) {
            // 2 bytes for x masked
            expected_buffer += make_array(MASK_X(events[i].x, 0b11111));
        }
        if (num_events % 6 != 0) {
            // 2 bytes for last 6 x
            expected_buffer +=
                make_array(MASK_X(events[6 * (num_events / 6)].x, (1 << 5) - (1 << (6 - (num_events % 6)))));
        }
        // 1 bytes for max count p RLE encoded
        expected_buffer += RLE_P(events[0].p, 0b1111111);
        // 2 extra bytes for count
        expected_buffer += make_array(EXPLODE16(num_events));

        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_EQ(expected_buffer.size(), encoder(events.data(), events.data() + num_events, buffer.data()));
        for (size_t i = 0; i < expected_buffer.size(); ++i) {
            EXPECT_EQ(expected_buffer[i], buffer[i]);
        }
    }
}

TEST(Decoder, decode_more_than_127_events_same_t_rle_same_y_rle_x_masked) {
    Decoder decoder;
    bool pack_ys_xs_and_ps = false;
    bool pack_xs_and_ps    = false;

    auto num_expected_events_list = {128, 133, 254, 255, 256, 312, 1891, 18923, 65535};
    for (auto &num_expected_events : num_expected_events_list) {
        std::vector<EventCD> expected_events;
        for (int i = 0; i < num_expected_events; ++i) {
            expected_events += EventCD{(unsigned short)(i % 1278), (unsigned short)(873), 1, 10000};
        }

        std::vector<std::uint8_t> buffer;
        // 4 bytes for header
        buffer += make_array(HEADER(num_expected_events, pack_ys_xs_and_ps, pack_xs_and_ps));
        // 8 bytes for absolute t
        buffer += make_array(ABSOLUTE_T(expected_events[0].t));
        // 1 byte for max count timestamps RLE encoded
        buffer += RLE_T(0, 0b1111);
        // 2 other byte for counts
        buffer += make_array(EXPLODE16(num_expected_events));
        // 2 bytes for max count y RLE encoded
        buffer += make_array(RLE_Y(expected_events[0].y, 0b11111));
        // 2 other bytes for counts
        buffer += make_array(EXPLODE16(num_expected_events));
        for (int i = 0; (i + 6) <= num_expected_events; i += 6) {
            // 2 bytes for x masked
            buffer += make_array(MASK_X(expected_events[i].x, 0b11111));
        }
        if (num_expected_events % 6 != 0) {
            // 2 bytes for last 6 x
            buffer += make_array(MASK_X(expected_events[6 * (num_expected_events / 6)].x,
                                        (1 << 5) - (1 << (6 - (num_expected_events % 6)))));
        }
        // 1 bytes for max count p RLE encoded
        buffer += RLE_P(expected_events[0].p, 0b1111111);
        // 2 extra bytes for count
        buffer += make_array(EXPLODE16(num_expected_events));

        EXPECT_EQ(num_expected_events * sizeof(EventCD),
                  decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(buffer.data())));

        std::vector<EventCD> events(100000);
        ASSERT_EQ(expected_events.size() * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));
        for (size_t i = 0; i < expected_events.size(); ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(EncoderDecoder, encode_decode_random_events) {
    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279), dy(0, 719), dp(0, 1), dt(0, 3), dn(8817, 23887);
    mt_rand.seed(42);

    for (size_t j = 0; j < 100; ++j) {
        Encoder encoder;
        Decoder decoder;

        const size_t num_expected_events = dn(mt_rand);
        std::vector<EventCD> expected_events(num_expected_events);
        for (size_t i = 0; i < num_expected_events; ++i) {
            expected_events[i].x = dx(mt_rand);
            expected_events[i].y = dy(mt_rand);
            expected_events[i].p = dp(mt_rand);
            if (i == 0) {
                expected_events[i].t = dt(mt_rand);
            } else {
                expected_events[i].t = expected_events[i - 1].t + dt(mt_rand);
            }
        }

        std::vector<EventCD> events(100000);
        std::vector<std::uint8_t> buffer(100000);
        ASSERT_TRUE(encoder(expected_events.data(), expected_events.data() + num_expected_events, buffer.data()) > 0);
        ASSERT_EQ(num_expected_events * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));

        for (size_t i = 0; i < num_expected_events; ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}

TEST(EncoderDecoder, encode_decode_random_coherent_events) {
    std::mt19937 mt_rand; // Mersenne twister
    std::uniform_int_distribution<int> dx(0, 1279), dy(0, 719), dp(0, 1), dt(0, 3), dn(8817, 23887);
    std::discrete_distribution<int> dxx{0, 10, 5, 2, 1}, dyy{10, 1}, dpp{10, 1}, dtt{100, 1};
    mt_rand.seed(42);

    for (size_t j = 0; j < 100; ++j) {
        Encoder encoder;
        Decoder decoder;

        const size_t num_expected_events = dn(mt_rand);
        std::vector<EventCD> expected_events(num_expected_events);
        for (size_t i = 0; i < num_expected_events; ++i) {
            if (i == 0 || dx(mt_rand) < 100) {
                expected_events[i].x = dx(mt_rand);
            } else {
                expected_events[i].x = expected_events[i - 1].x + dxx(mt_rand);
                if (expected_events[i].x >= 1280) {
                    expected_events[i].x = dx(mt_rand);
                }
            }
            if (i == 0 || dy(mt_rand) < 10) {
                expected_events[i].y = dy(mt_rand);
            } else {
                expected_events[i].y = expected_events[i - 1].y + dyy(mt_rand);
                if (expected_events[i].y >= 720) {
                    expected_events[i].y = dy(mt_rand);
                }
            }
            if (i == 0) {
                expected_events[i].p = dp(mt_rand);
            } else {
                expected_events[i].p = dpp(mt_rand) == 0 ? expected_events[i - 1].p : dp(mt_rand);
            }
            if (i == 0) {
                expected_events[i].t = dt(mt_rand);
            } else {
                expected_events[i].t = expected_events[i - 1].t + dtt(mt_rand);
            }
            ASSERT_GT(1280, expected_events[i].x);
            ASSERT_LE(0, expected_events[i].x);
            ASSERT_GT(720, expected_events[i].y);
            ASSERT_LE(0, expected_events[i].y);
            ASSERT_GE(1, expected_events[i].p);
            ASSERT_LE(0, expected_events[i].p);
        }

        std::vector<EventCD> events(100000);
        std::vector<std::uint8_t> buffer(1000000);
        ASSERT_TRUE(encoder(expected_events.data(), expected_events.data() + num_expected_events, buffer.data()) > 0);
        ASSERT_EQ(num_expected_events * sizeof(EventCD),
                  decoder(buffer.data(), buffer.data() + buffer.size(), events.data()));

        for (size_t i = 0; i < num_expected_events; ++i) {
            EXPECT_EQ(expected_events[i].x, events[i].x);
            EXPECT_EQ(expected_events[i].y, events[i].y);
            EXPECT_EQ(expected_events[i].p, events[i].p);
            EXPECT_EQ(expected_events[i].t, events[i].t);
        }
    }
}
