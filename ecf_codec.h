/**********************************************************************************************************************
 * Copyright (c) Prophesee S.A. - All Rights Reserved                                                                 *
 *                                                                                                                    *
 * Subject to Prophesee Metavision Licensing Terms and Conditions ("License T&C's").                                  *
 * You may not use this file except in compliance with these License T&C's.                                           *
 * A copy of these License T&C's is located at docs.prophesee.ai/licensing and in the "LICENSE" file accompanying     *
 * this file.                                                                                                         *
 **********************************************************************************************************************/

#ifndef ECF_CODEC_H
#define ECF_CODEC_H

#include <cstddef>
#include <cstdint>
#include <vector>
#include "ecf_codec_export.h"

namespace ECF {
struct EventCD {
    unsigned short x, y;
    short p;
    long long t;
};

class Decoder {
public:
    ECF_CODEC_EXPORT Decoder();
    ECF_CODEC_EXPORT size_t operator()(const std::uint8_t *cur_raw_data, const std::uint8_t *raw_data_end,
                                       EventCD *event_ptr);
    ECF_CODEC_EXPORT size_t getDecompressedSize(const std::uint32_t *header_ptr) const;

private:
    size_t decode_ts(const std::uint8_t *ptr);
    size_t decode_ys(const std::uint8_t *ptr);
    size_t decode_ps(const std::uint8_t *ptr);
    size_t decode_xs_masked(const std::uint8_t *ptr);
    size_t decode_xs_and_ps_packed(const std::uint8_t *ptr);
    size_t decode_ys_xs_and_ps_packed(const std::uint8_t *ptr);
    void transpose(EventCD *ptr);

    size_t num_events_;
    bool ys_xs_and_ps_packed_, xs_and_ps_packed_;
    std::vector<std::uint64_t> ts_;
    std::vector<std::uint16_t> ys_, xs_;
    std::vector<std::uint8_t> ps_;
};

class Encoder {
public:
    ECF_CODEC_EXPORT Encoder();
    ECF_CODEC_EXPORT size_t operator()(const EventCD *begin_ev_ptr, const EventCD *end_ev_ptr,
                                       std::uint8_t *raw_ev_ptr);
    ECF_CODEC_EXPORT size_t getCompressedSize() const;

private:
    size_t encode(std::uint8_t *raw_ev_ptr);
    size_t encode_ts(std::uint8_t *ptr);
    size_t encode_ys(std::uint8_t *ptr);
    size_t encode_ps(std::uint8_t *ptr);
    size_t encode_xs_masked(std::uint8_t *ptr);
    size_t encode_xs_and_ps_packed(std::uint8_t *ptr);
    size_t encode_ys_xs_and_ps_packed(std::uint8_t *ptr);
    void transpose(const EventCD *begin, const EventCD *end);
    bool should_pack_ys_xs_and_ps(double ratio) const;
    bool should_pack_xs_and_ps(double ratio) const;

    size_t count_;
    std::vector<std::uint64_t> ts_;
    std::vector<std::uint16_t> ys_, xs_;
    std::vector<std::uint8_t> ps_;
};
} // namespace ECF

#endif // ECF_CODEC_H
