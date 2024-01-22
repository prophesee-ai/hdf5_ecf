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

#ifndef ECF_H5FILTER_H
#define ECF_H5FILTER_H

#include <cstdlib>
#include <hdf5.h>
#include "ecf_codec.h"

#define H5Z_FILTER_ECF 0x8ECF

#define H5Z_ECF_PUSH_ERR(func, minor, str) H5Epush1(__FILE__, func, __LINE__, H5E_PLINE, minor, str)

static size_t H5Z_filter_ecf(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                             size_t *buf_size, void **buf) {
    if (NULL == buf) {
        return 0;
    }

    void *outbuf = NULL;
    void *inbuf  = NULL;
    inbuf        = *buf;

    size_t ret_value;
    if (flags & H5Z_FLAG_REVERSE) {
        ECF::Decoder decoder;

        size_t outbuf_size = decoder.getDecompressedSize(reinterpret_cast<const std::uint32_t *>(inbuf));
        if (NULL == (outbuf = malloc(outbuf_size))) {
            return 0;
        }

        const std::uint8_t *begin_in_ptr = reinterpret_cast<const std::uint8_t *>(inbuf),
                           *end_in_ptr   = reinterpret_cast<const std::uint8_t *>(inbuf) + nbytes;
        ECF::EventCD *out_ptr            = reinterpret_cast<ECF::EventCD *>(outbuf);
        outbuf_size                      = decoder(begin_in_ptr, end_in_ptr, out_ptr);

        free(*buf);
        *buf      = outbuf;
        outbuf    = NULL;
        ret_value = outbuf_size;
    } else {
        const size_t event_size = sizeof(ECF::EventCD);
        const size_t num_events = nbytes / event_size;
        ECF::Encoder encoder;

        size_t outbuf_size = encoder.getCompressedSize();
        if (NULL == (outbuf = malloc(outbuf_size))) {
            return 0;
        }

        const ECF::EventCD *begin_in_ptr = reinterpret_cast<const ECF::EventCD *>(inbuf),
                           *end_in_ptr   = reinterpret_cast<const ECF::EventCD *>(inbuf) + num_events;
        std::uint8_t *out_ptr            = reinterpret_cast<std::uint8_t *>(outbuf);
        outbuf_size                      = encoder(begin_in_ptr, end_in_ptr, out_ptr);

        free(*buf);
        *buf      = outbuf;
        *buf_size = outbuf_size;
        outbuf    = NULL;
        ret_value = outbuf_size;
    }

    if (outbuf != NULL) {
        free(outbuf);
    }
    return ret_value;
}

const H5Z_class2_t H5Z_ECF[1] = {{H5Z_CLASS_T_VERS, (H5Z_filter_t)(H5Z_FILTER_ECF), 1, 1,
                                  "HDF ECF filter; see http://www.hdfgroup.org/services/contributions.html", NULL, NULL,
                                  (H5Z_func_t)(H5Z_filter_ecf)}};

[[maybe_unused]] static int ecf_register_h5filter() {
    int retval = H5Zregister(H5Z_ECF);
    if (retval < 0) {
        H5Z_ECF_PUSH_ERR("Register ECF", H5E_CANTREGISTER, "Can't register ECF filter");
    }
    return retval;
}

#undef H5Z_ECF_PUSH_ERR

#endif /* ECF_H5FILTER_H */
