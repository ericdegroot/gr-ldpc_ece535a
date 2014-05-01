/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * Image sink block
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "image_sink_impl.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

namespace gr {
  namespace ldpc_ece535a {

    image_sink::sptr
    image_sink::make()
    {
      return gnuradio::get_initial_sptr
        (new image_sink_impl());
    }

    /*
     * The private constructor
     */
    image_sink_impl::image_sink_impl()
      : gr::sync_block("image_sink",
                       gr::io_signature::make(1, 1, sizeof(unsigned char)),
                       gr::io_signature::make(0, 0, 0)),
        d_fileSize(0)

    {}

    /*
     * Our virtual destructor.
     */
    image_sink_impl::~image_sink_impl()
    {
    }

    int
    image_sink_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];

        for (int i = 0; i < noutput_items; i++) {
          // Look for BMP header
          if (i < noutput_items - 18 && in[i] == 'B' && in[i + 1] == 'M' && in[i + 6] == 0
              && in[i + 7] == 0 && in[i + 8] == 0 && in[i + 9] == 0
              && (in[i + 14] == 12 || in[i + 14] == 40 || in[i + 14] == 52
                  || in[i + 14] == 56 || in[i + 14] == 64 || in[i + 14] == 108
                  || in[i + 14] == 124)) {
            if (d_fileSize > 0 && d_buffer.length() >= d_fileSize) {
              std::ofstream file("result.bmp", std::ios::out | std::ios::binary);
              if (file.is_open()) {
                file.write(d_buffer.data(), d_fileSize);
                file.close();

                std::cout << "File written" << std::endl;

                std::system("/usr/bin/display result.bmp &");
              }
            }

            d_buffer.erase();

            d_fileSize = (in[i + 5] << 24) | (in[i + 4] << 16) | (in[i + 3] << 8) | in[i + 2];

            std::cout << "BMP Header Found: fileSize=" << d_fileSize << std::endl;
          }

          d_buffer.push_back(in[i]);
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */

