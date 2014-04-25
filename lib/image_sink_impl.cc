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
              gr::io_signature::make(0, 0, 0))
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
          /*
          if (i < noutput_items - 14 && in[i] == 0x42 && in[i + 1] == 0x4D && in[i + 2] == 0x36 &&
              in[i + 3] == 0x4C && in[i + 4] == 0x02 && in[i + 5] == 0x00 && in[i + 6] == 0x00 &&
              in[i + 7] == 0x00 && in[i + 8] == 0x00 && in[i + 9] == 0x00 && in[i + 10] == 0x36 &&
              in[i + 11] == 0x00 && in[i + 12] == 0x00 && in[i + 13] == 0x00) {
          */
          if (i < noutput_items - 4 && in[i] == 0x42 && in[i + 1] == 0x4D && in[i + 2] == 0x46 &&
              in[i + 3] == 0x4B) {
            std::cout << "Header Found" << std::endl;

            if (d_buffer.length() > 0) {
              std::ofstream file("result.bmp", std::ios::out | std::ios::binary);
              if (file.is_open()) {
                file.write(d_buffer.data(), d_buffer.length());
                file.close();

                std::cout << "File written" << std::endl;
              }

              d_buffer.erase();
            }
          }

          d_buffer.push_back(in[i]);
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */

