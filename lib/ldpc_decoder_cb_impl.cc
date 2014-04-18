/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ldpc_decoder_cb_impl.h"

#include <gnuradio/io_signature.h>

namespace gr {
  namespace ldpc_ece535a {

    ldpc_decoder_cb::sptr
    ldpc_decoder_cb::make()
    {
      return gnuradio::get_initial_sptr
        (new ldpc_decoder_cb_impl());
    }

    /*
     * The private constructor
     */
    ldpc_decoder_cb_impl::ldpc_decoder_cb_impl()
      : gr::block("ldpc_decoder_cb",
		  gr::io_signature::make(1, 1, sizeof(gr_complex)),
                  gr::io_signature::make(1, 1, sizeof(unsigned char))),
        d_constellation(digital::constellation_bpsk::make())
    {}

    /*
     * Our virtual destructor.
     */
    ldpc_decoder_cb_impl::~ldpc_decoder_cb_impl()
    {
    }

    void
    ldpc_decoder_cb_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    int
    ldpc_decoder_cb_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      unsigned char *out = (unsigned char *) output_items[0];

      for (int i = 0; i < noutput_items; i++) {
        out[i] = d_constellation->decision_maker(&in[i]);
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */
