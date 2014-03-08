/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "ldpc_decoder_bb_impl.h"

namespace gr {
  namespace ldpc_ece535a {

    ldpc_decoder_bb::sptr
    ldpc_decoder_bb::make()
    {
      return gnuradio::get_initial_sptr
        (new ldpc_decoder_bb_impl());
    }

    /*
     * The private constructor
     */
    ldpc_decoder_bb_impl::ldpc_decoder_bb_impl()
      : gr::block("ldpc_decoder_bb",
		  gr::io_signature::make(1, 1, sizeof(unsigned char)),
		  gr::io_signature::make(1, 1, sizeof(unsigned char)))
    {}

    /*
     * Our virtual destructor.
     */
    ldpc_decoder_bb_impl::~ldpc_decoder_bb_impl()
    {
    }

    void
    ldpc_decoder_bb_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    int
    ldpc_decoder_bb_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const unsigned char *in = (const unsigned char *) input_items[0];
      unsigned char *out = (unsigned char *) output_items[0];

      // Do <+signal processing+>
      // Just copy input to output for now...
      std::memcpy(out, in, noutput_items * sizeof(unsigned char));

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */
