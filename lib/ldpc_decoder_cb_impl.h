/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */

#ifndef INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H
#define INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H

#include <ldpc_ece535a/ldpc_decoder_cb.h>

namespace gr {
  namespace ldpc_ece535a {

    class ldpc_decoder_cb_impl : public ldpc_decoder_cb
    {
     private:
      // Nothing to declare in this block.

     public:
      ldpc_decoder_cb_impl();
      ~ldpc_decoder_cb_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H */

