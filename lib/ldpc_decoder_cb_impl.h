/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */

#ifndef INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H
#define INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H

#include <ldpc_ece535a/ldpc_decoder_cb.h>

#include <gnuradio/digital/constellation.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <suitesparse/umfpack.h>

namespace gr {
  namespace ldpc_ece535a {

    using namespace boost::numeric;

    class ldpc_decoder_cb_impl : public ldpc_decoder_cb
    {
     private:
      digital::constellation_sptr d_constellation;
      ublas::matrix<int> d_H;
      unsigned int d_iterations;

      ublas::vector<int>
      decodeBitFlipping(const ublas::vector<double> &rx,
                        const ublas::matrix<int> &H,
                        const unsigned int iterations);

      int sign(double val);

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

