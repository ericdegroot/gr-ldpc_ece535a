/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */

#ifndef INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H
#define INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_IMPL_H

#include <ldpc_ece535a/ldpc_decoder_cb.h>

#include <boost/random.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace gr {
  namespace ldpc_ece535a {

    using namespace boost::numeric;

    class ldpc_decoder_cb_impl : public ldpc_decoder_cb
    {
     private:
      int d_method;
      int d_state;
      unsigned int d_M;
      unsigned int d_N;
      ublas::matrix<int> d_H;
      unsigned int d_iterations;
      unsigned int d_errors;

      boost::random::mt19937 d_gen;
      boost::random::uniform_int_distribution<> d_uDist;

      int checkFrame(const ublas::vector<int> &u,
                     const int threshold);
      void reorderHMatrix(ublas::matrix<int> &H,
                          ublas::matrix<int> &L,
                          ublas::matrix<int> &U);
      ublas::vector<int> decodeLogDomainSimple(const ublas::vector<double> &rx,
                                               const ublas::matrix<int> &H,
                                               const unsigned int iterations);
      ublas::vector<int> decodeBitFlipping(const ublas::vector<double> &rx,
                                           const ublas::matrix<int> &H,
                                           const unsigned int iterations);
      ublas::vector<int> decodeSumProductSoft(const ublas::vector<double> &rx,
                                              const ublas::matrix<int> &H,
                                              const unsigned int iterations);
      ublas::vector<int> decodeHard(const ublas::vector<double> &rx);

      int sign(double val);
      ublas::vector<int> mod2(const ublas::vector<int> &u);

     public:
      ldpc_decoder_cb_impl(const int method);
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

