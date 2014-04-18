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

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

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
        d_constellation(digital::constellation_bpsk::make()),
        d_H(8, 16), d_iterations(5)
    {
      // M = 8
      // N = 16
      // makeLdpc(M, N, 1, 1, 3)
      const int h_data[] = {
        1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,
        0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,
        0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,
        0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,
        0,0,0,0,0,1,0,1,0,1,0,0,0,0,1,0,
        0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,
        0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,
        0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1
      };
 
      for (unsigned int i = 0; i < d_H.size1(); i++) {
        for (unsigned int j = 0; j < d_H.size2(); j++) {
          d_H(i, j) = h_data[i * d_H.size2() + j];
        }
      }
    }

    /*
     * Our virtual destructor.
     */
    ldpc_decoder_cb_impl::~ldpc_decoder_cb_impl()
    {
    }

    void
    ldpc_decoder_cb_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items * d_H.size1();
    }

    int
    ldpc_decoder_cb_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      unsigned char *out = (unsigned char *) output_items[0];

      int input_available = ninput_items[0];
      int output_produced = 0;

      /*
      for (int i = 0; i < noutput_items; i++) {
        out[i] = d_constellation->decision_maker(&in[i]);
      }
      */

      while (input_available >= d_H.size2()) {
        ublas::vector<double> tx(d_H.size2());
        for (int i = 0; i < d_H.size2(); i++) {
          tx(i) = (*in).real();
          in++;
        }

        /*
        std::cout << "tx=[";
        for (int i = 0; i < tx.size(); i++)
          std::cout << tx(i) << ",";
        std::cout << "]" << std::endl;
        */

        input_available -= d_H.size2();

        ublas::vector<int> vhat = decodeBitFlipping(tx, d_H, d_iterations);

        /*
        std::cout << "vhat=[";
        for (int i = 0; i < vhat.size(); i++)
          std::cout << vhat(i) << ",";
        std::cout << "]" << std::endl;
        */

        int numBytes = vhat.size() / 8;
        for (int i = 0; i < numBytes; i++) {
          for (int j = 0; j < 8; j++) {
            if (vhat(i * 8 + j) == 1) {
              *out |= 1 << (7 - j);
            }
          }

          out++;
          output_produced++;
        }
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(ninput_items[0] - input_available);

      // Tell runtime system how many output items we produced.
      return output_produced;
    }

    ublas::vector<int>
    ldpc_decoder_cb_impl::decodeBitFlipping(const ublas::vector<double> &rx,
                                            const ublas::matrix<int> &H,
                                            const unsigned int iterations)
    {
      // Get matrix dimensions
      const unsigned int M = H.size1();
      const unsigned int N = H.size2();

      // Prior hard-decision
      ublas::vector<int> ci(rx.size());
      for (unsigned int i = 0; i < rx.size(); i++) {
        ci(i) = 0.5 * (sign(rx(i)) + 1);
      }

      //std::cout << "ci=" << std::endl;
      //printVector<int>(ci);

      // Initialization
      ublas::matrix<int> rji(ublas::zero_matrix<int>(M, N));

      // Associate the ci matrix with non-zero elements of H
      ublas::matrix<int> qij(M, N);
      for (unsigned int i = 0; i < M; i++) {
        ublas::row(qij, i) = ublas::element_prod(ublas::row(H, i), ci);
      }

      //std::cout << "qij=" << std::endl;
      //printMatrix<int>(qij);

      ublas::vector<int> vHat(N);

      // Iteration
      for (unsigned int n = 0; n < iterations; n++) {
        // std::cout << "Iteration : " << n << std::endl;

        // Horizontal step
        for (unsigned int i = 0; i < M; i++) {
          int qijSum = 0;
          for (unsigned int j = 0; j < N; j++) {
            if (H(i, j) != 0) {
              qijSum += qij(i, j);
            }
          }

          for (unsigned int j = 0; j < N; j++) {
            if (H(i, j) != 0) {
              rji(i, j) = (qijSum + qij(i, j)) % 2;
            }
          }
        }

        // Vertical step
        for (unsigned int j = 0; j < N; j++) {
          int rjiNumberOfOnes = 0;
          for (unsigned int i = 0; i < M; i++) {
            if (rji(i, j) != 0) {
              rjiNumberOfOnes++;
            }
          }

          int hNumberOfOnes = 0;
          for (unsigned int i = 0; i < M; i++) {
            if (H(i, j) != 0) {
              hNumberOfOnes++;
            }
          }
      
          for (unsigned int i = 0; i < M; i++) {
            if (H(i, j) != 0) {
              // Update qij, set '1' for majority of 1s else '0', excluding i
              if (rjiNumberOfOnes + ci(j) >= hNumberOfOnes - rjiNumberOfOnes + rji(i, j)) {
                qij(i, j) = 1;
              } else {
                qij(i, j) = 0;
              }
            }
          }

          // Bit decoding
          if (rjiNumberOfOnes + ci(j) >= hNumberOfOnes - rjiNumberOfOnes) {
            vHat(j) = 1;
          } else {
            vHat(j) = 0;
          }
        }
      }

      return vHat;
    }

    int
    ldpc_decoder_cb_impl::sign(double val)
    {
      return (val > 0) - (val < 0);
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */
