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

#define STATE_OUT_OF_SYNC 0
#define STATE_IN_SYNC 1
#define STATE_IN_SYNC_INVERTED 2

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
        d_state(STATE_OUT_OF_SYNC), d_M(32), d_N(64), d_H(d_M, d_N), d_iterations(5),
        d_uDist(1, d_M)
    {
      d_gen.seed(static_cast<unsigned int>(std::time(0)));

      /*
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
      */

      // M = 32
      // N = 64
      // makeLdpc(M, N, 1, 1, 3)
      const int h_data[] = {
        0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,
        1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,1,
        0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
        0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,
        0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,
        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,
        0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
        0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
        0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
        0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
      };
 
      for (unsigned int i = 0; i < d_M; i++) {
        for (unsigned int j = 0; j < d_N; j++) {
          d_H(i, j) = h_data[i * d_N + j];
        }
      }

      ublas::matrix<int> L(ublas::zero_matrix<int>(d_M, d_N - d_M));
      ublas::matrix<int> U(ublas::zero_matrix<int>(d_M, d_N - d_M));
      reorderHMatrix(d_H, L, U);
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
      ninput_items_required[0] = noutput_items * d_N;
    }

    int
    ldpc_decoder_cb_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      unsigned char *out = (unsigned char *) output_items[0];

      const int min_output_required = d_M / 8;
      const int frame_error_threshold = d_M; //0.1 * d_M;

      int input_consumed = 0;
      int output_produced = 0;

      while ((ninput_items[0] - input_consumed) >= d_N
             && (noutput_items - output_produced) >= min_output_required) {
        ublas::vector<double> tx(d_N);
        for (int i = 0; i < d_N; i++) {
          tx(i) = (in + i)->real()
            * (d_state == STATE_IN_SYNC_INVERTED ? -1 : 1);
        }

        //ublas::vector<int> vhat = decodeBitFlipping(tx, d_H, d_iterations);
        ublas::vector<int> vhat = decodeLogDomainSimple(tx, d_H, d_iterations);

        int sNotZero = checkFrame(vhat, frame_error_threshold);
        //if (sNotZero > frame_error_threshold) {
        if (sNotZero > 0) {
          // Check if phase lock has flipped
          vhat = decodeLogDomainSimple(-tx, d_H, d_iterations);
          if ((sNotZero = checkFrame(vhat, 0)) == 0) {
            if (d_state == STATE_IN_SYNC_INVERTED) {
              std::cout << "PHASE RESTORED" << std::endl;
              d_state = STATE_IN_SYNC;
            } else {
              std::cout << "PHASE INVERTED" << std::endl;
              d_state = STATE_IN_SYNC_INVERTED;
            }
          } else {
            if (d_state == STATE_IN_SYNC) {
              std::cout << "STATE_OUT_OF_SYNC" << std::endl;
              d_state = STATE_OUT_OF_SYNC;
            }

            // Skip ahead a random number of samples in search of the
            // next frame
            int skip = std::min(ninput_items[0] - input_consumed, d_uDist(d_gen));
            in += skip;
            input_consumed += skip;
          }
        } 

        if (sNotZero == 0) {
          if (d_state == STATE_OUT_OF_SYNC) {
            std::cout << "STATE_IN_SYNC" << std::endl;
            d_state = STATE_IN_SYNC;
          }

          // Send data bits to output byte stream
          for (int i = 0; i < min_output_required; i++) {
            *out = 0;

            for (int j = 0; j < 8; j++) {
              if (vhat(d_M + i * 8 + j) == 1) {
                *out |= 1 << (7 - j);
              }
            }

            out++;
          }

          in += d_N;

          input_consumed += d_N;
          output_produced += min_output_required;
        }
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(input_consumed);

      // Tell runtime system how many output items we produced.
      return output_produced;
    }

    int
    ldpc_decoder_cb_impl::checkFrame(const ublas::vector<int> &u,
                                     const int threshold)
    {
      int sNotZero = 0;
      for (unsigned int k = 0; k < d_H.size1(); k++) {
        int s = ublas::inner_prod(u, ublas::row(d_H, k)) % 2;
        if (s != 0) {
          sNotZero++;
        }

        if (sNotZero > threshold) {
          break;
        }
      }

      return sNotZero;
    }

    void
    ldpc_decoder_cb_impl::reorderHMatrix(ublas::matrix<int> &H,
                                         ublas::matrix<int> &L,
                                         ublas::matrix<int> &U)
    {
      // Get matrix dimensions
      const unsigned int M = H.size1();
      const unsigned int N = H.size2();

      // Set a new matrix F for LU decomposition
      ublas::matrix<int> F(H);

      // Re-order the M x (N - M) submatrix
      for (unsigned int i = 0; i < M; i++) {
        int chosenCol = 0;

        // Create diagonally structured matrix using 'First' strategy
        for (unsigned int j = i; j < N; j++) {
          if (F(i, j) != 0) {
            chosenCol = j;
            break;
          }
        }

        //std::cout << "chosenCol=" << chosenCol << std::endl;
        //printMatrix<int>(H);

        // Re-ordering columns of F
        const ublas::vector<int> tmp1(ublas::column(F, i));
        ublas::column(F, i) = ublas::column(F, chosenCol);
        ublas::column(F, chosenCol) = tmp1;

        // Re-ordering columns of H
        const ublas::vector<int> tmp2(ublas::column(H, i));
        ublas::column(H, i) = ublas::column(H, chosenCol);
        ublas::column(H, chosenCol) = tmp2;

        // Fill the LU matrices column by column
        ublas::subrange(L, i, M, i, i + 1) = ublas::subrange(F, i, M, i, i + 1);
        ublas::subrange(U, 0, i + 1, i, i + 1) = ublas::subrange(F, 0, i + 1, i, i + 1);

        // There will be no rows operation at the last row
        if (i < M - 1) {
          // Find the later rows with non-zero elements in column i
          for (unsigned int k = i + 1; k < M; k++) {
            if (F(k, i) != 0) {
              // Add current row to the later rows which have a 1 in column i
              ublas::row(F, k) = mod2(ublas::row(F, k) + ublas::row(F, i));
            }
          }
        }
      }  
    }

    ublas::vector<int>
    ldpc_decoder_cb_impl::decodeLogDomainSimple(const ublas::vector<double> &rx,
                                                const ublas::matrix<int> &H,
                                                const unsigned int iterations) {
      // Get matrix dimensions
      const unsigned int M = H.size1();
      const unsigned int N = H.size2();

      // Prior hard-decision
      ublas::vector<double> Lci(rx.size());
      for (unsigned int i = 0; i < rx.size(); i++) {
        Lci(i) = -rx(i);
      }

      // Initialization
      ublas::matrix<double> Lrji(ublas::zero_matrix<double>(M, N));
      ublas::matrix<double> Pibetaij(ublas::zero_matrix<double>(M, N));

      // Associate the L(ci) matrix with non-zero elements of H
      ublas::matrix<double> Lqij(M, N);
      for (unsigned int i = 0; i < M; i++) {
        ublas::row(Lqij, i) = ublas::element_prod(ublas::row(H, i), Lci);
      }

      ublas::vector<int> vHat(N);

      // Iteration
      for (unsigned int n = 0; n < iterations; n++) {
        //std::cout << "Iteration : " << n << std::endl;

        // Get the sign and magnitude of L(qij)
        ublas::matrix<int> alphaij(M, N);
        ublas::matrix<double> betaij(M, N);
        for (unsigned int i = 0; i < M; i++) {
          for (unsigned int j = 0; j < N; j++) {
            alphaij(i, j) = sign(Lqij(i, j));
            betaij(i, j) = std::abs(Lqij(i, j));
          }
        }

        // Horizontal step
        for (unsigned int i = 0; i < M; i++) {
          int prodOfalphaij = 1;
          for (unsigned int j = 0; j < N; j++) {
            if (H(i, j) != 0) {
              prodOfalphaij *= alphaij(i, j);
            }
          }

          // Get the minimum of betaij
          for (unsigned int j = 0; j < N; j++) {
            if (H(i, j) != 0) {
              // Minimum of betaij
              double minOfBetaij = std::numeric_limits<double>::max();
              for (unsigned int k = 0; k < N; k++) {
                if (j != k && H(i, k) != 0) {
                  if (betaij(i, k) < minOfBetaij) {
                    minOfBetaij = betaij(i, k);
                  }
                }
              }

              // Multiplication alphaij
              // Update L(rji)
              Lrji(i, j) = prodOfalphaij * alphaij(i, j) * minOfBetaij;
            }
          }
        }

        // Vertical step
        for (unsigned int j = 0; j < N; j++) {
          double sumOfLrji = 0.0;
          for (unsigned int i = 0; i < M; i++) {
            if (H(i, j) != 0) {
              sumOfLrji += Lrji(i, j);
            }
          }

          for (unsigned int i = 0; i < M; i++) {
            if (H(i, j) != 0) {
              // Update L(qij) by summation of L(rij)
              Lqij(i, j) = Lci(j) + sumOfLrji - Lrji(i, j);
            }
          }

          // Get L(Qij)
          double LQi = Lci(j) + sumOfLrji;

          // Decode L(Qi)
          if (LQi < 0) {
            vHat(j) = 1;
          } else {
            vHat(j) = 0;
          }
        }

        // Check result for early exit
        if (n + 1 < iterations && checkFrame(vHat, 0) == 0) {
          break;
        }
      }

      return vHat;
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

    ublas::vector<int>
    ldpc_decoder_cb_impl::mod2(const ublas::vector<int> &u)
    {
      ublas::vector<int> v(u.size());

      for (unsigned int i = 0; i < u.size(); i++) {
        v(i) = u(i) % 2;

        // modulus result can be negative if u(i) negative, depending on implementation
        if (v(i) < 0) {
          v(i) += 2;
        }
      }

      return v;
    }

  } /* namespace ldpc_ece535a */
} /* namespace gr */
