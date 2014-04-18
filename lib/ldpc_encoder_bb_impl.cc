/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC encoder block
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ldpc_encoder_bb_impl.h"

#include <gnuradio/io_signature.h>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace gr {
  namespace ldpc_ece535a {

    ldpc_encoder_bb::sptr
    ldpc_encoder_bb::make()
    {
      return gnuradio::get_initial_sptr
        (new ldpc_encoder_bb_impl());
    }

    /*
     * The private constructor
     */
    ldpc_encoder_bb_impl::ldpc_encoder_bb_impl()
      : gr::block("ldpc_encoder_bb",
		  gr::io_signature::make(1, 1, sizeof(unsigned char)),
		  gr::io_signature::make(1, 1, sizeof(unsigned char))),
        d_H(8, 16)
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
    ldpc_encoder_bb_impl::~ldpc_encoder_bb_impl()
    {
    }

    void
    ldpc_encoder_bb_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = 2 * noutput_items; // code rate 1/2
    }

    int
    ldpc_encoder_bb_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const unsigned char *in = (const unsigned char *) input_items[0];
      unsigned char *out = (unsigned char *) output_items[0];

      int min_input_required = d_H.size1() / 8;
      int input_available = ninput_items[0];
      int input_consumed = 0;

      while (input_available >= min_input_required) {
        ublas::vector<int> data(d_H.size1());

        // Populate data vector
        for (int i = 0; i < min_input_required; i++) {
          // Get data bits from current input byte
          for (int j = 0; j < 8; j++) {
            int b = *(in + i) & (1 << (7 - j));
            data(i * 8 + j) = b == 0 ? 0 : 1;
          }
        }

        // Encode message
        ublas::matrix<int> newH(d_H);
        const ublas::vector<int> c = makeParityCheck(data, newH);

        // Write check bytes to output
        for (int i = 0; i < min_input_required; i++) {
          // Write check bits to current output byte
          for (int j = 0; j < 8; j++) {
            if (c(i * 8 + j) == 1) {
              *out |= 1 << (7 - j);
            }
          }

          out++;
        }

        // Write data bytes to output
        for (int i = 0; i < min_input_required; i++) {
          *out = *in;

          out++;
          in++;
        }

        input_available -= min_input_required;
        input_consumed += min_input_required;
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(input_consumed);

      // Tell runtime system how many output items we produced.
      return input_consumed * 2; // code rate 1/2
    }

    ublas::vector<int>
    ldpc_encoder_bb_impl::solve(const ublas::matrix<int> &A,
                                const ublas::vector<int> &B)
    {
      const int n = A.size1();

      // Count non-zero entries
      int nz = 0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (A(i, j) != 0) {
            nz++;
          }
        }
      }

      int *Ti = new int[nz];
      int *Tj = new int[nz];
      double *Tx = new double[nz];

      int k = 0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (A(i, j) != 0) {
            Ti[k] = i;
            Tj[k] = j;
            Tx[k] = A(i, j);
            k++;
          }
        }
      }

      double *b = new double[n];

      for (int i = 0; i < n; i++) {
        b[i] = B(i);
      }

      int *Ap = new int[n + 1];
      int *Ai = new int[nz];
      double *Ax = new double[nz];

      int status = umfpack_di_triplet_to_col(n, n, nz, Ti, Tj, Tx, Ap, Ai, Ax, NULL);
      if (status != UMFPACK_OK) {
        std::cerr << "umfpack_di_triplet_to_col(): ERROR" << std::endl;
      }

      delete[] Ti;
      delete[] Tj;
      delete[] Tx;

      void *symbolic;
      status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
      if (status != UMFPACK_OK) {
        std::cerr << "umfpack_di_symbolic(): ERROR" << std::endl;
      }

      void *numeric;
      status = umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        std::cerr << "umfpack_di_numeric(): ERROR" << std::endl;
      }

      umfpack_di_free_symbolic(&symbolic);

      double *x = new double[n];
      status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        std::cerr << "umfpack_di_solve(): ERROR" << std::endl;
      }

      umfpack_di_free_numeric(&numeric);

      delete[] b;
      delete[] Ap;
      delete[] Ai;
      delete[] Ax;

      ublas::vector<int> x_vec(n);
      for (int i = 0; i < n; i++) {
        x_vec(i) = x[i];
      }

      delete[] x;

      return x_vec;
    }

    ublas::vector<int>
    ldpc_encoder_bb_impl::makeParityCheck(const ublas::vector<int> &dSource,
                                          ublas::matrix<int> &H)
    {
      // Get matrix dimensions
      const unsigned int M = H.size1();
      const unsigned int N = H.size2();

      // Set a new matrix F for LU decomposition
      ublas::matrix<int> F(H);

      // LU matrices
      ublas::matrix<int> L(ublas::zero_matrix<int>(M, N - M));
      ublas::matrix<int> U(ublas::zero_matrix<int>(M, N - M));

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

      // Find B.dsource
      const ublas::vector<int> z(mod2(ublas::prod(ublas::subrange(H, 0, M, N - M, N), dSource)));

      // Parity check vector found by solving sparse LU
      const ublas::vector<int> x1(solve(L, z));
      const ublas::vector<int> x2(solve(U, x1));
      const ublas::vector<int> c(mod2(x2));

      return c;
    }

    ublas::vector<int>
    ldpc_encoder_bb_impl::mod2(const ublas::vector<int> &u)
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
