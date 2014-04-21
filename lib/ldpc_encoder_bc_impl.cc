/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC encoder block
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ldpc_encoder_bc_impl.h"

#include <gnuradio/io_signature.h>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <lapacke/lapacke.h>

namespace gr {
  namespace ldpc_ece535a {

    ldpc_encoder_bc::sptr
    ldpc_encoder_bc::make()
    {
      return gnuradio::get_initial_sptr
        (new ldpc_encoder_bc_impl());
    }

    /*
     * The private constructor
     */
    ldpc_encoder_bc_impl::ldpc_encoder_bc_impl()
      : gr::block("ldpc_encoder_bc",
		  gr::io_signature::make(1, 1, sizeof(unsigned char)),
		  gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_M(8), d_N(16), d_H(d_M, d_N),
        d_L(ublas::zero_matrix<int>(d_M, d_N - d_M)),
        d_U(ublas::zero_matrix<int>(d_M, d_N - d_M))
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

      for (unsigned int i = 0; i < d_M; i++) {
        for (unsigned int j = 0; j < d_N; j++) {
          d_H(i, j) = h_data[i * d_N + j];
        }
      }

      reorderHMatrix(d_H, d_L, d_U);
    }

    /*
     * Our virtual destructor.
     */
    ldpc_encoder_bc_impl::~ldpc_encoder_bc_impl()
    {
    }

    void
    ldpc_encoder_bc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      // code rate 1/2
      ninput_items_required[0] = std::ceil(noutput_items / 16.0);
    }

    int
    ldpc_encoder_bc_impl::general_work (int noutput_items,
					gr_vector_int &ninput_items,
					gr_vector_const_void_star &input_items,
					gr_vector_void_star &output_items)
    {
      const unsigned char *in = (const unsigned char *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      int min_input_required = d_M / 8;
      int min_output_required = d_N;

      int input_consumed = 0;
      int output_produced = 0;

      while ((noutput_items - output_produced) >= min_output_required
             && (ninput_items[0] - input_consumed) >= min_input_required) {
        ublas::vector<int> data(d_M);

        // Populate data vector
        for (int i = 0; i < min_input_required; i++) {
          // Get data bits from current input byte
          for (int j = 0; j < 8; j++) {
            int b = *(in + i) & (1 << (7 - j));
            data(i * 8 + j) = b == 0 ? 0 : 1;
          }

          in++;
        }

        // Encode message
        const ublas::vector<int> c = makeParityCheck(data, d_H, d_L, d_U);

        // Write check to output
        for (unsigned int i = 0; i < c.size(); i++) {
          *out = c(i) == 1 ? 1 : -1;
          out++;
        }

        // Write data to output
        for (unsigned int i = 0; i < data.size(); i++) {
          *out = data(i) == 1 ? 1 : -1;
          out++;
        }

        input_consumed += min_input_required;
        output_produced += d_N;
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each(input_consumed);

      // Tell runtime system how many output items we produced.
      return output_produced; // code rate 1/2
    }

    ublas::vector<int>
    ldpc_encoder_bc_impl::solve(const ublas::matrix<int> &A,
                                const ublas::vector<int> &B)
    {
      const int n = A.size1();
      const int nrhs = 1;
      const int lda = n;
      const int ldb = nrhs;

      int *ipiv = new int[n];

      double *a = new double[lda * n];
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          a[i * n + j] = A(i, j);
        }
      }

      double *b = new double[ldb * n];
      for (int i = 0; i < n; i++) {
        b[i] = B(i);
      }

      ublas::vector<int> x_vec(n);

      const int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
      if (info > 0) {
        std::cerr << "The diagonal element of the triangular factor of A, "
                  << "U(" << info << "," << info << ") is zero, so that A is singular; "
                  << "the solution could not be computed." << std::endl;
      } else if (info < 0) {
        std::cerr << "Illegal value for argument " << -info << "." << std::endl;
      } else {
        for (int i = 0; i < n; i++) {
          x_vec(i) = b[i];
        }
      }

      delete[] a;
      delete[] b;
      delete[] ipiv;

      return x_vec;
    }

    void ldpc_encoder_bc_impl::reorderHMatrix(ublas::matrix<int> &H, ublas::matrix<int> &L, ublas::matrix<int> &U) {
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
    ldpc_encoder_bc_impl::makeParityCheck(const ublas::vector<int> &dSource,
                                          const ublas::matrix<int> &H,
                                          const ublas::matrix<int> &L,
                                          const ublas::matrix<int> &U)
    {
      // Get matrix dimensions
      const unsigned int M = H.size1();
      const unsigned int N = H.size2();

      // Find B.dsource
      const ublas::vector<int> z(mod2(ublas::prod(ublas::subrange(H, 0, M, N - M, N), dSource)));

      // Parity check vector found by solving sparse LU
      const ublas::vector<int> x1(solve(L, z));
      const ublas::vector<int> x2(solve(U, x1));
      const ublas::vector<int> c(mod2(x2));

      return c;
    }

    ublas::vector<int>
    ldpc_encoder_bc_impl::mod2(const ublas::vector<int> &u)
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
