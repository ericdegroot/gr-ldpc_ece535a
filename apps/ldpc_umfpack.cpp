/*
ECE 535A LDPC Project Testing

required packages: boost-devel suitesparse-devel

g++ -Wall -o ldpc_umfpack ldpc_umfpack.cpp -lumfpack -lamd
./ldpc_umfpack

 */
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random.hpp>

#include <suitesparse/umfpack.h>

using namespace boost::numeric;

template <typename T>
void initMatrix(ublas::matrix<T> &A, const T data[]) {
  for (unsigned int i = 0; i < A.size1(); i++) {
    for (unsigned int j = 0; j < A.size2(); j++) {
      A(i, j) = data[i * A.size2() + j];
    }
  }
}

template <typename T>
void initVector(ublas::vector<T> &u, const T data[]) {
  for (unsigned int i = 0; i < u.size(); i++) {
    u(i) = data[i];
  }
}

template <typename T>
void printMatrix(const ublas::matrix<T> &A) {
  for (unsigned int i = 0; i < A.size1(); i++) {
    for (unsigned int j = 0; j < A.size2(); j++) {
      std::cout << A(i, j) << "\t";
    }

    std::cout << std::endl;
  } 
}

template <typename T>
void printVector(const ublas::vector<T> &u) {
  for (unsigned int i = 0; i < u.size(); i++) {
    std::cout << u(i) << "\t";
  }
  
  std::cout << std::endl;
}

ublas::vector<int> mod2(const ublas::vector<int> &u) {
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

int sign(double val) {
  return (val > 0) - (val < 0);
}

ublas::vector<int> solve(const ublas::matrix<int> &A, const ublas::vector<int> &B) {
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

ublas::vector<int> makeParityCheck(const ublas::vector<int> &dSource, ublas::matrix<int> &H) {
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

  // Find B.dsource
  const ublas::vector<int> z(mod2(ublas::prod(ublas::subrange(H, 0, M, N - M, N), dSource)));

  //std::cout << "z=" << std::endl;
  //printVector<int>(z);

  //std::cout << "L=" << std::endl;
  //printMatrix<int>(L);

  //std::cout << "U=" << std::endl;
  //printMatrix<int>(U);

  // Parity check vector found by solving sparse LU
  const ublas::vector<int> x1(solve(L, z));
  //std::cout << "x1=" << std::endl;
  //printVector<int>(x1);

  const ublas::vector<int> x2(solve(U, x1));
  //std::cout << "x2=" << std::endl;
  //printVector<int>(x2);

  const ublas::vector<int> c(mod2(x2));

  return c;
}

ublas::vector<int> decodeBitFlipping(const ublas::vector<double> &rx, const ublas::matrix<int> &H, const unsigned int iterations) {
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
    std::cout << "Iteration : " << n << std::endl;

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

double biterr(const ublas::vector<int> &u, const ublas::vector<int> &v) {
  int numErr = 0;
  for (unsigned int i = 0; i < u.size(); i++) {
    if (u(i) != v(i)) {
      numErr++;
    }
  }

  return static_cast<double>(numErr) / u.size();
}

const int M1 = 5;
const int N1 = 10;
const int hData1[] = {
  1, 0, 1, 0, 1, 0, 0, 0, 1, 0,
  0, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
  0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 1, 1, 1, 1, 0
};

const int M2 = 50;
const int N2 = 100;
const int hData2[] = {
  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,
  0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,
  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0
};

const int M3 = 8;
const int N3 = 16;
const int hData3[] = {
  1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,
  0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,
  0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,
  0,0,0,0,0,1,0,1,0,1,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,
  0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,
  0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1
};

int main() {
  boost::random::mt19937 gen;
  gen.seed(static_cast<unsigned int>(std::time(0)));

  const unsigned int M = M3;
  const unsigned int N = N3;

  const unsigned int iterations = 5;
  const unsigned int frames = 30;

  const double EbN0[4] = { 0, 0.5, 1.0, 1.5 };
  double ber1[4];

  ublas::matrix<int> H(M, N);

  initMatrix<int>(H, hData3);

  boost::random::uniform_int_distribution<> uniformDistribution(0, 1);
  boost::random::normal_distribution<> normalDistribution;

  for (unsigned int i = 0; i < 4; i++) {
    ber1[i] = 0;

    // Make random data (0/1)
    ublas::matrix<int> dSource(M, frames);
    for (unsigned int j = 0; j < frames; j++) {
      for (unsigned int k = 0; k < M; k++) {
        dSource(k, j) = uniformDistribution(gen);
      }
    }

    //std::cout << "dSource=" << std::endl;
    //printMatrix<int>(dSource);

    for (unsigned int j = 0; j < frames; j++) {
      std::cout << "Frame : " << j << std::endl;

      // Encoding message
      ublas::matrix<int> newH(H);
      const ublas::vector<int> c = makeParityCheck(ublas::column(dSource, j), newH);
      
      //std::cout << "newH=" << std::endl;
      //printMatrix<int>(newH);

      //std::cout << "c=" << std::endl;
      //printVector<int>(c);

      ublas::vector<int> u(c.size() + M);
      ublas::subrange(u, 0, c.size()) = c;
      ublas::subrange(u, c.size(), c.size() + M) = ublas::column(dSource, j);

      //std::cout << "u=" << std::endl;
      //printVector<int>(u);

      // BPSK modulation
      ublas::vector<int> bpskMod(u.size());
      for (unsigned int k = 0; k < u.size(); k++) {
        bpskMod(k) = 2 * u(k) - 1;
      }

      //std::cout << "bpskMod=" << std::endl;
      //printVector<int>(bpskMod);

      // Additional white gaussian noise
      double N0 = 1 / (std::exp(EbN0[i] * std::log(10) / 10));

      //std::cout << "N0=" << N0 << std::endl;

      ublas::vector<double> tx(bpskMod);
      for (unsigned int k = 0; k < tx.size(); k++) {
        tx(k) += std::sqrt(N0) * normalDistribution(gen);
      }

      std::cout << "tx=" << std::endl;
      printVector<double>(tx);

      ublas::vector<int> vhat1 = decodeBitFlipping(tx, H, iterations); // newH

      std::cout << "vhat1=" << std::endl;
      printVector<int>(vhat1);

      double rat1 = biterr(vhat1, u);
      //std::cout << "rat1=" << rat1 << std::endl;
      ber1[i] += rat1;
    }

    ber1[i] /= frames;
  }

  std::cout << std::endl;

  std::cout << "EbN0=[";
  for (unsigned int i = 0; i < 4; i++) {
    std::cout << EbN0[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "ber1=[";
  for (unsigned int i = 0; i < 4; i++) {
    std::cout << ber1[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "semilogy(EbN0, ber1, 'o-');" << std::endl;

  return 0;
}
