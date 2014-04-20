/*
ECE 535A LDPC Project Testing

Required packages:

  Fedora: boost-devel suitesparse-devel
  Ubuntu: libsuitesparse-dev libboost-all-dev

Build:

  g++ -Wall -o ldpc_umfpack ldpc_umfpack.cpp -lumfpack -lamd

Run:

  ./ldpc_umfpack

 */
#include <iostream>
#include <limits>

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
      std::cout << A(i, j);

      if (j + 1 < A.size2()) {
        std::cout << ", ";
      }
    }

    if (i + 1 < A.size1()) {
      std::cout << ";";
    }

    std::cout << std::endl;
  } 
}

template <typename T>
void printVector(const ublas::vector<T> &u) {
  for (unsigned int i = 0; i < u.size(); i++) {
    std::cout << u(i);

    if (i + 1 < u.size()) {
      std::cout << ", ";
    }
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

void reorderHMatrix(ublas::matrix<int> &H, ublas::matrix<int> &L, ublas::matrix<int> &U) {
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

ublas::vector<int> makeParityCheck(const ublas::vector<int> &dSource, const ublas::matrix<int> &H, const ublas::matrix<int> &L, const ublas::matrix<int> &U) {
  // Get matrix dimensions
  const unsigned int M = H.size1();
  const unsigned int N = H.size2();

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

ublas::vector<int> decodeLogDomainSimple(const ublas::vector<double> &rx, const ublas::matrix<int> &H, const unsigned int iterations) {
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
    std::cout << "Iteration : " << n << std::endl;

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
  }

  return vHat;
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

ublas::vector<int> decodeBPSK(const ublas::vector<double> &rx) {
  ublas::vector<int> vHat(rx.size());
  for (unsigned int i = 0; i < rx.size(); i++) {
    vHat(i) = 0.5 * (sign(rx(i)) + 1);
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

/*
const int dSourceData3[] = {
  0,0,1,1,0,0,1,0,0,0,1,1,1,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1,1,0,
  0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,0,0,1,0,
  1,0,0,1,1,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,
  1,1,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1,1,0,
  1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,1,1,0,0,1,0,0,0,
  0,1,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,
  1,0,1,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,0,
  1,0,1,0,0,1,0,1,1,1,0,1,0,0,0,1,1,1,1,0,0,0,1,0,1,0,1,0,1,1,
  1,1,1,1,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,1,
  1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,0,1,1,1,1,1,1,0,1,0,0,1,1,1,
  0,0,1,1,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,1,0,1,0,1,0,0,0,1,1,0,
  0,1,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,0,1,0,1,1,1,1,0,0,1,1,0,
  0,0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,
  1,0,1,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,
  1,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,1,0,0,
  0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,0,0,0,1,0,1,0,1,0,0,
  0,0,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,0,0,1,
  1,0,0,1,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,
  1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,0,
  1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,
  1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,0,1,
  1,0,0,1,1,1,0,0,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,
  0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,
  1,1,1,0,1,1,1,1,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,
  0,1,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,
  0,0,1,1,1,1,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,1,
  0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,0,1,1,
  0,1,0,1,1,0,1,1,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,1,0,1,
  1,1,1,0,1,0,0,0,0,0,1,1,0,1,1,0,1,0,0,1,0,0,0,1,1,0,0,0,0,1,
  0,0,1,0,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,
  0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,0,
  1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,
  0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,
  1,1,1,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,
  0,1,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,0,1,1,
  1,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,
  1,0,0,1,1,1,0,0,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,
  0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,1,1,1,
  1,1,1,0,0,0,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1,0,
  1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,
  0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1,1,1,
  1,0,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,1,0,
  0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,1,1,0,1,1,
  0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,1,1,0,1,1,1,0,0,1,1,1,
  1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,0,
  1,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,0,
  0,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,0,0,0,
  1,0,0,1,0,1,0,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0,1,
  1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,1,0,1,1,0,1,0,1,1,1,1,
  0,1,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,0,1,0,1,1
};
*/

#define EBN0_SIZE 29

int main() {
  boost::random::mt19937 gen;
  gen.seed(static_cast<unsigned int>(std::time(0)));

  const unsigned int M = M3;
  const unsigned int N = N3;

  const unsigned int iterations = 5;
  const unsigned int frames = 30;

  const double EbN0[EBN0_SIZE] = { -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0  };
  double ber0[EBN0_SIZE];
  double ber1[EBN0_SIZE];
  double ber2[EBN0_SIZE];

  ublas::matrix<int> H(M, N);

  initMatrix<int>(H, hData3);

  // LU matrices
  ublas::matrix<int> L(ublas::zero_matrix<int>(M, N - M));
  ublas::matrix<int> U(ublas::zero_matrix<int>(M, N - M));

  reorderHMatrix(H, L, U);

  boost::random::uniform_int_distribution<> uniformDistribution(0, 1);
  boost::random::normal_distribution<> normalDistribution;

  ublas::matrix<double> N0(EBN0_SIZE, frames);

  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    ber0[i] = 0.0;
    ber1[i] = 0.0;
    ber2[i] = 0.0;

    // Make random data (0/1)
    ublas::matrix<int> dSource(M, frames);
    for (unsigned int j = 0; j < frames; j++) {
      for (unsigned int k = 0; k < M; k++) {
        dSource(k, j) = uniformDistribution(gen);
        //dSource(k, j) = dSourceData3[j * frames + k];
      }
    }

    //std::cout << "dSource=" << std::endl;
    //printMatrix<int>(dSource);

    for (unsigned int j = 0; j < frames; j++) {
      std::cout << "Frame : " << j << std::endl;

      // Encoding message
      const ublas::vector<int> c = makeParityCheck(ublas::column(dSource, j), H, L, U);
      
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
      N0(i, j) = 1 / (std::exp(EbN0[i] * std::log(10) / 10));

      //std::cout << "N0=" << N0 << std::endl;

      ublas::vector<double> tx(bpskMod);
      for (unsigned int k = 0; k < tx.size(); k++) {
        tx(k) += std::sqrt(N0(i, j)) * normalDistribution(gen);
      }

      //std::cout << "tx=" << std::endl;
      //printVector<double>(tx);

      ublas::vector<int> vhat0 = decodeBPSK(tx);
      double rat0 = biterr(vhat0, u);
      ber0[i] += rat0;

      ublas::vector<int> vhat1 = decodeBitFlipping(tx, H, iterations); // newH

      //std::cout << "vhat1=" << std::endl;
      //printVector<int>(vhat1);

      double rat1 = biterr(vhat1, u);
      //std::cout << "rat1=" << rat1 << std::endl;
      ber1[i] += rat1;

      ublas::vector<int> vhat2 = decodeLogDomainSimple(tx, H, iterations);
      double rat2 = biterr(vhat2, u);
      ber2[i] += rat2;
    }

    ber0[i] /= frames;
    ber1[i] /= frames;
    ber2[i] /= frames;
  }

  std::cout << std::endl;

  std::cout.precision(4);

  std::cout << "N0=[" << std::endl;
  printMatrix<double>(N0);
  std::cout << "];" << std::endl;

  std::cout << "EbN0=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << EbN0[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "ber0=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << ber0[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "plot(EbN0, ber0, 'or--');" << std::endl;

  std::cout << "hold;" << std::endl;

  std::cout << "ber1=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << ber1[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "plot(EbN0, ber1, 'og-');" << std::endl;

  std::cout << "ber2=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << ber2[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "plot(EbN0, ber2, 'ob-');" << std::endl;

  std::cout << "grid on;" << std::endl;
  std::cout << "hold off;" << std::endl;

  std::cout << "legend('BPSK', 'BitFlip', 'LogDomainSimple');" << std::endl;
  std::cout << "xlabel('EbN0');" << std::endl;
  std::cout << "ylabel('BER');" << std::endl;

  return 0;
}
