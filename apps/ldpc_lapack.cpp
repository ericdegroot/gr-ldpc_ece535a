/*
ECE 535A LDPC Project Testing

Required packages:

  Fedora: boost-devel lapack-devel
  Ubuntu: libboost-all-dev liblapack-dev

Build:

  g++ -Wall -o ldpc_lapack ldpc_lapack.cpp -llapack -llapacke

Run:

  ./ldpc_lapack

 */
#include <iostream>
#include <fstream>
#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/random.hpp>

#include <lapacke/lapacke.h>

#include "test_data.h"

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

template <typename T>
void writeMatrix(const ublas::matrix<T> &A, const char *fileName) {
  std::ofstream fout(fileName, std::ofstream::out);

  for (unsigned int i = 0; i < A.size1(); i++) {
    for (unsigned int j = 0; j < A.size2(); j++) {
      fout << A(i, j);

      if (j + 1 < A.size2()) {
        fout << ",";
      }
    }

    fout << std::endl;
  }

  fout.close();
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

int checkFrame(const ublas::vector<int> &u, const ublas::matrix<int> &H, const int threshold) {
  int sNotZero = 0;
  for (unsigned int k = 0; k < H.size1(); k++) {
    int s = ublas::inner_prod(u, ublas::row(H, k)) % 2;
    if (s != 0) {
      sNotZero++;
    }

    if (sNotZero > threshold) {
      break;
    }
  }

  return sNotZero;
}

ublas::vector<int> decodeSumProductSoft(const ublas::vector<double> &rx, const ublas::matrix<int> &H, const unsigned int iterations) {
  const unsigned int m = H.size1();
  const unsigned int n = H.size2();

  ublas::vector<double> r(-rx);

  // Initialization
  ublas::matrix<double> M(ublas::zero_matrix<double>(m, n));
  for (unsigned int j = 0; j < m; j++) {
    for (unsigned int i = 0; i < n; i++) {
      if (H(j, i) != 0) {
        M(j, i) = r(i);
      }
    }
  }

  ublas::vector<int> vHat(n);
  ublas::matrix<double> E(m, n);

  for (unsigned int h = 0; h < iterations; h++) {
    // Step 1: Check messages
    for (unsigned int j = 0; j < m; j++) {
      for (unsigned int i = 0; i < n; i++) {
        if (H(j, i) != 0) {
          double T = 1.0;
          for (unsigned int k = 0; k < n; k++) {
            if (H(j, k) != 0 && k != i) {
              T *= std::tanh(M(j, k) / 2.0);
            }
          }

          E(j, i) = std::log((1.0 + T) / (1.0 - T));
        }
      }
    }

    // Test
    for (unsigned int i = 0; i < n; i++) {
      double L = 0.0;
      for (unsigned int j = 0; j < m; j++) {
        if (H(j, i) != 0) {
          L += E(j, i) + r(i);
        }
      }

      if (L <= 0) {
        vHat(i) = 1;
      } else {
        vHat(i) = 0;
      }
    }
    
    // Finished?
    if (checkFrame(vHat, H, 0) == 0) {
      break;
    }

    // Step 2: Bit messages
    for (unsigned int j = 0; j < m; j++) {
      for (unsigned int i = 0; i < n; i++) {
        if (H(j, i) != 0) {
          double T = 0.0;
          for (unsigned int k = 0; k < m; k++) {
            if (k != j && H(k, i) != 0) {
              T += E(k, i) + r(i);
            }          
          }

          M(j, i) = T;
        }
      }
    }
  }

  return vHat;
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

    // Finished?
    if (n + 1 < iterations && checkFrame(vHat, H, 0) == 0) {
      break;
    }    
  }

  return vHat;
}

ublas::vector<int> decodeBitFlipping(const ublas::vector<double> &rx, const ublas::matrix<int> &H, const unsigned int iterations) {
  // Get matrix dimensions
  const unsigned int M = H.size1();
  const unsigned int N = H.size2();

  // Prior hard-decision
  ublas::vector<int> y(rx.size());
  for (unsigned int i = 0; i < rx.size(); i++) {
    if (rx(i) < 0.0) {
      y(i) = 0;
    } else {
      y(i) = 1;
    }
  }

  ublas::vector<int> ci(y);

  // Initialization
  ublas::matrix<int> E(ublas::zero_matrix<int>(M, N));

  // Iteration
  for (unsigned int n = 0; n < iterations; n++) {
    // Calculate check node messages
    for (unsigned int i = 0; i < M; i++) {
      for (unsigned int j = 0; j < N; j++) {
        int ciMod2Sum = 0;
        for (unsigned int k = 0; k < N; k++) {
          if (k != j && H(i, k) != 0) {
            ciMod2Sum += ci(k);
          }
        }
        
        E(i, j) = ciMod2Sum % 2;
      }
    }

    // Compare check node messages with y
    for (unsigned int j = 0; j < N; j++) {
      int disagreements = 0;
      for (unsigned int i = 0; i < M; i++) {
        if (H(i, j) != 0 && E(i, j) != y(j)) {
          disagreements++;
        }
      }

      // Do the majority of messages disagree with y?
      if (disagreements > M / 2) {
        ci(j) = (y(j) + 1) % 2;
      }
    }

    // Finished?
    if (n + 1 < iterations && checkFrame(ci, H, 0) == 0) {
      break;
    }    
  }

  return ci;
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

int checkMessage(const ublas::vector<int> &u, const ublas::matrix<int> &H) {
  int sNotZero = 0;
  for (unsigned int k = 0; k < H.size1(); k++) {
    int s = ublas::inner_prod(u, ublas::row(H, k)) % 2;
    if (s != 0) {
      sNotZero++;
    }
  }

  return sNotZero;
}

#define EBN0_SIZE 35

int main() {
  boost::random::mt19937 gen;
  gen.seed(static_cast<unsigned int>(std::time(0)));

  const unsigned int M = M2;
  const unsigned int N = N2;

  const unsigned int iterations = 5;
  const unsigned int frames = 30;

  const double EbN0[EBN0_SIZE] = { -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0  };
  double ber0[EBN0_SIZE];
  double ber1[EBN0_SIZE];
  double ber2[EBN0_SIZE];
  double ber3[EBN0_SIZE];

  ublas::vector<int> fer0(EBN0_SIZE);
  ublas::vector<int> fer1(EBN0_SIZE);
  ublas::vector<int> fer2(EBN0_SIZE);
  ublas::vector<int> fer3(EBN0_SIZE);

  ublas::matrix<int> H(M, N);

  initMatrix<int>(H, hData2);

  // LU matrices
  ublas::matrix<int> L(ublas::zero_matrix<int>(M, N - M));
  ublas::matrix<int> U(ublas::zero_matrix<int>(M, N - M));

  reorderHMatrix(H, L, U);

  /*
  writeMatrix<int>(H, "H.csv");
  writeMatrix<int>(L, "Lpre.csv");
  writeMatrix<int>(U, "Upre.csv");
  */

  boost::random::uniform_int_distribution<> uniformDistribution(0, 1);
  boost::random::normal_distribution<> normalDistribution;

  //ublas::matrix<int> C(EBN0_SIZE * frames, M);

  /*
  ublas::matrix<int> dSource(M, frames);
  initMatrix<int>(dSource, dSourceData2);
  writeMatrix<int>(dSource, "dSource.csv");
  */

  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    // std::cout << "EbN0 : " << i << std::endl;

    ber0[i] = 0.0;
    ber1[i] = 0.0;
    ber2[i] = 0.0;
    ber3[i] = 0.0;

    fer0(i) = 0;
    fer1(i) = 0;
    fer2(i) = 0;
    fer3(i) = 0;

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
      // std::cout << "Frame : " << j << std::endl;

      // Encoding message
      const ublas::vector<int> c = makeParityCheck(ublas::column(dSource, j), H, L, U);
      //ublas::row(C, i * (EBN0_SIZE + 1) + j) = c;

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

      ublas::vector<int> vhat3 = decodeSumProductSoft(tx, H, iterations);
      double rat3 = biterr(vhat3, u);
      ber3[i] += rat3;

      // Check messages
      if (checkMessage(vhat0, H) > 0) {
        if (rat0 <= 0.0) {
          //std::cout << "False negative: rat0=" << rat0 << std::endl;
        }

        fer0(i)++;
      } else if (rat0 > 0) {
        //std::cout << "False positive: rat0=" << rat0 << std::endl;
      }

      if (checkMessage(vhat1, H) > 0) {
        if (rat1 <= 0.0) {
          //std::cout << "False negative: rat1=" << rat1 << std::endl;
        }

        fer1(i)++;
      } else if (rat1 > 0) {
        //std::cout << "False positive: rat1=" << rat1 << std::endl;
      }

      if (checkMessage(vhat2, H) > 0) {
        if (rat2 <= 0.0) {
          //std::cout << "False negative: rat2=" << rat2 << std::endl;
        }

        fer2(i)++;
      } else if (rat2 > 0) {
        //std::cout << "False positive: rat2=" << rat2 << std::endl;
      }

      if (checkMessage(vhat3, H) > 0) {
        if (rat3 <= 0.0) {
          //std::cout << "False negative: rat3=" << rat3 << std::endl;
        }

        fer3(i)++;
      } else if (rat3 > 0) {
        //std::cout << "False positive: rat3=" << rat3 << std::endl;
      }
    }

    ber0[i] /= frames;
    ber1[i] /= frames;
    ber2[i] /= frames;
    ber3[i] /= frames;
  }

  std::cout << std::endl;

  /*
  std::cout << "H=[" << std::endl;
  printMatrix<int>(H);
  std::cout << "];" << std::endl;
  */

  std::cout.precision(4);

  std::cout << "EbN0=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << EbN0[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "figure(1);" << std::endl;

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

  std::cout << "ber3=[";
  for (unsigned int i = 0; i < EBN0_SIZE; i++) {
    std::cout << ber3[i] << " ";
  }
  std::cout << "];" << std::endl;

  std::cout << "plot(EbN0, ber3, 'om-');" << std::endl;

  std::cout << "grid on;" << std::endl;
  std::cout << "hold off;" << std::endl;

  std::cout << "title('Bit Error Rate');" << std::endl;
  std::cout << "legend('BPSK', 'BitFlip', 'LogDomainSimple', 'SumProduct');" << std::endl;
  std::cout << "xlabel('EbN0');" << std::endl;
  std::cout << "ylabel('BER');" << std::endl;

  std::cout << "figure(2);" << std::endl;

  std::cout << "fer0=[" << std::endl;
  printVector<int>(fer0);
  std::cout << "];" << std::endl;  

  std::cout << "plot(EbN0, fer0, 'or:');" << std::endl;

  std::cout << "hold;" << std::endl;

  std::cout << "fer1=[" << std::endl;
  printVector<int>(fer1);
  std::cout << "];" << std::endl;  

  std::cout << "plot(EbN0, fer1, 'og-');" << std::endl;

  std::cout << "fer2=[" << std::endl;
  printVector<int>(fer2);
  std::cout << "];" << std::endl;  

  std::cout << "plot(EbN0, fer2, 'ob-');" << std::endl;

  std::cout << "fer3=[" << std::endl;
  printVector<int>(fer3);
  std::cout << "];" << std::endl;  

  std::cout << "plot(EbN0, fer3, 'om-');" << std::endl;

  std::cout << "grid on;" << std::endl;
  std::cout << "hold off;" << std::endl;

  std::cout << "title('Frame Errors');" << std::endl;
  std::cout << "legend('BPSK', 'BitFlip', 'LogDomainSimple', 'SumProduct');" << std::endl;
  std::cout << "xlabel('EbN0');" << std::endl;
  std::cout << "ylabel('FER');" << std::endl;

  /*
  writeMatrix<int>(L, "Lpost.csv");
  writeMatrix<int>(U, "Upost.csv");
  //writeMatrix<int>(dSource, "dSource.csv");
  writeMatrix<int>(C, "C.csv");
  */

  return 0;
}
