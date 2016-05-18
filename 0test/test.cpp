#include <iostream>
#include <array>
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor

using IndexPair = Eigen::IndexPair<int>;
template<int naxes>
using Contraction = std::array<IndexPair, naxes>;

int main() {
  Eigen::Tensor<double, 4> g(5, 5, 5, 5);
  Eigen::Tensor<double, 2> D(5, 5);
  Eigen::Tensor<double, 2> f(5, 5);
  g.setRandom();
  D.setRandom();
  Contraction<2> ctr({IndexPair(1, 1), IndexPair(3, 0)});
  f = g.contract(D, ctr); // f_ij = g_ikjl * D_lk

  /* this is what I would like to avoid: */
  Eigen::Matrix<double, -1, -1> F(5, 5);
  for(int i=0; i<5; ++i)
    for(int j=0; j<5; ++j)
      F(i,j) = f(i,j);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, -1, -1> > eigensolver(F);
  Eigen::Matrix<double, -1, -1> eigenvectors = eigensolver.eigenvectors();
  Eigen::Matrix<double, -1,  1> eigenvalues  = eigensolver.eigenvalues();
  std::cout << eigenvectors << std::endl;
  std::cout << eigenvalues  << std::endl;
}
