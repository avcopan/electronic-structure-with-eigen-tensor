#include <unsupported/Eigen/CXX11/Tensor>    // Eigen::Tensor
#include <iostream>                          // std::cout
#include <tuple>                             // std::tie
#include "helpers.hpp"

using IndexPair = Eigen::IndexPair<int>;
template<int naxes>
using Contraction = std::array<IndexPair, naxes>;

int main() {
  Eigen::Tensor<double, 2> m(5, 5);
  m.setRandom();
  m = m % matricks::transpose( m ); // make a symmetric, positive semi-definite matrix

  std::cout << m << std::endl;
  
  Eigen::Tensor<double, 2> minv = matricks::inverse( m );
  Eigen::Tensor<double, 2> msqrt = matricks::sqrth( m );

  std::cout << m % minv << std::endl;
  std::cout << msqrt << std::endl;
  std::cout << msqrt % msqrt % minv << std::endl;

  Eigen::Tensor<double, 2> U(5, 5);
  Eigen::Tensor<double, 1> e(5);
  std::tie(e, U) = matricks::eigh( m );

  std::cout << e << std::endl;
  std::cout << U << std::endl;

  std::cout << matricks::trace( msqrt % msqrt % minv ) << std::endl;
}
