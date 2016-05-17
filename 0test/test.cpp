#include <iostream>
#include <array>
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor

int main() {
  Eigen::Tensor<double, 4> g(7, 7, 7, 7);
  Eigen::Tensor<double, 2> C(7, 2);
  g.setRandom();
  C.setRandom();
  typedef Eigen::IndexPair<int> pair;
  std::array<pair, 1> p({pair(0, 0)});

  // v_ijkl = g_mnrs * C_mi * C_nj * C_rk * C_sl
  Eigen::Tensor<double, 4> v = g.contract(C, p).contract(C, p).contract(C, p).contract(C, p);

  std::cout << v << std::endl;
}
