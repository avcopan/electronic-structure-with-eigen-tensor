#ifndef RMP2_HPP_
#define RMP2_HPP_

#include <iostream>  // std::cout
#include <cstdio>    // std::printf
#include <liboptions/liboptions.h>        // psi::Options
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../4rhf/rhf.hpp"

namespace rmp2 {

/* typedefs */
using IndexPair   = Eigen::IndexPair<int>;
template<class C> using Shared      = boost::shared_ptr<C>;
template<int n>   using Contraction = std::array<IndexPair, n>;
template<int n>   using Tuple       = std::array<int      , n>;

class RMP2 {
private:
  int nbf_;
  int naocc_;
  int navir_;
  double hfenergy_;
  double cenergy_;
  psi::Options options_;
  Shared<rhf::RHF> rhf_;
  Eigen::Tensor<double, 1> e_;
  Eigen::Tensor<double, 4> g_;
public:
  RMP2(Shared<rhf::RHF>, psi::Options&);
  double compute_energy();
};

}

#endif // RMP2_HPP_
