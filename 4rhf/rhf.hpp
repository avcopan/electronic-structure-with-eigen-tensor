#ifndef RHF_HPP_
#define RHF_HPP_

#include <iostream>                       // std::cout
#include <cstdio>                         // std::printf
#include <tuple>                          // std::tie
#include <cmath>                          // std::abs
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../3integrals/integrals.hpp"
#include "../2helpers/helpers.hpp"

namespace rhf {

/* typedefs */
using IndexPair   = Eigen::IndexPair<int>;
template<class C> using Shared      = boost::shared_ptr<C>;
template<int n>   using Contraction = std::array<IndexPair, n>;
template<int n>   using Tuple       = std::array<int      , n>;

class RHF {
private:
  int nbf_;
  int ndocc_;
  double energy_;
  double vnu_;
  Eigen::Tensor<double, 2> c_;
  Eigen::Tensor<double, 1> e_;
  psi::Options options_;
  integrals::Integrals integrals_;
public:
  RHF(Shared<psi::Wavefunction> wfn, psi::Options& options);
  double compute_energy();
};

}

#endif // RHF_HPP_
