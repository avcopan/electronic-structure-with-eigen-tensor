#ifndef RHF_HPP_
#define RHF_HPP_

#include <iostream>
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../3integrals/integrals.hpp"
#include "../2helpers/helpers.hpp"

namespace rhf {

/* typedefs */
using IndexPair   = Eigen::IndexPair<int>;
template<class Class> using Shared      = boost::shared_ptr<Class>;
template<int naxes>   using Contraction = std::array<IndexPair, naxes>;

class RHF {
private:
  int nbf_;
  int ndocc_;
  double vnu_;
  psi::Options options_;
  integrals::Integrals integrals_;
public:
  RHF(Shared<psi::Wavefunction> wfn, psi::Options& options);
  double compute_energy();
};

}

#endif // RHF_HPP_
