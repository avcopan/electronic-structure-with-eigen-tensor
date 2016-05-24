#ifndef UHF_HPP_
#define UHF_HPP_

#include <iostream> // std::cout
#include <cstdio>   // std::prinft
#include <cmath>    // std::abs
#include <tuple>    // std::tie
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../3integrals/integrals.hpp"
#include "../2helpers/helpers.hpp"

namespace uhf {

/* typedefs */
using IndexPair = Eigen::IndexPair<int>;
template<class C> using Shared      = boost::shared_ptr<C>;
template<int n>   using Contraction = std::array<IndexPair, n>;
template<int n>   using Tuple       = std::array<int      , n>;

class UHF {
private:
  int nbf_;
  int naocc_;
  int nbocc_;
  double energy_;
  double vnu_;
  Eigen::Tensor<double, 2> aC_;
  Eigen::Tensor<double, 2> bC_;
  Eigen::Tensor<double, 1> ae_;
  Eigen::Tensor<double, 1> be_;
  Eigen::Tensor<double, 2> aD_;
  Eigen::Tensor<double, 2> bD_;
  psi::Options options_;
  Shared<integrals::Integrals> integrals_;
public:
  UHF(Shared<psi::Wavefunction> wfn, psi::Options& options);
  double compute_energy();
  double get_energy() { return energy_; }
  int get_naocc() { return naocc_; }
  int get_nbocc() { return nbocc_; }
  std::pair<Eigen::Tensor<double, 1>, Eigen::Tensor<double, 2> >
    compute_natural_orbitals();
  Shared<integrals::Integrals>
    get_integrals()  { return integrals_; }
  Eigen::Tensor<double, 1>
    get_a_orbital_energies() { return ae_; }
  Eigen::Tensor<double, 1>
    get_b_orbital_energies() { return be_; }
};

} // uhf

#endif // UHF_HPP_
