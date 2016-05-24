#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include <boost/shared_ptr.hpp>           // boost::shared_ptr
#include <libmints/mints.h>               // psi::OneBodyAOInt, psi::TwoBodyAOInt
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../2helpers/helpers.hpp"

namespace integrals {

/* typedefs */
using IndexPair   = Eigen::IndexPair<int>;
template<int n>   using Contraction = std::array<IndexPair, n>;
template<class Class> using Shared = boost::shared_ptr<Class>;

/* classes */
class Integrals {
private:
  int nbf_;
  Eigen::Tensor<double, 2> ac_; // alpha MO coefficient matrix
  Eigen::Tensor<double, 2> bc_; // beta  MO coefficient matrix
  Eigen::Tensor<double, 2> s_;  // overlap integrals
  Eigen::Tensor<double, 2> t_;  // electronic kinetic energy operator
  Eigen::Tensor<double, 2> v_;  // electron-nuclear repulsion operator
  Eigen::Tensor<double, 4> g_;  // two-electron integrals in physicist's notation, <mu nu | rh si>
  
public:
  /* constructor */
  Integrals(Shared<psi::Wavefunction> wfn, psi::Options& options);
  int get_nbf();
  Eigen::Tensor<double, 2> get_ao_overlap();
  Eigen::Tensor<double, 2> get_ao_kinetic();
  Eigen::Tensor<double, 2> get_ao_potential();
  Eigen::Tensor<double, 2> get_ao_corehamiltonian();
  Eigen::Tensor<double, 2> get_ao_orthogonalizer();
  Eigen::Tensor<double, 4> get_ao_eri_physnotation();
  void set_mo_a_coefficients(Eigen::Tensor<double, 2>);
  void set_mo_b_coefficients(Eigen::Tensor<double, 2>);
  Eigen::Tensor<double, 2> get_mo_a_coefficients();
  Eigen::Tensor<double, 2> get_mo_b_coefficients();
  Eigen::Tensor<double, 2> get_mo_aa_overlap(); // for testing purposes
  Eigen::Tensor<double, 2> get_mo_bb_overlap(); // for testing purposes
  Eigen::Tensor<double, 2> get_mo_aa_kinetic();
  Eigen::Tensor<double, 2> get_mo_bb_kinetic();
  Eigen::Tensor<double, 2> get_mo_aa_potential();
  Eigen::Tensor<double, 2> get_mo_bb_potential();
  Eigen::Tensor<double, 2> get_mo_aa_corehamiltonian();
  Eigen::Tensor<double, 2> get_mo_bb_corehamiltonian();
  Eigen::Tensor<double, 4> get_mo_aaaa_eri_physnotation();
  Eigen::Tensor<double, 4> get_mo_abab_eri_physnotation();
  Eigen::Tensor<double, 4> get_mo_bbbb_eri_physnotation();
};

/* functions */
Eigen::Tensor<double, 2> compute_oei(psi::OneBodyAOInt*);
Eigen::Tensor<double, 4> compute_tei(psi::TwoBodyAOInt*);

} // end namespace ao_integrals

#endif // INTEGRALS_HPP_
