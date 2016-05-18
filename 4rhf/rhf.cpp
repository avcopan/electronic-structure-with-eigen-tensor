#include "rhf.hpp"

namespace rhf {

RHF::RHF(Shared<psi::Wavefunction> wfn, psi::Options& options) : options_(options), integrals_(wfn, options) {
  nbf_ = integrals_.get_nbf();
  // grab the nuclear repulsion energy from the molecule object
  Shared<psi::Molecule> mol = wfn->molecule();
  vnu_  = mol->nuclear_repulsion_energy();
  // calculate the number of electrons (= 2 * the number of occupied orbitals) as the difference
  // between the sum of nuclear charges and the overall molecular charge
  int nelec = - mol->molecular_charge();
  for (int A = 0; A < mol->natom(); ++A) nelec += mol->Z(A);
  ndocc_ = nelec / 2;
  // if it isn't a closed-shell molecule, complain
  if (nelec % 2 != 0) throw std::invalid_argument("RHF code requires a closed-shell molecule.");

}

double RHF::compute_energy() {
  // grab integrals
  Eigen::Tensor<double, 2> x = integrals_.get_ao_orthogonalizer();
  Eigen::Tensor<double, 2> h = integrals_.get_ao_corehamiltonian();
  Eigen::Tensor<double, 4> g = integrals_.get_ao_eri_physnotation();
  // declare intermediates
  Eigen::Tensor<double, 2> d(nbf_, nbf_); d.setZero();
  Eigen::Tensor<double, 2> f(nbf_, nbf_);
  Eigen::Tensor<double, 2> c(nbf_, nbf_);
  Eigen::Tensor<double, 2> tf(nbf_, nbf_);
  Eigen::Tensor<double, 2> tc(nbf_, nbf_);

  Contraction<2> jcontraction({IndexPair(3, 0), IndexPair(1, 1)});
  Contraction<2> kcontraction({IndexPair(2, 0), IndexPair(1, 1)});

  for(int iter = 0; iter < 1; ++iter)  { // options_.get_int("MAXITER"); ++iter) {
    f  = h;
    f += g.contract(d, jcontraction);
    f -= g.contract(d, kcontraction) * 0.5;
    tf = x % f % x;
  }

  std::cout << f << std::endl;

  return 0.0;
}

}
