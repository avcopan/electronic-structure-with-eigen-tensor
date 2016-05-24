#include "uhf.hpp"

namespace uhf {

UHF::UHF(Shared<psi::Wavefunction> wfn, psi::Options& options) : options_(options), integrals_(new integrals::Integrals(wfn, options)) {
  Shared<psi::Molecule> mol = wfn->molecule();
  int nelec = - mol->molecular_charge();
  for (int A = 0; A < mol->natom(); ++A) nelec += mol->Z(A);
  int nsocc = mol->multiplicity() - 1;
  int ndocc = (nelec - nsocc) / 2;
  naocc_ = ndocc + nsocc;
  nbocc_ = ndocc;
  nbf_ = integrals_->get_nbf();
  vnu_ = mol->nuclear_repulsion_energy();
  energy_ = 0.0;
}

double UHF::compute_energy() {
  // grab integrals
  Eigen::Tensor<double, 2> X = integrals_->get_ao_orthogonalizer();
  Eigen::Tensor<double, 2> H = integrals_->get_ao_corehamiltonian();
  Eigen::Tensor<double, 4> G = integrals_->get_ao_eri_physnotation();
  // declare intermediates
  Eigen::Tensor<double, 2> aF, bF, aV, bV, taC, tbC, taF, tbF, oaC, obC;
  double energy, denergy;
  aD_ = Eigen::Tensor<double, 2>(nbf_, nbf_); aD_.setZero();
  bD_ = Eigen::Tensor<double, 2>(nbf_, nbf_); bD_.setZero();

  Contraction<2> jcontraction({IndexPair(3, 0), IndexPair(1, 1)});
  Contraction<2> kcontraction({IndexPair(2, 0), IndexPair(1, 1)});

  for(int iter = 0; iter < options_.get_int("MAXITER"); ++iter) {
    aV = G.contract(aD_ + bD_, jcontraction) - G.contract(aD_, kcontraction);
    bV = G.contract(aD_ + bD_, jcontraction) - G.contract(bD_, kcontraction);
    aF = H + aV;
    bF = H + bV;
    taF = X % aF % X;
    tbF = X % bF % X;
    std::tie(ae_, taC) = matricks::eigh(taF);
    std::tie(be_, tbC) = matricks::eigh(tbF);
    aC_ = X % taC;
    bC_ = X % tbC;
    oaC = aC_.slice(Tuple<2>{0, 0}, Tuple<2>{nbf_, naocc_});
    obC = bC_.slice(Tuple<2>{0, 0}, Tuple<2>{nbf_, nbocc_});
    aD_ = (oaC % matricks::transpose(oaC));
    bD_ = (obC % matricks::transpose(obC));
    energy = matricks::trace( (H + aV * 0.5).eval() % aD_ ) + matricks::trace( (H + bV * 0.5).eval() % bD_ ) + vnu_;
    denergy = energy - energy_; energy_ = energy;
    std::printf("@UHF %-3d %20.15f %20.15f\n", iter, energy, denergy);
    if(std::abs(denergy) < options_.get_double("E_CONVERGENCE")) break;
  }

  integrals_->set_mo_a_coefficients(aC_);
  integrals_->set_mo_b_coefficients(bC_);

  return 0.0;
}

std::pair<Eigen::Tensor<double, 1>, Eigen::Tensor<double, 2> > UHF::compute_natural_orbitals() {
  Eigen::Tensor<double, 2> X, M, tD, tC, C;
  Eigen::Tensor<double, 1> n;

  X  = integrals_->get_ao_orthogonalizer();
  M  = matricks::sqrth( integrals_->get_ao_overlap() );
  tD = M % (aD_ + bD_).eval() % M;
  std::tie(n, tC) = matricks::eigh(tD);
  C = X % tC;
  return std::make_pair(n, C);
}

} // uhf
