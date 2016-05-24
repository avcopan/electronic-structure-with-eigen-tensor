#include "rhf.hpp"

namespace rhf {

RHF::RHF(Shared<psi::Wavefunction> wfn, psi::Options& options) : options_(options), integrals_(new integrals::Integrals(wfn, options)) {
  Shared<psi::Molecule> mol = wfn->molecule();
  int nelec = - mol->molecular_charge();                     // calculate the number of electrons as the difference
  for (int A = 0; A < mol->natom(); ++A) nelec += mol->Z(A); // between the total nuclear charge and the molecular charge
  // if it isn't a closed-shell molecule, complain
  if (nelec % 2 != 0) throw std::invalid_argument("RHF code requires a closed-shell molecule.");
  naocc_  = nelec / 2;
  nbf_    = integrals_->get_nbf();
  vnu_    = mol->nuclear_repulsion_energy();
  energy_ = 0.0;
}

double RHF::compute_energy() {
  // grab integrals
  Eigen::Tensor<double, 2> X = integrals_->get_ao_orthogonalizer();
  Eigen::Tensor<double, 2> H = integrals_->get_ao_corehamiltonian();
  Eigen::Tensor<double, 4> G = integrals_->get_ao_eri_physnotation();
  // declare intermediates
  Eigen::Tensor<double, 2> V, F, tF, tC, oC;
  Eigen::Tensor<double, 2> D(nbf_, nbf_); D.setZero();
  double energy, denergy;

  Contraction<2> jcontraction({IndexPair(3, 0), IndexPair(1, 1)});
  Contraction<2> kcontraction({IndexPair(2, 0), IndexPair(1, 1)});

  for(int iter = 0; iter < options_.get_int("MAXITER"); ++iter) {
    V  = G.contract(D, jcontraction) - G.contract(D, kcontraction) * 0.5;
    F  = H + V;
    tF = X % F % X;
    std::tie(e_, tC) = matricks::eigh(tF);
    C_ = X % tC;
    oC = C_.slice(Tuple<2>{0, 0}, Tuple<2>{nbf_, naocc_});
    D  = (oC % matricks::transpose(oC)) * 2.0;
    energy = matricks::trace( (H + V * 0.5).eval() % D ) + vnu_;
    denergy = energy - energy_; energy_ = energy;
    std::printf("@RHF %-3d %20.15f %20.15f\n", iter, energy, denergy);
    if(std::abs(denergy) < options_.get_double("E_CONVERGENCE")) break;
  }

  integrals_->set_mo_a_coefficients(C_);

  return 0.0;
}

}
