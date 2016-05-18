#include "rhf.hpp"

namespace rhf {

RHF::RHF(Shared<psi::Wavefunction> wfn, psi::Options& options) : options_(options), integrals_(wfn, options) {
  Shared<psi::Molecule> mol = wfn->molecule();
  int nelec = - mol->molecular_charge();                     // calculate the number of electrons as the difference
  for (int A = 0; A < mol->natom(); ++A) nelec += mol->Z(A); // between the total nuclear charge and the molecular charge
  // if it isn't a closed-shell molecule, complain
  if (nelec % 2 != 0) throw std::invalid_argument("RHF code requires a closed-shell molecule.");
  ndocc_  = nelec / 2;
  nbf_    = integrals_.get_nbf();
  vnu_    = mol->nuclear_repulsion_energy();
  energy_ = 0.0;
}

double RHF::compute_energy() {
  // grab integrals
  Eigen::Tensor<double, 2> x = integrals_.get_ao_orthogonalizer();
  Eigen::Tensor<double, 2> h = integrals_.get_ao_corehamiltonian();
  Eigen::Tensor<double, 4> g = integrals_.get_ao_eri_physnotation();
  // declare intermediates
  Eigen::Tensor<double, 2> d (nbf_, nbf_), v (nbf_, nbf_), f (nbf_, nbf_  ),
                           tf(nbf_, nbf_), tc(nbf_, nbf_), oc(nbf_, ndocc_);
  double energy, denergy;
  d.setZero();

  Contraction<2> jcontraction({IndexPair(3, 0), IndexPair(1, 1)});
  Contraction<2> kcontraction({IndexPair(2, 0), IndexPair(1, 1)});

  for(int iter = 0; iter < options_.get_int("MAXITER"); ++iter) {
    v  = g.contract(d, jcontraction) - g.contract(d, kcontraction) * 0.5;
    f  = h + v;
    tf = x % f % x;
    std::tie(e_, tc) = matricks::eigh(tf);
    c_ = x % tc;
    oc = c_.slice(Tuple<2>{0, 0}, Tuple<2>{nbf_, ndocc_});
    d  = (oc % matricks::transpose(oc)) * 2.0;
    energy = matricks::trace( (h + v * 0.5).eval() % d ) + vnu_;
    denergy = energy - energy_; energy_ = energy;
    std::printf("@RHF %-3d %20.15f %20.15f\n", iter, energy, denergy);
    if(std::abs(denergy) < options_.get_double("E_CONVERGENCE")) break;
  }

  return 0.0;
}

}
