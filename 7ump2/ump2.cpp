#include "ump2.hpp"

namespace ump2 {

UMP2::UMP2(Shared<uhf::UHF> uhf, psi::Options& options) : uhf_(uhf), options_(options) {
  hfenergy_ = uhf_->get_energy();
  nbf_      = uhf_->get_integrals()->get_nbf();
  naocc_    = uhf_->get_naocc();
  nbocc_    = uhf_->get_nbocc();
  ae_       = uhf_->get_a_orbital_energies();
  be_       = uhf_->get_b_orbital_energies();
  aag_      = uhf_->get_integrals()->get_mo_aaaa_eri_physnotation();
  abg_      = uhf_->get_integrals()->get_mo_abab_eri_physnotation();
  bbg_      = uhf_->get_integrals()->get_mo_bbbb_eri_physnotation();
}

double UMP2::compute_energy() {
  cenergy_ = 0.0;
  for(int i=0; i < naocc_; ++i)
    for(int j=0; j < naocc_; ++j)
      for(int a=naocc_; a < nbf_; ++a)
        for(int b=naocc_; b < nbf_; ++b)
          cenergy_ += (1./4) * (aag_(i, j, a, b) - aag_(i, j, b, a)) * (aag_(i, j, a, b) - aag_(i, j, b, a))
                             / (ae_(i) + ae_(j) - ae_(a) - ae_(b));
  for(int i=0; i < naocc_; ++i)
    for(int J=0; J < nbocc_; ++J)
      for(int a=naocc_; a < nbf_; ++a)
        for(int B=nbocc_; B < nbf_; ++B)
          cenergy_ +=    abg_(i, J, a, B) * abg_(i, J, a, B)
                      / (ae_(i) + be_(J) - ae_(a) - be_(B));
  for(int I=0; I < nbocc_; ++I)
    for(int J=0; J < nbocc_; ++J)
      for(int A=nbocc_; A < nbf_; ++A)
        for(int B=nbocc_; B < nbf_; ++B)
          cenergy_ += (1./4) * (bbg_(I, J, A, B) - bbg_(I, J, B, A)) * (bbg_(I, J, A, B) - bbg_(I, J, B, A))
                             / (be_(I) + be_(J) - be_(A) - be_(B));
  std::printf("@UMP2 correlation energy: %20.15f\n", cenergy_);
  std::printf("@UMP2 total energy:       %20.15f\n", cenergy_ + hfenergy_);
}

}
