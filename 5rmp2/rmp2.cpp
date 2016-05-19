#include "rmp2.hpp"

namespace rmp2 {

RMP2::RMP2(Shared<rhf::RHF> rhf, psi::Options& options) : rhf_(rhf), options_(options) {
  hfenergy_ = rhf_->get_energy();
  nbf_      = rhf_->get_integrals()->get_nbf();
  naocc_    = rhf_->get_naocc();
  navir_    = nbf_ - naocc_;
  e_        = rhf_->get_orbital_energies();
  g_        = rhf_->get_integrals()->get_mo_eri_physnotation();
}

double RMP2::compute_energy() {
/* // I tried doing this using tensor operations rather than for loops and it's a total nightmare
  // 1. build t_ijab = g_ijab / (e_i + e_j - e_a - e_b)
  Eigen::Tensor<double, 4> t(naocc_, naocc_, navir_, navir_);
  std::pair<Tuple<1>, Tuple<1> > o(Tuple<1>{{     0}}, Tuple<1>{{naocc_}});
  std::pair<Tuple<1>, Tuple<1> > v(Tuple<1>{{naocc_}}, Tuple<1>{{navir_}});
  std::pair<Tuple<4>, Tuple<4> > oovv(Tuple<4>{{     0,      0, naocc_, naocc_}},
                                      Tuple<4>{{naocc_, naocc_, navir_, navir_}});
  Tuple<4> oxxx{{naocc_,1,1,1}}, xoxx{{1,naocc_,1,1}},
           xxvx{{1,1,navir_,1}}, xxxv{{1,1,1,navir_}};
  t = g_.slice(oovv.first, oovv.second) / ( e_.slice(o.first, o.second).reshape(oxxx).broadcast(oovv.second)
                                          + e_.slice(o.first, o.second).reshape(xoxx).broadcast(oovv.second)
                                          - e_.slice(v.first, v.second).reshape(xxvx).broadcast(oovv.second)
                                          - e_.slice(v.first, v.second).reshape(xxxv).broadcast(oovv.second) );
  // 2. determine energy as E_c = t_ijab * (g_ijab - 1/2 * g_ijba)
  Tuple<4> pqsr{{0,1,3,2}};
  Eigen::Tensor<double,1> Ec = ( t * (g_.slice(oovv.first, oovv.second) - g_.slice(oovv.first, oovv.second) * 0.5) ).sum();
*/
  cenergy_ = 0.0;
  for(int i=0; i < naocc_; ++i)
    for(int j=0; j < naocc_; ++j)
      for(int a=naocc_; a < nbf_; ++a)
        for(int b=naocc_; b < nbf_; ++b)
          cenergy_ += g_(i, j, a, b) * ( 2 * g_(i, j, a, b) - g_(i, j, b, a) ) / (e_(i) + e_(j) - e_(a) - e_(b));
  std::printf("@RMP2 correlation energy: %20.15f\n", cenergy_);
  std::printf("@RMP2 total energy:       %20.15f\n", cenergy_ + hfenergy_);
}

}
