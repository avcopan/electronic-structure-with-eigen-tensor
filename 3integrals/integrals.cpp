#include "integrals.hpp"

namespace integrals {

Integrals::Integrals(Shared<psi::Wavefunction> wfn, psi::Options& options) {
  // 1. find out which basis set we're using
  std::string basisname = options.get_str("BASIS");
  // 2. build integral factory
  Shared<psi::Molecule>        mol   = wfn->molecule();
  Shared<psi::BasisSet>        basis = psi::BasisSet::pyconstruct_orbital(mol, "BASIS", basisname);
  Shared<psi::IntegralFactory> integralfactory(new psi::IntegralFactory(basis, basis, basis, basis));
  // 3. compute integrals and store them in Eigen arrays
  nbf_ = basis->nbf();
  s_   = compute_oei(integralfactory->ao_overlap()  );
  t_   = compute_oei(integralfactory->ao_kinetic()  );
  v_   = compute_oei(integralfactory->ao_potential());
  Eigen::Tensor<double, 4> g = compute_tei(integralfactory->eri());
  std::array<int,4> phys({0,2,1,3});
  g_ = g.shuffle(phys);
}

int Integrals::get_nbf() {
  return nbf_;
}

Eigen::Tensor<double, 2> Integrals::get_ao_overlap() {
  return s_;
}

Eigen::Tensor<double, 2> Integrals::get_ao_kinetic() {
  return t_;
}

Eigen::Tensor<double, 2> Integrals::get_ao_potential() {
  return v_;
}

Eigen::Tensor<double, 4> Integrals::get_ao_eri_physnotation() {
  return g_;
}

Eigen::Tensor<double, 2> Integrals::get_ao_corehamiltonian() {
  return t_ + v_;
}

Eigen::Tensor<double, 2> Integrals::get_ao_orthogonalizer() {
  return matricks::inverse( matricks::sqrth(s_) );
}

void Integrals::set_mo_coefficients(Eigen::Tensor<double, 2> c) {
  c_ = c;
}

Eigen::Tensor<double, 2> Integrals::get_mo_coefficients() {
  return c_;
}

Eigen::Tensor<double, 2> Integrals::get_mo_overlap() {
  return matricks::transpose(c_) % s_ % c_;
}

Eigen::Tensor<double, 2> Integrals::get_mo_kinetic() {
  return matricks::transpose(c_) % t_ % c_;
}

Eigen::Tensor<double, 2> Integrals::get_mo_potential() {
  return matricks::transpose(c_) % v_ % c_;
}

Eigen::Tensor<double, 2> Integrals::get_mo_corehamiltonian() {
  return matricks::transpose(c_) % (t_ + v_).eval() % c_;
}

Eigen::Tensor<double, 4> Integrals::get_mo_eri_physnotation() {
  Contraction<1> ctr({IndexPair(0, 0)});
  return g_.contract(c_, ctr).contract(c_, ctr).contract(c_, ctr).contract(c_, ctr);
}

Eigen::Tensor<double, 2> compute_oei(psi::OneBodyAOInt* O) {
  boost::shared_ptr<psi::BasisSet> bs1 = O->basis1(); int dim1 = bs1->nbf();
  boost::shared_ptr<psi::BasisSet> bs2 = O->basis2(); int dim2 = bs2->nbf();
 
  Eigen::Tensor<double, 2> A(dim1, dim2);
  double* pA = A.data();

  const double* buffer = O->buffer();

  for (int M = 0; M < bs1->nshell(); ++M) {
    for (int N = 0; N < bs2->nshell(); ++N) {

      O->compute_shell(M, N);

      int index = 0;
      for (int m = 0; m < bs1->shell(M).nfunction(); ++m)
        for (int n = 0; n < bs2->shell(N).nfunction(); ++n)
          pA[(bs1->shell(M).function_index() + m)*1   +
             (bs2->shell(N).function_index() + n)*dim1] = buffer[index++];
    }
  }

  return A;
}

Eigen::Tensor<double, 4> compute_tei(psi::TwoBodyAOInt* T) {

  boost::shared_ptr<psi::BasisSet> bs1 = T->basis1(); int dim1 = bs1->nbf();
  boost::shared_ptr<psi::BasisSet> bs2 = T->basis2(); int dim2 = bs2->nbf();
  boost::shared_ptr<psi::BasisSet> bs3 = T->basis3(); int dim3 = bs3->nbf();
  boost::shared_ptr<psi::BasisSet> bs4 = T->basis4(); int dim4 = bs4->nbf();

  Eigen::Tensor<double, 4> A(dim1, dim2, dim3, dim4);
  double* pA = A.data();

  const double* buffer = T->buffer();

  for (int M = 0; M < bs1->nshell(); ++M) {
    for (int N = 0; N < bs2->nshell(); ++N) {
      for (int P = 0; P < bs3->nshell(); ++P) {
        for (int Q = 0; Q < bs4->nshell(); ++Q) {

          T->compute_shell(M,N,P,Q);

          int index = 0;
          for (int m = 0; m < bs1->shell(M).nfunction(); ++m)
            for (int n = 0; n < bs2->shell(N).nfunction(); ++n)
              for (int p = 0; p < bs3->shell(P).nfunction(); ++p)
                for (int q = 0; q < bs4->shell(Q).nfunction(); ++q)
                  pA[(bs1->shell(M).function_index() + m)*1             +
                     (bs2->shell(N).function_index() + n)*dim1          +
                     (bs3->shell(P).function_index() + p)*dim1*dim2     +
                     (bs4->shell(Q).function_index() + q)*dim1*dim2*dim3] = buffer[index++];

        }
      }
    }
  }

  return A;
}

} // end namespace ao_integrals
