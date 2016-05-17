#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include <boost/shared_ptr.hpp>           // boost::shared_ptr
#include <libmints/mints.h>               // psi::OneBodyAOInt, psi::TwoBodyAOInt
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor

namespace integrals {

/* typedefs */
template<class Class> using Shared = boost::shared_ptr<Class>;
using EigenTensor = Eigen::Tensor<double, 4>;
using EigenMatrix = Eigen::Matrix<double,-1,-1>;

/* classes */
class Integrals {
private:
  EigenMatrix s_; // overlap integrals
  EigenMatrix t_; // electronic kinetic energy operator
  EigenMatrix v_; // electron-nuclear repulsion operator
  EigenTensor g_; // two-electron integrals in physicist's notation, <mu nu | rh si>
  
public:
  /* constructor */
  Integrals(Shared<psi::Wavefunction> wfn, psi::Options& options);
  EigenMatrix ao_overlap()   { return s_; }
  EigenMatrix ao_kinetic()   { return t_; }
  EigenMatrix ao_potential() { return v_; }
  EigenTensor ao_repulsion() { return g_; }
};

/* functions */
EigenMatrix compute_oei(psi::OneBodyAOInt*);
EigenTensor compute_tei(psi::TwoBodyAOInt*);

} // end namespace ao_integrals

#endif // INTEGRALS_HPP_
