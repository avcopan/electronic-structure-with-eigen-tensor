#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include <boost/shared_ptr.hpp>           // boost::shared_ptr
#include <libmints/mints.h>               // psi::OneBodyAOInt, psi::TwoBodyAOInt
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../2helpers/helpers.hpp"

namespace integrals {

/* typedefs */
template<class Class> using Shared = boost::shared_ptr<Class>;

/* classes */
class Integrals {
private:
  int nbf_;
  Eigen::Tensor<double, 2> s_; // overlap integrals
  Eigen::Tensor<double, 2> t_; // electronic kinetic energy operator
  Eigen::Tensor<double, 2> v_; // electron-nuclear repulsion operator
  Eigen::Tensor<double, 4> g_; // two-electron integrals in physicist's notation, <mu nu | rh si>
  
public:
  /* constructor */
  Integrals(Shared<psi::Wavefunction> wfn, psi::Options& options);
  int get_nbf() { return nbf_; }
  Eigen::Tensor<double, 2> get_ao_orthogonalizer();
  Eigen::Tensor<double, 4> get_ao_eri_physnotation() { return g_;      }
  Eigen::Tensor<double, 2> get_ao_overlap()          { return s_;      }
  Eigen::Tensor<double, 2> get_ao_kinetic()          { return t_;      }
  Eigen::Tensor<double, 2> get_ao_potential()        { return v_;      }
  Eigen::Tensor<double, 2> get_ao_corehamiltonian()  { return t_ + v_; }
};

/* functions */
Eigen::Tensor<double, 2> compute_oei(psi::OneBodyAOInt*);
Eigen::Tensor<double, 4> compute_tei(psi::TwoBodyAOInt*);

} // end namespace ao_integrals

#endif // INTEGRALS_HPP_
