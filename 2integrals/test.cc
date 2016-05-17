#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "integrals.hpp"

using EigenTensor  = Eigen::Tensor<double, 4>;
using EigenMatrix  = Eigen::Matrix<double,-1,-1>;

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  integrals::Integrals integrals(wfn, options);
  EigenMatrix S = integrals.ao_overlap();
  EigenMatrix T = integrals.ao_kinetic();
  EigenMatrix V = integrals.ao_potential();
  EigenTensor g = integrals.ao_repulsion();
  std::cout << S << std::endl;
  std::cout << T << std::endl;
  std::cout << V << std::endl;
  std::cout << g << std::endl;

  return wfn;
}

