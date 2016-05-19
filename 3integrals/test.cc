#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "integrals.hpp"

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  integrals::Integrals integrals(wfn, options);
  Eigen::Tensor<double, 2> S = integrals.get_ao_overlap();
  Eigen::Tensor<double, 2> X = integrals.get_ao_orthogonalizer();
  Eigen::Tensor<double, 2> T = integrals.get_ao_kinetic();
  Eigen::Tensor<double, 2> V = integrals.get_ao_potential();
  Eigen::Tensor<double, 4> g = integrals.get_ao_eri_physnotation();
  std::cout << S << std::endl;
  std::cout << T << std::endl;
  std::cout << V << std::endl;
  //std::cout << g << std::endl;
  std::cout << X % X % S << std::endl;
  int nbf = integrals.get_nbf();
  std::cout << nbf << std::endl;
  integrals.set_mo_coefficients(X);
  std::cout << integrals.get_mo_coefficients() << std::endl;

  return wfn;
}

