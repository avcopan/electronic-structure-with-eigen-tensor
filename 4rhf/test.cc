#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "rhf.hpp"

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
  rhf::RHF rhf(wfn, options);
  rhf.compute_energy();

  return wfn;
}

