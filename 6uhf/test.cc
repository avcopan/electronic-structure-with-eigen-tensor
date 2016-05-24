#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "uhf.hpp"

template<class C> using Shared = boost::shared_ptr<C>;

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  Shared<uhf::UHF> uhf(new uhf::UHF(wfn, options));
  uhf->compute_energy();
  Eigen::Tensor<double, 1> n;
  Eigen::Tensor<double, 2> C;
  std::tie(n, C) = uhf->compute_natural_orbitals();
  for(int p=0; p < uhf->get_integrals()->get_nbf(); ++p)
    std::printf("@UHF NO %-2d occupation: %20.15f\n", p+1, n(p));

  return wfn;
}

