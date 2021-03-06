#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "rhf.hpp"

template<class C> using Shared = boost::shared_ptr<C>;

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  Shared<rhf::RHF> rhf(new rhf::RHF(wfn, options));
  rhf->compute_energy();

  return wfn;
}

