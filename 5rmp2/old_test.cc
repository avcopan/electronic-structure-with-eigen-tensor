#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "../4rhf/rhf.hpp"
#include "rmp2.hpp"

template<class C> using Shared = boost::shared_ptr<C>;

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  Shared<rhf::RHF>   rhf(new rhf::RHF(wfn, options));
  Shared<rmp2::RMP2> rmp2(new rmp2::RMP2(rhf, options));
  rmp2->compute_energy();

  return wfn;
}

