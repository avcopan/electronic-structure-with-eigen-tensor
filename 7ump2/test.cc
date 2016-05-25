#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "../6uhf/uhf.hpp"
#include "ump2.hpp"

template<class C> using Shared = boost::shared_ptr<C>;

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::SharedWavefunction test(psi::SharedWavefunction wfn, psi::Options& options)
{

  /* Your code goes here */
  Shared<uhf::UHF>   uhf(new uhf::UHF(wfn, options));
  uhf->compute_energy();
  Shared<ump2::UMP2> ump2(new ump2::UMP2(uhf, options));
  ump2->compute_energy();

  return wfn;
}

