#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

extern "C"
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C"
psi::PsiReturnType test(psi::Options& options)
{

  /* Your code goes here */
  std::cout << "Hello world" << std::endl;

  return psi::Success;
}

