#ifndef TYPE_HPP
#define TYPE_HPP

#include <string>
#include <typeinfo>

namespace demangler{

std::string demangle(const char* name);

template <class T>
std::string type(const T& t) {
  return demangle(typeid(t).name());
}

}

#endif
