#ifndef H_APPLY_TRIGEFFS
#define H_APPLY_TRIGEFFS

#include <cstddef>
#include <string>
#include "core/named_func.hpp"
 
namespace Higfuncs{
  // analysis trigger and its efficiency
  extern const NamedFunc get_0l_trigeff;
  extern const NamedFunc get_0l_fakemet_trigeff;
  extern const NamedFunc get_1el_trigeff;
  extern const NamedFunc get_1mu_trigeff;
  extern const NamedFunc get_2el_trigeff;
  extern const NamedFunc get_2mu_trigeff;
    }

#endif
