#ifndef SPECTRUM_FRAGMENT_UTILS_DISPATCH_H
#define SPECTRUM_FRAGMENT_UTILS_DISPATCH_H

// Both canonicalisation variants share an include guard and function names, so each is
// included in its own namespace. Callers pick one with FragmentCollection::unary_layout():
//   unary_utils: fragments alternate unary/binary slots; null (0) unary slots have no location entry
//   plain_utils: binary-only fragments (faster; wrong for fragments with unary slots)
namespace unary_utils {
#include "fragment_utils.h"
}
#undef SPECTRUM_FRAGMENT_UTILS_H
namespace plain_utils {
#include "fragment_utils_no_unaries.h"
}

// The random-selection helpers are the same in both variants.
using plain_utils::get_random_elements;
using plain_utils::direct_random_non_cyclic_selection;

#endif //SPECTRUM_FRAGMENT_UTILS_DISPATCH_H
