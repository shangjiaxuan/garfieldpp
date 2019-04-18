#ifndef G_UTILITIES_H
#define G_UTILITIES_H

#include <algorithm>
#include <functional>
#include <string>

namespace Garfield {

inline void ltrim(std::string& line) {
  line.erase(line.begin(),
             find_if(line.begin(), line.end(),
                     std::not1(std::ptr_fun<int, int>(std::isspace))));
}
}

#endif
