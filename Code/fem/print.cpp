#include <cstdio>
#include <string>
#include <cassert>
#include "print.h"

using namespace std;

template< typename... Args >
string Print::formatted( const string &format, Args... args ) {
  size_t size = 1 + snprintf(nullptr, 0, format.c_str(), args ...);
  unique_ptr<char[]> buf(new char[size]);
  snprintf(buf.get(), size, format.c_str(), args ...);
  return string(buf.get(), buf.get() + size);
}
