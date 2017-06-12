#pragma GCC diagnostic ignored "-Wformat"
#include <cstdio>
#include <string>
#include <cassert>
#include <ostream>

class Print {
public:
  template<typename... Args>
  static std::string formatted(const std::string &format, Args... args) {
    size_t size = 1 + std::snprintf(nullptr, 0, format.c_str(), args ...);
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size-1);
  }
};
