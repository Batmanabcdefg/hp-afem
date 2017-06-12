#include <iostream>
#include <string>
#include <getopt.h>

class Option {
protected:
  std::string key;
  std::string _val;
public:
  Option(std::string key, std::string val) : key(key), _val(val) {}
  std::string val() { return _val; }
  void setVal(std::string val) { _val = val; }
};

class Options {
protected:
  int argc;
  char **argv;

  std::map<std::string, Option> opts;
public:
  int get(std::string key) {
    if( !opts.count(key)) {
      throw runtime_error("No valid option!");
    }
    return opts[key].val();
  }

  void set(std::string key, std::string val) {
    if( !opts.count(key)) {
      opts[key] = Option(key, val);
    } else {
      opts[key].setVal(val);
    }
  }

  Options(int argc, char **argv);
};
