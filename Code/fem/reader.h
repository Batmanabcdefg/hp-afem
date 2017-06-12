#include <iostream>
#include <string>
#include <getopt.h>

class Reader {
protected:
  int argc;
  char **argv;
public:
  static struct option options[];
  Reader(int argc, char **argv);
};
