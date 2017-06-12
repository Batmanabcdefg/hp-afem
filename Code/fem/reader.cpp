#include "reader.h"

using namespace std;

struct option Reader::options[] = {
  {"basis",     required_argument, 0, 0},
  {"mesh",      required_argument, 0, 0},
  {"rhs",       required_argument, 0, 0},
  {"degree",    required_argument, 0, 0},
  {"refines",   required_argument, 0, 0},
  {NULL, 0, NULL, 0}
};

Reader::Reader(int argc, char **argv) : argc(argc), argv(argv) {
  int c;
  int optindex;

  while((c = getopt_long(argc, argv, "", options, &optindex)) != -1) {

    cout << options[optindex].name << " " << optarg << endl;
  }
}
