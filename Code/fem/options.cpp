#include "options.h"

using namespace std;

Option defaults[] = {
  Option("basis", "Bases/degree6/"),
  Option("mesh", "Meshes/unitsquare.mesh"),
  Option("rhs", "Meshes/unitsquare.rhs"),
};

struct option Options::options[] = {
  {"basis",     required_argument, 0, 0},
  {"mesh",      required_argument, 0, 0},
  {"rhs",       required_argument, 0, 0},
  {"degree",    required_argument, 0, 0},
  {"refines",   required_argument, 0, 0},
  {NULL, 0, NULL, 0}
};

Options::Options(int argc, char **argv) : argc(argc), argv(argv) {

  int c;
  int optindex;

  while((c = getopt_long(argc, argv, "", options, &optindex)) != -1) {

    cout << options[optindex].name << " " << optarg << endl;
  }
}
