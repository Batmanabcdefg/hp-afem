#include <iostream>
#include <string>
#include <argp.h>

#include "../print.h"
#include "../partition.h"
#include "../errorestimator/refine.h"

using namespace std;

namespace HpIdeal {

class Options {
 public:
  int initial_degree = 2;
  scalar theta = 0.8;
  string bases =        "../../CombinedBases/degree20/";
  string meshfile =     "../../Meshes/lshaped6.mesh";
  string rhsfile =      "../../Meshes/lshaped6_ones.rhs";
  string paramstring = "";

  string rhs() {
    auto slash = rhsfile.find_last_of("/")+1;
    auto dot = rhsfile.find_last_of(".");
    return rhsfile.substr(slash, dot-slash-1);
  }

  string outfile() {
    if (paramstring.size() == 0) {
      return rhs() + "_ideal_default.out";
    }
    return rhs() + "_ideal_" + this->paramstring + ".out";
  }

  static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    Options *arguments = (Options *)state->input;
    switch (key)
    {
      case 't':
        arguments->theta = atof(arg);
        arguments->paramstring = arguments->paramstring + "_t" + arg;
        break;
      case 'i':
        arguments->initial_degree = atoi(arg);
        break;
      case 'b':
        arguments->bases = arg;
        break;
      case 'm':
        arguments->meshfile = arg;
        break;
      case 'r':
        arguments->rhsfile = arg;
        break;
      default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
  }
};

}

static struct argp_option options[] = {
  {"theta",     't', "(0..1]",      0, "Value of theta"},
  {"initdeg",   'i', "natural",     0, "Initial degree"},
  {"basesdir",  'b', "DIR",         0, "Bases directory"},
  {"meshfile",  'm', "FILE",        0, "File with initial h-triangulation"},
  {"rhsfile",   'r', "FILE",        0, "File with forcing function on `meshfile`"},
  { 0 }
};

static struct argp argp = { options, HpIdeal::Options::parse_opt, "hay", "heyu" };

bool touches_corner(Element* elt) {
  for (auto v : elt->verts()) {
    if (v->x == 0.0 && v->y == 0.0) return true;
    if (v->x == 1.0 && v->y == 1.0) return true;
    if (v->x == -1.0 && v->y == 1.0) return true;
    if (v->x == 1.0 && v->y == -1.0) return true;
    if (v->x == -1.0 && v->y == -1.0) return true;
    if (v->x == 0.0 && v->y == -1.0) return true;
    if (v->x == 1.0 && v->y == 0.0) return true;
  }
  return false;
}

int main(int argc, char **argv) {
  HpIdeal::Options options;
  argp_parse(&argp, argc, argv, 0, 0, &options);

  const int refine_error_est_h = 2;
  const int refine_error_est_p = 0;
  using RefineEstimator = ErrorEstimator::Refine;

  string outfile = "output/" + options.outfile();// lshaped_" + to_string(time(nullptr)) + ".out";

  ofstream outstream(outfile, ofstream::app);
  outstream << "Problem settings:" << endl;
  outstream << "mesh file: " << options.meshfile << endl;
  outstream << "rhs file: " << options.rhsfile << endl;
  outstream << endl;

  Partition p(options.bases, options.meshfile, options.rhsfile);

  p.handler().increaseTo(p.leaves(), Degree::degreeToDim(options.initial_degree));

  outstream << "Initial triangulation settings:" << endl;
  outstream << "initial degree: " << options.initial_degree << endl;
  outstream << "number of triangles: " << p.leaves().size() << endl;
  outstream << endl;

  //p.handler().redetermineOn(p.leaves());
  p.solve();

  RefineEstimator initialError(p, p.sol());
  outstream << "Error estimator settings:" << endl;
  outstream << "h-refines: " << refine_error_est_h << endl;
  outstream << "p-refines: " << refine_error_est_p << endl;

  scalar epsilon = initialError.error();
  scalar fineps = epsilon/5000000000;
  outstream << "initial error: " << epsilon << endl;
  outstream << endl;

  outstream << "hp-IDEAL settings:" << endl;
  outstream << "theta: " << options.theta << endl;
  outstream << "--------------" << endl;
  outstream.flush();

  int i = 0;
  while (true) {
    p.handler().redetermineOn(p.leaves());
    p.solve();

    RefineEstimator err(p, p.sol());
    epsilon = err.error();
    Errors errors = err.sqerrors();
    int numdofs = p.handler().recomputeNumDOFs();
    cerr << numdofs << " " << epsilon << endl;

    string solfile = Print::formatted("output/testideal_%d_%d.sol", i++, numdofs);
    cerr << "\t sol: " << solfile << endl;
    p.printSolution(solfile);

    // compute the local error on each element and store in a vector
    std::vector<std::pair<scalar, Element *>> errvec;
    for( auto const p : errors.on(p.leaves())) {
      errvec.push_back(make_pair(p.second, p.first));
    }

    // sort this vector by error in descending order
    std::sort(errvec.begin(), errvec.end());
    std::reverse(errvec.begin(), errvec.end());

    if (epsilon < fineps) break;

    // else, we refine a subset (Dorfler marking); beware of squared errors
    ElementSet marked;
    scalar untilnow = 0;
    for( auto const p : errvec) {
      untilnow += p.first;
      marked.insert(p.second);
      if (untilnow >= epsilon * options.theta * options.theta * epsilon) break;
    }

    ElementSet to_refine;
    ElementSet to_enrich;
    for (auto elt : marked) {
      // if touches a corner: h-refine
      if (touches_corner(elt)) {
        to_refine.insert(elt);
      // else: p-enrich
      } else {
        to_enrich.insert(elt);
      }
    }
    p.refine(to_refine);
    p.handler().increaseDegreeBy(to_enrich, 1);
  }
  return 0;
}
