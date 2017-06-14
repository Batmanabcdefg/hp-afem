#include <iostream>
#include <string>
#include <argp.h>

#include "../fem/solver.h"
#include "../print.h"
#include "../partition.h"
#include "../nearbest.h"
#include "../errorestimator/refine.h"
#include "../reduce.h"

using namespace std;


namespace HpAFEM {

class Options {
 public:
  int initial_degree = 2;
  scalar theta = 0.8;
  scalar mu = 0.5;
  scalar omega = 4;
  string bases =        "../../CombinedBases/degree20/";
  string meshfile =     "../../Meshes/lshaped6.mesh";
  string rhsfile =      "../../Meshes/lshaped6_ones.rhs";
  bool hafem_do = false;
  int hafem_iter = 0;
  bool print_rhs = false;
  string paramstring = "";
  bool analyze = false;

  string rhs() {
    auto slash = rhsfile.find_last_of("/")+1;
    auto dot = rhsfile.find_last_of(".");
    return rhsfile.substr(slash, dot-slash-1);
  }

  string outfile() {
    if (paramstring.size() == 0) {
      return rhs() + "_default.out";
    }
    return rhs() + "_" + this->paramstring + ".out";
  }

  static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    Options *arguments = (Options *)state->input;
    switch (key)
    {
      case 't':
        arguments->theta = atof(arg);
        arguments->paramstring = arguments->paramstring + "_t" + arg;
        break;
      case 'u':
        arguments->mu = atof(arg);
        arguments->paramstring = arguments->paramstring + "_u" + arg;
        break;
      case 'w':
        arguments->omega = atof(arg);
        arguments->paramstring = arguments->paramstring + "_w" + arg;
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
      case 'h':
        arguments->hafem_do = atoi(arg) > 0;
        arguments->hafem_iter = atoi(arg);
        arguments->paramstring = arguments->paramstring + "_h" + arg;
        break;
      case 'a':
        arguments->analyze = atoi(arg) > 0;
        break;
      case 'p':
        arguments->print_rhs = atoi(arg) > 0;
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
  {"mu",        'u', "(0,1)",       0, "Reduction parameter"},
  {"omega",     'w', "(2,infty)",   0, "Coarsening factor"},
  {"initdeg",   'i', "natural",     0, "Initial degree"},
  {"basesdir",  'b', "DIR",         0, "Bases directory"},
  {"meshfile",  'm', "FILE",        0, "File with initial h-triangulation"},
  {"rhsfile",   'r', "FILE",        0, "File with forcing function on `meshfile`"},
  {"hafem",     'h', "natural",     0, "hp-AFEM iteration above which to do h-AFEM"},
  {"printrhs",  'p', "bool",        0, "print FEM RHS to file after each iteration"},
  {"analyze",   'a', "bool",        0, "Analyze global stiffness matrix"},
  { 0 }
};

static struct argp argp = { options, HpAFEM::Options::parse_opt, "hay", "heyu" };

int main(int argc, char **argv) {
  HpAFEM::Options options;
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

  p.refineLeavesUniformly();
  p.refineLeavesUniformly();

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

  outstream << "hp-AFEM settings:" << endl;
  outstream << "theta: " << options.theta << endl;
  outstream << "omega: " << options.omega << endl;
  outstream << "mu: " << options.mu << endl;
  outstream << "--------------" << endl;
  outstream.flush();

  int i = 0;
  while( epsilon > fineps) {
    cerr << endl;

    RefineEstimator prestartError(p, p.sol());

    scalar prestart_h1_error = prestartError.error();
    outstream << "1 " << prestart_h1_error << endl;
    cerr << "Going into NearBest" << endl;
    //find a good approximation to the current numerical solution
    NearBest(p, p.sol(), options.omega*epsilon, 1000000, outstream);
    cerr << "nu hier" <<endl;
    p.printHpMesh("output/endmesh_" + to_string(++i) + ".mesh");

    //make conform
    cerr << "Making conform again" << endl;
    p.makeConform();
    // resolve for the test on #DoF a few lines below
    p.handler().redetermineOn(p.leaves());

    if (options.analyze) {
      std::unique_ptr<FEM::Solver> solver = p.solve();
      auto residual = solver->systemMatrix() * solver->sol().global() - solver->systemRhs();
      scalar residual_norm = sqrt(residual.dot(residual));

      ofstream outtstream(outfile, ofstream::app);
      outtstream << "7 " << residual_norm << endl;
      outtstream.flush();
    } else {
      p.solve();
    }

    RefineEstimator initialError(p, p.sol());

    scalar start_h1_error = initialError.error();
    outstream << "2 " << p.handler().maxDegree() << " " << start_h1_error << endl;
    
    //refine some more
    cerr << "Going to refine to get the necessary error" << endl;
    if(options.hafem_do && i >= options.hafem_iter) {
      Reduce<RefineEstimator> reduce(p, p.sol(), 0, options.theta, false, outstream);
      cerr << "final error was " << reduce.error() << endl;
    } else {
      Reduce<RefineEstimator> reduce(p, p.sol(), options.mu*epsilon,
                                     options.theta, false, outstream);
      cerr << "final error was " << reduce.error() << endl;
    }

    //find approximation to solution
    epsilon *= options.mu;
  }
  
  int Nmax = 1000000;
  cerr << endl;
  cerr << "Calling NearBest one last time" << endl;
  NearBest(p, p.sol(), options.omega*epsilon, Nmax, outstream);
  cerr << "nu hier" <<endl;

  cerr << "Conforming and recomputing DOFs" << endl;
  p.makeConform();
  p.handler().redetermineOn(p.leaves());

  cerr << "solving" << endl;
  p.solve();

  //linearize
  cerr << "Going to linearize" << endl;
  p.linearizeSolution();

  string solfile = "output/testLshaped_final_" + to_string(Nmax) + ".sol"; 
  cerr << "Saving to " << solfile << endl;
  p.printLinearInterpolant(solfile + ".lin");
  p.printSolution(solfile);
}
