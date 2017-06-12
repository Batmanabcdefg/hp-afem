/**
 * 1. Read a solution file from the command line arguments.
 * 2. Create in-memory partition and solution from this.
 * 3. Estimate the error.
 * 4. Output this.
 */

#include <iostream>
#include <Eigen/Sparse>
#include <argp.h>

#include "../spectra/include/SymEigsSolver.h"
#include "../spectra/include/MatOp/SparseSymMatProd.h"
#include "../fem/solver.h"
#include "../partition.h"
#include "../errorestimator/refine.h"

using namespace std;

namespace EstimateError {
class Options {
 public:
  string infile = "output/testLshaped_15.0_3855.sol";
  int hrefines = 2;
  int prefines = 0;
  string bases =        "../../CombinedBases/degree20/";
  string meshfile =     "../../Meshes/lshaped6.mesh";
  string rhsfile =      "../../Meshes/lshaped6_ones.rhs";
  bool analyze = false;

  string to_string() {
    return meshfile + " " + rhsfile + " " + infile + " " + std::to_string(hrefines) + " " + std::to_string(prefines);
  }

  string outfile() {
    if (analyze) {
      return "aanalyze3_" + infile.substr(infile.find_last_of("/")+1) + ".out";
    }

    if (rhsfile.find("/") != string::npos) {
      return "calibration_" + infile.substr(infile.find_last_of("/")+1) + ".out";
    }
    return "calibration_" + infile + ".out";
  }

  static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    Options *arguments = (Options *)state->input;
    switch (key)
    {
      case 'i':
        arguments->infile = arg;
        break;
      case 'h':
        arguments->hrefines = atoi(arg);
        break;
      case 'p':
        arguments->prefines = atoi(arg);
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
      case 'a':
        arguments->analyze = atoi(arg) > 0;
        break;
      default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
  }
};
}

static struct argp_option options[] = {
  {"infile",   'i', "FILE",    0, "Solution file"},
  {"hrefines", 'h', "natural", 0, "Number of h-refines in error estimator"},
  {"prefines", 'p', "natural", 0, "Number of p-refines in error estimator"},
  {"basesdir",  'b', "DIR",         0, "Bases directory"},
  {"meshfile",  'm', "FILE",        0, "File with initial h-triangulation"},
  {"rhsfile",   'r', "FILE",        0, "File with forcing function on `meshfile`"},
  {"analyze",   'a', "bool",        0, "Analyze global stiffness matrix"},
  { 0 }
};

static struct argp argp = { options, EstimateError::Options::parse_opt, "hay", "heyu" };

// Exposes the protected constructor; this is an ad-hoc partition.
class AdhocPartition : public Partition {
 public:
  AdhocPartition(std::string basisdir, std::string meshfn)
    : Partition(basisdir, meshfn)
  {}

  void setOnesRhs() {
    _rhs = FEM::Rhs(3);
    for (Element* leaf : leaves()) {
      Vector rhs = Vector::Ones(3);
      _rhs.insert_vector(leaf, rhs, true);
    }

    resetRoots();
    rebuildElementFinder();
    assert(isConform(true));
    assert(isTriTypesCorrect(true));
    handler().redetermine();
  }
};

scalar cond2spd(const StiffnessMatrix& M) {
  int itermax = 1000;
  scalar rtol = 1e-6;

  Eigen::SimplicialLDLT<StiffnessMatrix> sldlt;
  sldlt.compute(M);

  int n = M.rows();
  Vector x = Vector::Random(n);
  x = x/sqrt(x.dot(x));

  Vector y = M * x;
  scalar lM_old = x.dot(y);
  scalar lM = -1;
  int iter = 0;

  while (true) {
    iter += 1;
    x = y/sqrt(y.dot(y));
    y = M * x;
    lM = x.dot(y);
    cerr << "it " << iter << ": lM = " << lM << endl;
    auto crit = abs((lM - lM_old)/lM);
    if (crit < rtol && iter > 3) {
      break;
    } else {
      lM_old = lM;
    }
    if (iter >= itermax) {
      cerr << "no convergence: got " << crit << " instead of " << rtol << endl;
      break;
    }
  }

  auto vM = x;

  x = Vector::Random(n);
  x = x/sqrt(x.dot(x));
  y = sldlt.solve(x);
  scalar lm_old = x.dot(y);
  scalar lm = -1;
  iter = 0;

  while (true) {
    iter = iter + 1;
    x = y/sqrt(y.dot(y));
    y = sldlt.solve(x);
    lm = x.dot(y);
    cerr << "it " << iter << ": lm = " << 1.0/lm << endl;
    auto crit = abs((lm - lm_old)/lm);
    if (crit < rtol && iter > 3) {
      break;
    } else {
      lm_old = lm;
    }
    if (iter >= itermax) {
      cerr << "no convergence: got " << crit << " instead of " << rtol << endl;
      break;
    }
  }
  auto vm = x;
  lm = 1.0/lm;
  return lM/lm;
}


scalar cond2spectra(const StiffnessMatrix& A) {
  StiffnessMatrix Asq = A;
  Spectra::SparseSymMatProd<scalar> op(A);
  Spectra::SymEigsSolver<scalar, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<scalar>> eigs(&op, 1, A.rows());
  eigs.init();
  eigs.compute();
  scalar maxev = eigs.eigenvalues().real().maxCoeff();

  Spectra::SymEigsSolver<scalar, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<scalar>> eigs2(&op, 1, A.rows());
  eigs2.init();
  eigs2.compute();
  scalar minev = eigs2.eigenvalues().real().minCoeff();
  return abs(maxev)/abs(minev);
}

int main(int argc, char **argv) {
  EstimateError::Options options;
  argp_parse(&argp, argc, argv, 0, 0, &options);

  string outfile = "output/" + options.outfile();
  cerr << outfile << endl;

  // Load hp-mesh into partition
  AdhocPartition rp(options.bases, options.infile);
  if (options.rhsfile.find("_ones") != string::npos) {
    rp.setOnesRhs();
    if (options.analyze) {
      std::unique_ptr<FEM::Solver> solver = rp.solve();
      auto residual = solver->systemMatrix() * solver->sol().global() - solver->systemRhs();
      scalar residual_norm = sqrt(residual.dot(residual));

      ofstream outtstream(outfile, ofstream::app);
      //outtstream << solver->systemMatrix() << endl << endl;
      outtstream << "8 " << solver->systemMatrix().nonZeros() << endl;
      outtstream << "7 " << residual.size() << " " << residual_norm << endl;
      outtstream << "9 " << cond2spectra(solver->systemMatrix()) << endl;



      outtstream.flush();
    } else {
      rp.setSolution(FEM::Solution::fromSolFile(options.infile, rp.handler().asSet(), rp.rhs()));

      ErrorEstimator::Refine estimated_error(rp, rp.sol(), options.hrefines, options.prefines);
      ofstream outstream(outfile, ofstream::app);
      cerr << options.to_string() << " " << estimated_error.error() << endl;
      outstream << options.to_string() << " " << estimated_error.error() << endl;
      outstream.close();
    }
  }

  return 0;
}
