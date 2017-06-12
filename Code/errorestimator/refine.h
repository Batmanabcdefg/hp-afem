#include "base.h"
#include "../fem/solver.h"
#include "../dofhandler.h"
#include "../triangleset.h"

/**
 *  RefineErrorEstimator.h
 *
 *  An error estimator based on the idea that just refining all elements a
 *  couple of times (in either h- or p-sense) will yield a solution so much
 *  better that we can view it as the "exact" solution.  Because it is fairly
 *  easy to compute the norm-difference of two piecewise polynomials (see 
 *  Approximator.h), this is easy to implement.
 *
 *  h is the number of times we will run refineLeavesUniformly, i.e., 
 *  the number of times we will multiply the number of leaves with four.
 *  p is the input of increaseDegreeBy(), i.e., the number of degrees we will 
 *  add to the elements.
 */

namespace ErrorEstimator {
class Refine : public Base {
public:
  Refine(Partition &partition, const FEM::Solution &sol)
    : Refine(partition, sol, 2, 0) {}
  Refine(Partition &partition, const FEM::Solution &sol, int h, int p)
    : ErrorEstimator::Base(partition, sol)
  {
    // make local copies of the DOFHandler and current solution
    DOFHandler curhandler = partition.handler();
    FEM::Solution cursol = sol;

    // make local copy of the current leaves so that we can restore later
    ElementSet leaves = partition.leaves();

    // do some p-refinements
    if( p > 0) curhandler.increaseDegreeBy(leaves, p);

    // do some h-refinements
    for( int i = 0; i < h; i++) {
      partition.refineLeavesUniformly(curhandler);
    }

    // solve the system on this refined partition
    FEM::Solution exact = FEM::Solver(curhandler, partition.rhs()).sol();

    // get the (squared) errors
    _sqerrors = Errors(cursol.squaredH1NormsOfDifferenceWith(exact));

    partition.printElementScalarSet(_sqerrors.on(leaves),
                                    "output/errors_" + std::to_string(partition.handler().recomputeNumDOFs()) + "_" 
                                                     + std::to_string(error()) + ".error");

    // get the partition back in its original state, at least leaves-wise
    for( auto &elt: leaves) partition.trim(elt);
  }
};
}
