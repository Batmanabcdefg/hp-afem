/**
 *  Solvable.h
 *
 *  This class is the heart of the FEM-part of this program.  It handles 
 *  everything pertaining to solving the Finite-Element approximation of the 
 *  PDE.
 *
 *  Definitions:
 *  - Transfering a solution: this moves a solution from a source set of 
 *    ancestor nodes to a target set of descendant nodes, invalidating the 
 *    solution at the source nodes.  It must be noted that the target nodeset
 *    is required to be a partition, i.e., a full binary tree.  In the process, 
 *    it is possible that degrees 
 *    of freedom are copied from parent to child nodes.  This has the nasty 
 *    side-effect that parents share degrees of freedom with their children and 
 *    should be done sparingly.  Currently, this is used when creating a linear
 *    interpolant ready for plotting.
 *  - Coping a solution: this does pretty much the same as transferring, but
 *    does not invalidate the solution at the source nodes, nor at any node
 *    between source and target.  This should be used even more sparingly,
 *    and is only used in the internals of NearBest.
 */

#pragma once
#include <map>
#include <fstream>
#include <iostream>
#include "elementtree.h"
#include "elementfinder.h"
#include "dofhandler.h"
#include "system.h"
#include "fem/solution.h"
#include "fem/solver.h"
#include "fem/rhs.h"

class Solvable : public ElementTree {
 public:
  Solvable() : ElementTree(), _handler(finder()) {}

  void setSolution(const FEM::Solution &newSol) { _sol = newSol; }
  void resetSolution() { setSolution(FEM::Solution()); }
  std::unique_ptr<FEM::Solver> solve();

  /* expose some members */
  FEM::Solution &sol() { return _sol; }
  FEM::Rhs &rhs() { return _rhs; }
  DOFHandler &handler() { return _handler; }

  /* printing various things */
  std::ostream &printLinearInterpolant(std::ostream &os = std::cout);
  void          printLinearInterpolant(std::string filename);
  std::ostream &printVerts(std::ostream &os = std::cout);
  std::ostream &printDOFs(std::ostream &os = std::cout);
  void          printDOFs(std::string filename);
  std::ostream &printSolution(std::ostream &os = std::cout);
  void          printSolution(std::string filename);
  std::ostream &printLeaves(std::ostream &os = std::cout);
  void          printLeaves(std::string filename);
  std::ostream &printRhsMesh(std::ostream &os = std::cout);
  void          printRhsMesh(std::string filename);
  std::ostream &printHpMesh(std::ostream &os = std::cout);
  void          printHpMesh(std::string filename);
  std::ostream &printElementScalarSet(const ElementScalarSet &on, std::ostream &os = std::cout);
  void          printElementScalarSet(const ElementScalarSet &on, std::string filename);

 protected:
  FEM::Solution _sol;
  FEM::Rhs _rhs;
  DOFHandler _handler;
};
