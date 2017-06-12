#pragma once

#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>

#include "errorestimator/base.h"
#include "fem/solution.h"
#include "partition.h"

//E is the error estimator; theta the Dorfler marking constant
template<class E>
class Reduce {
private:
  Partition &_partition;
  FEM::Solution &_sol;
  scalar _delta;
  scalar _theta;
  scalar _error = -1.;
  static int _i;

public:
  Reduce(Partition &partition, FEM::Solution &sol, scalar delta, scalar theta) :
    Reduce(partition, sol, delta, theta, false) {}

  Reduce(Partition &partition, FEM::Solution &sol, scalar delta, scalar theta, bool print_rhs) :
    Reduce(partition, sol, delta, theta,
        std::ofstream("output/testLshaped_" + std::to_string(theta) + "_" + scalarname + ".errlog", std::ios_base::app))
  {}

  Reduce(Partition &partition, FEM::Solution &sol, scalar delta, scalar theta, bool print_rhs, std::ostream& os) :
    _partition(partition), _sol(sol), _delta(delta), _theta(theta) 
  {
    using namespace std;
    cerr << "In Reduce with delta=" << delta << " and theta=" << theta << endl;
    _i++;
    int j = 0;
    while(true) {
      _partition.handler().redetermineOn(_partition.leaves());
      _partition.solve();

      // instantiate the error estimator
      E err = E(_partition, _sol);
      _error = err.error();

      // compute the local error on each element and store in a vector
      Errors errors = err.sqerrors();
      std::vector<std::pair<scalar, Element *>> errvec;
      for( auto const p : errors.on(_partition.leaves())) {
        errvec.push_back(make_pair(p.second, p.first));
      }

      // sort this vector by error in descending order
      std::sort(errvec.begin(), errvec.end());
      std::reverse(errvec.begin(), errvec.end());

      // produce some logs
      cerr << "\tError was: " << _error << "; needed " << _delta << endl;
      int numdofs = _partition.handler().recomputeNumDOFs();
      os << "0 " << time(nullptr) << ": " << _i << " " << j << " " << numdofs << " " << _error << endl;
      string solfile = Print::formatted("output/testLshaped_%d.%d_%d.sol", _i, j++, numdofs);
      cerr << "\t sol: " << solfile << endl;
      _partition.printSolution(solfile);
      if (print_rhs) {
        string rhsfile = Print::formatted("output/testLshaped_%d.%d_%d.sol", _i, j-1, numdofs);
        ofstream os(rhsfile, ofstream::trunc);
        _partition.rhs().printForFile(_partition.leaves(), os);
        os.flush();
        os.close();
      }

      // if we are small enough, we are done
      if(_error <= delta) break;

      // else, we refine a subset (Dorfler marking); beware of squared errors
      ElementSet marked;
      scalar untilnow = 0;
      for( auto const p : errvec) {
        untilnow += p.first;
        marked.insert(p.second);
        if( untilnow >= _error * theta * theta * _error) break;
      }
      cerr << "needed " << marked.size() << " out of " << errvec.size() << " to get to theta=" << _theta << endl;

      _partition.refine(marked);
      ElementSet marked2;
      for (const auto& elt : marked) {
        // We know for sure that these are not nullpointers.
        marked2.insert(elt->left());
        marked2.insert(elt->right());
      }
      _partition.refine(marked2);
    }
  }

  scalar error() const { return _error; }
};

template<class E> int Reduce<E>::_i = 0;
