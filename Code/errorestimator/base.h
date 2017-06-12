#pragma once
#include "../partition.h"
#include "../fem/solution.h"
#include "../errors.h"

namespace ErrorEstimator {
class Base {
protected:
  Partition &partition;
  const FEM::Solution &sol;
  Errors _sqerrors;

  Base(Partition &partition, const FEM::Solution &sol) : 
    partition(partition), sol(sol) {}

public:
  virtual Errors sqerrors() { return _sqerrors; }
  virtual scalar error() {
    auto scalarset = _sqerrors.on(partition.roots());
    scalar ret = 0;
    for(auto const p : scalarset) ret += p.second;
    return std::sqrt(ret);
  }
};
}
