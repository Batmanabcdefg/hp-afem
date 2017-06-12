#pragma once
#include <iostream>
#include "basis.h"
#include "solvable.h"
#include "matchable.h"

class NearBest;

class Partition : public Matchable {
public:
  Partition(std::string basisdir, std::string meshfn, std::string rhsfn);

  friend class NearBest;

 protected:
  Partition(std::string basisdir, std::string meshfn);
};
