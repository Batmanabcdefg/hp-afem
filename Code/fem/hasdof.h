#pragma once

#include "degree.h"

class HasDOF {
public:
  virtual int dof()     const = 0;
          int degree()  const { return Degree::dofToDegree(dof()); }
          int dim()     const { return Degree::dofToDim(dof()); }
};
