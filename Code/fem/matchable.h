#pragma once
#include "refinable.h"

class Matchable : public Refinable {
protected:
  void match(Element *elt);
  int _isMatching = 0;          //-1: no, 0: not initialized, 1: yes
  int _isTriTypesCorrect = 0;   //-1: no, 0: not initialized, 1: yes

  void resetRoots() override;

  /**
   *  Make constructor protected so this object does not get instantiated
   */
  Matchable(std::string basisdir) : Refinable(std::move(basisdir)) {}
public:
  // determines if we satisfy the matching condition
  bool isMatching(bool force = false);    

  // determines if the triangle types are correct
  bool isTriTypesCorrect(bool force = false);

  // construct a conforming refinement satisfying the matching condition
  void makeMatching();  

  // determines the triangle types for this partition
  void makeTriTypesCorrect();
};
