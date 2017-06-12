#include <ginac/ginac.h>
//#include <omp.h>
#include <iostream>
#include <fstream>
#include "hierarchicalbasis.h"
#include "saver.h"
#include "tritype.h"

using namespace std;
using namespace GiNaC;
using namespace Precomputations;

int main(int argc, char **argv) {
  size_t starttt = 0;
  size_t ncores = 8;
  size_t startdeg = 3;
  size_t enddeg = 4;
  if( argc >= 3) {
    starttt = atoi(argv[1]);
    ncores = atoi(argv[2]);
  }
  if( argc >= 5) {
    startdeg = atoi(argv[3]);
    enddeg = atoi(argv[4]);
  }
  Triangle ref;
  numeric one = 1;
  numeric zero = 0;
  vector<numeric> half = {one/(one+one), one/(one+one)};
  vector<numeric> origin = {zero, zero};
  vector<numeric> right = {one, zero};
  vector<numeric> top = {zero, one};
  vector<numeric> bot = {zero, -one};
  vector<vector<numeric> > vertsLeft = {half, origin, right};
  vector<vector<numeric> > vertsRight = {half, top, origin};
  vector<vector<numeric> > vertsK2 = {origin, bot, right};
  Triangle triLeft(vertsLeft);
  Triangle triRight(vertsRight);
  Triangle triK2(vertsK2);

  string root = "Bases2/";

  for( size_t deg = startdeg; deg < enddeg; deg++) {
    for( size_t t = starttt; t < 8; t+= ncores) {
      Tritype tt(t);
      Basis b = HierarchicalBasis(deg, tt);
      
      /*
      Saver::save(root + "degree" + to_string(deg) +
                  "/basis_" + to_string(t) + ".mat", b);

      Basis cLeft = HierarchicalBasis(deg, tt.left());
      cLeft.transformTo(triLeft);
      Basis cRight = HierarchicalBasis(deg, tt.right());
      cRight.transformTo(triRight);
      */
      Basis cK2 = HierarchicalBasis(deg, tt);
      cK2.transformTo(triK2);

      /*
      Saver::save(root + "degree" + to_string(deg) + 
                  "/eltmat0_" + to_string(t) + ".mat", b.A0());
      Saver::save(root + "degree" + to_string(deg) + 
                  "/eltmat1_" + to_string(t) + ".mat", b.A1());
      Saver::save(root + "degree" + to_string(deg) + 
                  "/eltmat2_" + to_string(t) + ".mat", b.A2());
      Saver::save(root + "degree" + to_string(deg) + 
                  "/massmat_" + to_string(t) + ".mat", b.M());
      Saver::save(root + "degree" + to_string(deg) + 
                  "/eltvec_" + to_string(t) + ".mat", b.b());
                  */
      /*
      Saver::save(root + "degree" + to_string(deg) + 
                  "/transfermat0_" + to_string(t) + 
                  ".mat", b.transferMatrixs(cLeft));
      Saver::save(root + "degree" + to_string(deg) + 
                  "/transfermat1_" + to_string(t) + 
                  ".mat", b.transferMatrixs(cRight));
                  */
      /*
      Saver::saveSymbolic(root + "degree" + to_string(deg) + 
                  "/FPhi_k0_" + to_string(t) + 
                  ".symmat", b.genFPhi_ks(cLeft));
      Saver::saveSymbolic(root + "degree" + to_string(deg) + 
                  "/FPhi_k1_" + to_string(t) + 
                  ".symmat", b.genFPhi_ks(cRight));
      Saver::saveSymbolic(root + "degree" + to_string(deg) + 
                  "/FPhi0_" + to_string(t) + 
                  ".symmat", b.genFPhis(cLeft));
      Saver::saveSymbolic(root + "degree" + to_string(deg) + 
                  "/FPhi1_" + to_string(t) + 
                  ".symmat", b.genFPhis(cRight));
                  */
    }
  }
  return 0;
}
