#include <ginac/ginac.h>
#include <cassert>
#include <fstream>
#include "../../Precomputations/triangle.h"
#include "../../Precomputations/tritype.h"
#include "../../Precomputations/monomialbasis.h"
#include "../../Precomputations/hierarchicalbasis.h"
#include "../../Precomputations/ginacsupportforeigen.h"
#include "../../Precomputations/saver.h"
#include "../matrix.h"
#include "../degree.h"

using namespace GiNaC;
using namespace std;

class Polynomial {
protected:
  symbol _x;
  symbol _y;
  ex _ex;

public:
  Polynomial(const string &fn) : _x("x"), _y("y") {
    ifstream file(fn);
    assert(file.good());

    string poly;
    file >> poly;
    _ex = ex(poly, lst{_x, _y}).expand();
    cout << poly << " " << _ex << endl;

    file.close();
  }

  void laplace() { 
    _ex = (-_ex.diff(_x, 2) - _ex.diff(_y, 2)).expand();
  }

  const ex &poly() const { return _ex; }

  int deg() const {
    int d = 0;
    for( int j = 0; j <= _ex.degree(_x); j++) {
      for( int k = 0; k <= _ex.degree(_y); k++) {
        if( Eigen::ex2double(_ex.coeff(_x,j).coeff(_y,k)) != 0.0) {
          d = max(d, j+k);
        }
      }
    }
    return d;
  }

  Eigen::VectorXd monomial_coeffs() const {
    Eigen::VectorXd vec = Eigen::VectorXd::Zero(Degree::degreeToDim(deg()));
    int i = 0;
    for( int j = 0; j <= deg(); j++) {
      for( int k = 0; k <= deg(); k++) {
        if( j + k > deg()) continue;
        vec(i) = Eigen::ex2double(_ex.coeff(_x,j).coeff(_y,k));
        i += 1;
      }
    }

    return vec;
  }
};

class Elt {

};

class Mesh {
  vector<vector<numeric>> _v;
  vector<Precomputations::Triangle> _t;
  vector<Precomputations::Tritype> _tt;
public:
  Mesh(string fn) {
    ifstream file(fn);
    assert(file.good());

    int nVerts, nElements;
    file >> nVerts;
    for( auto i = 0; i < nVerts; i++) {
      double dx, dy;
      file >> dx >> dy;
      numeric x = dx+0.0, y = dy+0.0;
      _v.push_back({x, y});
    }

    file >> nElements;
    for( auto i = 0; i < nElements; i++) {
      int i0, i1, i2, type;
      file >> i0 >> i1 >> i2 >> type;
      _t.push_back(Precomputations::Triangle(
            vector<vector<numeric>>{_v[i0], _v[i1], _v[i2]}));
      _tt.push_back(type);
    }

    file.close();
  }

  Precomputations::Triangle tri(int i) const { return _t[i]; }
  Precomputations::Tritype type(int i) const { return _tt[i]; }
  vector<numeric> v(int i) const { return _v[i]; }
  int ntris() const { return _t.size(); }
  int nverts() const { return _v.size(); }
};

int main(int argc, char** argv) {
  int filedeg = 4;
  if (argc > 1) {
    filedeg = atoi(argv[1]);
  }
  Polynomial sol("../Meshes/unitsquare_degree" + to_string(filedeg) + ".sol");
  Mesh mesh("../Meshes/unitsquare.mesh");

  sol.laplace();
  int deg = sol.deg();

  auto p = Precomputations::Triangle().points(deg);
  auto mb = Precomputations::MonomialBasis(deg);

  Eigen::VectorXd mc = sol.monomial_coeffs();

  Eigen::MatrixXd endmat(mesh.ntris(), mc.rows());
  cout << endmat.rows() << " " << endmat.cols() << endl;

  for( int i = 0; i < mesh.ntris(); i++) {
    auto refbasis = Precomputations::HierarchicalBasis(deg, mesh.type(i));
    auto locbasis = refbasis;
    locbasis.transformTo(mesh.tri(i));

    auto mb2ref = mb.transferMatrix(refbasis);
    auto ref2loc = refbasis.transferMatrix(locbasis);
    auto between = mb2ref * mc;
    auto res = ref2loc * between;
    endmat.row(i) = res.transpose();
  }
  cout << endmat << endl;
  Precomputations::Saver::save("../Meshes/unitsquare_degree" + to_string(filedeg) + ".rhs", endmat);
  return 0;
}
