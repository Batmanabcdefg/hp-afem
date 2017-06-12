#include "saver.h"

using namespace std;

class Partition;

void Saver::save(string filename, ostream &( *printer)(ostream &)) {
  ofstream os(filename, ofstream::trunc);
  printer(os);
  os.flush();
  os.close();
}
