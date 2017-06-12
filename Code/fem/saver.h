#pragma once
#include <iostream>
#include <fstream>

class Partition;

class Saver {
public:
  static void save(std::string filename, std::ostream &( *printer)(std::ostream &));
};
