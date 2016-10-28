#ifndef _NANOBD_CELLCOLLECTION_HPP
#define _NANOBD_CELLCOLLECTION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include "Cell.hpp"

class CellCollection : public std::vector<Cell> {
  public:
    CellCollection() : std::vector<Cell>() {};
    int cells_per_dimensionx;
    int cells_per_dimensiony;
    int cells_per_dimensionz;
    double cell_sidex;
    double cell_sidey;
    double cell_sidez;
    void init(double sidex_, double sidey_, double sidez_, double rc_, double rs_, int verlet_);
};

#endif
