#include "CellCollection.hpp"

long neighbor(double x_, double y_, double z_, double sidex_, double sidey_, double sidez_,
              double cell_sidex_, double cell_sidey_, double cell_sidez_, 
              int cells_per_dimensionx_, int cells_per_dimensiony_,int cells_per_dimensionz_ ) {
  if (x_ > sidex_) x_ -= sidex_;
  if (x_ < 0.0)   x_ += sidex_;
  if (y_ > sidey_) y_ -= sidey_;
  if (y_ < 0.0)   y_ += sidey_;
  if (z_ > sidez_) z_ -= sidez_;
  if (z_ < 0.0)   z_ += sidez_;
  return int(x_ / cell_sidex_) +
         int(y_ / cell_sidey_) * cells_per_dimensionx_ +
         int(z_ / cell_sidez_) * cells_per_dimensionx_ * cells_per_dimensiony_;
}

void CellCollection::init(double sidex_, double sidey_, double sidez_, double rc_, double rs_, int verlet_) {
  this->clear();

  if (verlet_ < 2) cells_per_dimensionx = 1;
  else cells_per_dimensionx = int(sidex_ / (rc_ + rs_));
  cell_sidex = sidex_ / cells_per_dimensionx;
  if (verlet_ < 2) cells_per_dimensiony = 1;
  else cells_per_dimensiony = int(sidey_ / (rc_ + rs_));
  cell_sidey = sidey_ / cells_per_dimensiony;
  if (verlet_ < 2) cells_per_dimensionz = 1;
  else cells_per_dimensionz = int(sidez_ / (rc_ + rs_));
  cell_sidez = sidez_ / cells_per_dimensionz;
  /*std::cout << "Cells: " << cells_per_dimension << " x " <<
                            cells_per_dimension << " x " <<
                            cells_per_dimension << std::endl; */
  // initialize CellCollection with the right number of elements
  int total_cells = cells_per_dimensionx * cells_per_dimensiony * cells_per_dimensionz;
  for (long k = 0; k < total_cells; k++) {
    Cell c;
    this->push_back(c);
  }

  // apply stencil to each cell to get neighbors
  for (int i = 0; i < cells_per_dimensionx; i++) {
    for (int j = 0; j < cells_per_dimensiony; j++) {
      for (int k = 0; k < cells_per_dimensionz; k++) {
	double x = 0.5 * cell_sidex + cell_sidex * i;
	double y = 0.5 * cell_sidey + cell_sidey * j;
	double z = 0.5 * cell_sidez + cell_sidez * k;
	Cell c;
	c.neighbor_ids[0]  = neighbor(x + cell_sidex, y            , z            , sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[1]  = neighbor(x            , y + cell_sidey, z            , sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[2]  = neighbor(x + cell_sidex, y + cell_sidey, z            , sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[3]  = neighbor(x + cell_sidex, y - cell_sidey, z            , sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[4]  = neighbor(x + cell_sidex, y            , z - cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[5]  = neighbor(x            , y + cell_sidey, z - cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[6]  = neighbor(x + cell_sidex, y + cell_sidey, z - cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[7]  = neighbor(x + cell_sidex, y - cell_sidey, z - cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[8]  = neighbor(x + cell_sidex, y            , z + cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[9]  = neighbor(x            , y + cell_sidey, z + cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[10] = neighbor(x + cell_sidex, y + cell_sidey, z + cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[11] = neighbor(x + cell_sidex, y - cell_sidey, z + cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
	c.neighbor_ids[12] = neighbor(x            , y            , z + cell_sidez, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
        long cell_id = neighbor(x, y, z, sidex_, sidey_, sidez_, cell_sidex, cell_sidey, cell_sidez, cells_per_dimensionx, cells_per_dimensiony, cells_per_dimensionz);
        this->at(cell_id) = c;
      }
    }
  }
}
