#ifndef GRID_READ_PINGO_H
#define GRID_READ_PINGO_H

#include <cstdio>
#include <vector>

size_t input_darray(FILE *gfp, size_t n_values, std::vector<double> &array);
int grid_read_pingo(FILE *gfp);

#endif
