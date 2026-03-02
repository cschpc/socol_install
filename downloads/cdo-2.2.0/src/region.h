/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REGION_H
#define REGION_H

struct Regions
{
  std::vector<double> x, y;
  std::vector<size_t> segmentSize;
  std::vector<size_t> segmentOffset;
  size_t numSegments = 0;
};

void read_regions_from_file(const char *filename, Regions &regions);
void read_regions_from_dcw(const char *codeNames, Regions &regions);

#endif  //  REGION_H
