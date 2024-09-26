#pragma once
#ifndef Generator_h
#define Generator_h

#include "Point.h"
#include <vector>

std::vector<Com_Methods::Point> grid_generator(double start, double end, int parts, double sparse = 1.0);

#endif