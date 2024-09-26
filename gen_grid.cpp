#include "gen_grid.h"
#include "Point.h"

std::vector<Com_Methods::Point> grid_generator(double start, double end, int parts, double sparse){

	if (fabs(start - end) < 1e-7 || parts < 2 || sparse < 1e-7)
		throw std::exception("Invalid input");

	std::vector<Com_Methods::Point> res;
	double iter = start;
	res.push_back(Com_Methods::Point(iter, 0.0, 0.0));

	double temp = 1;
	for (int i = 1; i < parts; i++)
		temp += pow(sparse, i);
	double step = (end - start) / temp;

	for (int i = 0; i < parts; i++)
		res.push_back(Com_Methods::Point(iter += step * pow(sparse, i), 0.0, 0.0));

	return res;
}