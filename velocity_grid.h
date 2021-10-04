#pragma once

#include <vector>
using namespace std;

vector<double> Make_1D_v_grid(size_t v_size, double Vmax);

class VelocityGrid{
public:
	VelocityGrid(size_t v_size, double Vmax_);

	VelocityGrid(const vector<double>& v);

	size_t GetSize() const;

	double GetGridStep() const;

	vector<double> Get1DGrid() const;

	double GetMax() const;

private:
	const vector<double> v_grid;
};
