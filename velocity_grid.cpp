#include "velocity_grid.h"
#include <iostream>

vector<double> Make_1D_v_grid(size_t v_size, double Vmax){
	/* This function build 1D grid size of "v_size"
	 * with boundaries [-Vmax, Vmax].
	 * Attention!!! Grid includes point v = 0. So v_size must be odd!
	 * function throws exception "invalid_argument"
	 */
	if(v_size % 2 == 0){
		throw invalid_argument("v_size isn't odd");
	}
	vector<double> result;
	result.reserve(v_size);
	for(size_t i = 0; i < v_size; ++i){
		result.push_back( (-static_cast<double>(v_size-1)*0.5 + i) * 2.0 * Vmax / (v_size - 1) );
	}
	return result;
}

VelocityGrid::VelocityGrid(size_t v_size, double Vmax_) :
				v_grid(Make_1D_v_grid(v_size, Vmax_)){}

VelocityGrid::VelocityGrid(const vector<double>& v) : v_grid(v) {}

size_t VelocityGrid::GetSize() const{
	return v_grid.size();
}

double VelocityGrid::GetGridStep() const{
	return v_grid[1] - v_grid[0];
}

vector<double> VelocityGrid::Get1DGrid() const{
	return v_grid;
}

double VelocityGrid::GetMax() const{
	return *v_grid.rbegin();
}
