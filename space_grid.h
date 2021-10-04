#pragma once

#include <vector>
using namespace std;

vector<double> Make_1D_x_grid(size_t v_size, double distance);

enum class BC_Type{
	PerfectReflection,
	ConstantTemperatureWall,
	empty
};

struct Walls{
	Walls(const pair<BC_Type, BC_Type>& walls_BC_Type, const pair<double, double>& walls_temperature);

	pair<BC_Type, BC_Type> walls_BC;
	pair<double, double> walls_T;
};

class SpaceGrid{
public:
	SpaceGrid(size_t x_size, double distance, const pair<BC_Type, BC_Type>& walls_BC_Type, const pair<double, double>& walls_temperature);

	SpaceGrid();

	size_t GetSize() const;

	double GetGridStep() const;

	vector<double> Get1DGrid() const;

	double GetDistance() const;

	Walls GetWalls() const;

private:
	const vector<double> x_grid;
	Walls walls;
};
