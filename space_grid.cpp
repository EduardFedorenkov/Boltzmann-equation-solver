#include "space_grid.h"

vector<double> Make_1D_x_grid(size_t x_size, double distance){
	/* This function build 1D grid size of "v_size"
	 * with boundaries [-Vmax, Vmax].
	 * Attention!!! Grid includes point v = 0. So v_size must be odd!
	 * function throws exception "invalid_argument"
	 */
	vector<double> result;
	result.reserve(x_size);
	for(size_t i = 0; i < x_size; ++i){
		result.push_back( (0.5 + i) * distance / x_size );
	}
	return result;
}

Walls::Walls(const pair<BC_Type, BC_Type>& walls_BC_Type, const pair<double, double>& walls_temperature) : walls_BC(walls_BC_Type),
			walls_T(walls_temperature){};

SpaceGrid::SpaceGrid(size_t x_size, double distance, const pair<BC_Type, BC_Type>& walls_BC_Type, const pair<double, double>& walls_temperature) :
				x_grid(Make_1D_x_grid(x_size, distance)), walls( Walls(walls_BC_Type, walls_temperature) ){}

SpaceGrid::SpaceGrid() : x_grid(vector<double>({0})),
		walls(Walls(make_pair(BC_Type::empty, BC_Type::empty), make_pair(0.0,0.0))) {}

size_t SpaceGrid::GetSize() const{
	return x_grid.size();
}

double SpaceGrid::GetGridStep() const{
	return x_grid[1] - x_grid[0];
}

vector<double> SpaceGrid::Get1DGrid() const{
	return x_grid;
}

double SpaceGrid::GetDistance() const{
	return *x_grid.rbegin();
}

Walls SpaceGrid::GetWalls() const{
	return walls;
}
