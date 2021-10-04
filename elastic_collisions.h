#pragma once

#include "velocity_grid.h"
#include "particles.h"
#include <armadillo>

using namespace std;
using namespace arma;

vec3 PolarAxisRotation(const vec3& old_axis, const vec3& new_axis, vec3& point);

vector<vec3> ScatteringSphere(size_t N, const vec3& s_pole, const vec3& n_pole);

pair<vec3, vec3> PostCollisionVelocities(const vec3& initial_vel,
		const vec3& direction, const double reletive_velocity);

pair<array<size_t, 3>, bool> FindNearestNode(const vec& vel_grid, const vec3& point);

mat BuildCollisionMatrix(const VelocityGrid& v, const Particle& p);

bool IndexCheckRange(const array<int, 3>& idx, size_t size);

cube DFShiftProcedure(const cube& df, const uvec3& position);

cube ComputeCollisionsIntegral(const cube& df, const mat& coll_mat, const VelocityGrid& v, bool Do_corrections);

mat55 TreatedConservationLawsMatrix(const cube& df, const VelocityGrid& v);

vec5 ComputeCollIntegralMoments(const cube& st, const VelocityGrid& v);

void TreatCollIntegral(cube& st, const cube& df, const VelocityGrid& v);
