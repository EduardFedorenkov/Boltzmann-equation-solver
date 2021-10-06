#pragma once

#include "space_grid.h"
#include "velocity_grid.h"
#include "Particles.h"
#include <array>
#include <set>
#include <armadillo>

using namespace std;
using namespace arma;


enum class DistributionType{
	Maxwell,              // f ~ exp(-(v - u)^2); v, u - type of vector<double>(3)
	TestDistribution_1,   // f ~ exp(-(v - u1)^2)*u2^2; v, u1, u2 - type of vector<double>(3)
};

class DistributionFunction{
public:
	DistributionFunction(const DistributionType distribution_type,
			const Particle& p, const SpaceGrid& x,
			const VelocityGrid& v,
			const vector<double>& density, const vector<double>& temperature,
			const vector<vec3>& mean_vel);

	DistributionFunction(const DistributionType distribution_type,
			const Particle& p, const SpaceGrid& x,
			const VelocityGrid& v,
			vector<double>& density, vector<double>& temperature);

	DistributionFunction(const field<cube>& df, const Particle& p, const SpaceGrid& x,
			const VelocityGrid& v);

	DistributionFunction(const cube& df, const Particle& p, const VelocityGrid& v);

	DistributionFunction(const DistributionType distribution_type, const Particle& p, const VelocityGrid& v,
			const double density, const double temperature);

	VelocityGrid GetVelGrid() const;

	SpaceGrid GetSpaceGrid() const;

	vector<double> ComputeDensity() const;

	double ComputeFallingDensity(bool Is_left_wall) const;

	vector<double> ComputeTemperature(const vector<vec3>& first_moment) const;

	vector<double> ComputeTemperatureStatic() const;

	vector<vec3> ComputeMeanVelocity() const;

	field<cube> GetFullDistribution() const;

	cube GetDistrSlice(const size_t index) const;

	Particle GetParticle() const;

	double GetOneCollisionTime() const;

	void SaveMatrix_x_vx(const size_t vy_position, const size_t vz_position, const size_t time_index) const;

	void SaveMatrixes(const set<size_t>& xy_collection, const size_t z_position,
			const size_t space_position, const size_t time_step_num) const;

	void ChangeDFSliceValues(const cube& df_slice_val, size_t slice_index);

	void ChangeDFbyTransport(size_t slice_index, double time_step);

	void RungeKutta2_ElasticCollisons(const mat& coll_mat, double time_step, size_t slice_index, bool Do_treatment);

	// ----Can be in private section----

	cube Maxwell(double density, double temperature) const;

	cube TestDistribution_1(double density, double temperature) const;

	cube MaxwellReflectedHalf(double density, double temperature, bool Is_left_wall) const;

	cube ConstTemperatureWallReflection(bool Is_left_wall) const;

	cube WallPerfectReflection(bool Is_left_wall) const;

	cube FluxbyThreePoints(const cube& df_left, const cube& df_mid, const cube& df_right) const;

	cube ComputeFlax(size_t slice_index) const;

private:
	const Particle particle;            // particle data (mass, cross section)
	SpaceGrid space_grid;                // vector size of size_s
	VelocityGrid velocity_grid;          // vector size of size_v
	field<cube> distribution_function;   /* 3D DF preserved in 3D representation
	                                      *size = [size_s X size_v X size_v X size_v].
	                                      */
};
