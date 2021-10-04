#include "elastic_collisions.h"

vec3 PolarAxisRotation(const vec3& old_axis, const vec3& new_axis, vec3& point){
	vec3 ax = normalise(cross(old_axis, new_axis));
	if(norm(ax) >= datum::eps){
		double cos_theta = norm_dot(new_axis, old_axis);
		double theta = acos(cos_theta);
		mat33 M;
		M(0,0) = cos_theta + (1.0 - cos_theta)*Sqr(ax[0]);
		M(0,1) = (1.0 - cos_theta)*ax[0]*ax[1] - sin(theta)*ax[2];
		M(0,2) = (1.0 - cos_theta)*ax[0]*ax[2] + sin(theta)*ax[1];
		M(1,0) = (1.0 - cos_theta)*ax[0]*ax[1] + sin(theta)*ax[2];
		M(1,1) = cos_theta + (1.0 - cos_theta)*Sqr(ax[1]);
		M(1,2) = (1.0 - cos_theta)*ax[1]*ax[2] - sin(theta)*ax[0];
		M(2,0) = (1.0 - cos_theta)*ax[0]*ax[2] - sin(theta)*ax[1];
		M(2,1) = (1.0 - cos_theta)*ax[1]*ax[2] + sin(theta)*ax[0];
		M(2,2) = cos_theta + (1.0 - cos_theta)*Sqr(ax[2]);
		point = M * point;
	}
	return point;
}

vector<vec3> ScatteringSphere(size_t N, const vec3& s_pole, const vec3& n_pole){
	vector<vec3> sphere_points;
	sphere_points.reserve(N);
	double gr = (1.0 + sqrt(5)) / 2.0;
	double ga = 2.0*datum::pi*(1.0 - 1.0 / gr);

	vec3 e_z = {0,0,1};
	vec3 spiral_axis = s_pole - n_pole;
	for(size_t i = 0; i < N; ++i){
		vec3 point;
		point(0) = sin(acos(1.0 - i * (2.0 / (N-1))))*cos(i*ga);
		point(1) = sin(acos(1.0 - i * (2.0 / (N-1))))*sin(i*ga);
		point(2) = 1.0 - i * (2.0 / (N-1));
		PolarAxisRotation(e_z, spiral_axis, point);
		sphere_points.push_back(point);
	}
	return sphere_points;
}

pair<vec3, vec3> PostCollisionVelocities(const vec3& initial_vel,
		const vec3& direction, const double reletive_velocity){
	vec3 velocity_out_1, velocity_out_2;
	velocity_out_1 = (direction*reletive_velocity + initial_vel)*0.5;
	velocity_out_2 = (-direction*reletive_velocity + initial_vel)*0.5;
	return make_pair(velocity_out_1,velocity_out_2);
}

pair<array<size_t, 3>, bool> FindNearestNode(const vec& vel_grid, const vec3& point){
	double grid_step = vel_grid(1) - vel_grid(0);
	array<size_t, 3> result = {0, 0, 0};
	for(size_t i = 0; i < 3; ++i){
		result[i] = round((point(i) - vel_grid(0)) / grid_step);
		if(result[i] < -0 || result[i] >= vel_grid.size()){
			return make_pair(result, false);
		}
	}
	return make_pair(result, true);
}

mat BuildCollisionMatrix(const VelocityGrid& v, const Particle& p){
	size_t v_size = v.GetSize();
	size_t N_angle = 500;
	double phase_volume = pow(v.GetGridStep(),3);
	vec vel_1D(v.Get1DGrid());
	mat sourse(v_size * v_size * v_size, v_size * v_size * v_size);
	sourse.fill(0);
	mat runoff(sourse);
	double factor = p.hard_speres_cross_section * phase_volume * 4 * datum::pi / N_angle;
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				vec3 second_particle_velocity = {vel_1D(m), vel_1D(l), vel_1D(k)};
				double reletive_velocity = sqrt(Sqr(second_particle_velocity(0))
						+ Sqr(second_particle_velocity(1))
						+ Sqr(second_particle_velocity(2)));
				vector<vec3> sphere_points;
				if(array<size_t,3>({m, l, k}) != array<size_t,3>({0, 0, 0})){
					sphere_points = ScatteringSphere(N_angle, second_particle_velocity, {0,0,0});
				}
				for(auto& point : sphere_points){
					pair<vec3, vec3> post_collision_velocities = PostCollisionVelocities(second_particle_velocity, point, reletive_velocity);
					auto Node_1 = FindNearestNode(vel_1D, post_collision_velocities.first);
					auto Node_2 = FindNearestNode(vel_1D, post_collision_velocities.second);
					if (Node_1.second && Node_2.second){
						sourse(v_size*v_size*Node_1.first[2] + v_size*Node_1.first[1] + Node_1.first[0],
								v_size*v_size*Node_2.first[2] + v_size*Node_2.first[1] + Node_2.first[0]) +=
								factor * reletive_velocity;
						runoff(v_size*v_size*k + v_size*l + m, (v_size * v_size * v_size - 1) / 2) +=
								factor * reletive_velocity;
					}
				}
			}
		}
	}
	return sourse - runoff;
}

bool IndexCheckRange(const array<int, 3>& idx, size_t size){
	bool result;
	for(const int item : idx){
		if(item >= 0 && static_cast<size_t>(item)< size){
			result = true;
		}else{
			return false;
		}
	}
	return result;
}

cube DFShiftProcedure(const cube& df, const uvec3& position){
	size_t v_size = df.n_slices;
	size_t mid_point = (v_size - 1) / 2;
	cube result(size(df));
	result.fill(0);
	array<int, 3> slice_index = {0,0,0};
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				slice_index[0] = m + position(0) - mid_point;
				slice_index[1] = l + position(1) - mid_point;
				slice_index[2] = k + position(2) - mid_point;
				if(IndexCheckRange(slice_index, v_size)){
					result(m,l,k) = df(slice_index[0], slice_index[1], slice_index[2]);
				}
			}
		}
	}
	return result;
}

cube ComputeCollisionsIntegral(const cube& df, const mat& coll_mat, const VelocityGrid& v, bool Do_corrections){
	size_t v_size = df.n_slices;
	cube collisions_integral(size(df));
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				cube df_shift = DFShiftProcedure(df, {m, l, k});
				collisions_integral(m,l,k) = as_scalar(vectorise(df_shift).t()*coll_mat*vectorise(df_shift));
			}
		}
	}
	if(Do_corrections){
		TreatCollIntegral(collisions_integral, df, v);
	}
	return collisions_integral;
}

mat55 TreatedConservationLawsMatrix(const cube& df, const VelocityGrid& v){
	mat::fixed<5,5> m;
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();
	size_t v_size = v.GetSize();

	m(0,0) = accu(df);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t i = 0; i < v_size; ++i){
				m(0,1) = m(1,0) += df(i,l,k) * vel_1D(i);
				m(0,2) = m(2,0) += df(i,l,k) * vel_1D(l);
				m(0,3) = m(3,0) += df(i,l,k) * vel_1D(k);
				m(1,1) += df(i,l,k) * Sqr(vel_1D(i));
				m(2,2) += df(i,l,k) * Sqr(vel_1D(l));
				m(3,3) += df(i,l,k) * Sqr(vel_1D(k));

				m(2,1) = m(1,2) += df(i,l,k) * vel_1D(i) * vel_1D(l);
				m(3,1) = m(1,3) += df(i,l,k) * vel_1D(i) * vel_1D(k);
				m(3,2) = m(2,3) += df(i,l,k) * vel_1D(l) * vel_1D(k);

				m(0,4) = m(4,0) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k);
				m(1,4) = m(4,1) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(i);
				m(2,4) = m(4,2) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(l);
				m(3,4) = m(4,3) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(k);
				m(4,4) += Sqr( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k);
			}
		}
	}
	//m.for_each([](double& x){ if(abs(x) < datum::eps){x = 0;} });
	return m;
}

vec5 ComputeCollIntegralMoments(const cube& st, const VelocityGrid& v){
	vec5 st_moments;
	size_t v_size = v.GetSize();
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();

	st_moments(0) = accu(st);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				st_moments(1) += st(m,l,k) * vel_1D(m);
				st_moments(2) += st(m,l,k) * vel_1D(l);
				st_moments(3) += st(m,l,k) * vel_1D(k);
				st_moments(4) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(m)) ) * st(m,l,k);
			}
		}
	}
	//st_moments.for_each([](double& x){if(abs(x) < datum::eps){ x = 0;}});
	return -st_moments;
}

void TreatCollIntegral(cube& st, const cube& df, const VelocityGrid& v){
	vec coef = solve(TreatedConservationLawsMatrix(df,v),ComputeCollIntegralMoments(st,v));
	size_t v_size = v.GetSize();
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				st(m,l,k) += ( coef(0) + coef(1) * vel_1D(m) + coef(2) * vel_1D(l) + coef(3) * vel_1D(k)
						+ coef(4) * ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(m)) ) ) * df(m,l,k);
			}
		}
	}
}
