#pragma once

#include <armadillo>

using namespace std;
using namespace arma;

template<typename T>
T Sqr(T x){
	return x*x;
}

enum class ParticlesType{
	H,
	H_2,
	Electron,
	Proton,
	H_2Plus,
	test
};

enum class CollisionModel{
	hard_spheres,
	elastic_from_database,
	excitation_from_database,
	test
};

struct Particle{
	Particle() : particle_type(ParticlesType::test), collision_model(CollisionModel::test), mass(Sqr(datum::c_0*100)), hard_speres_cross_section(1.0){};

	Particle(ParticlesType p, CollisionModel cm) : particle_type(p), collision_model(cm){
		// Mass
		if(particle_type == ParticlesType::H || particle_type == ParticlesType::Proton){
			mass =  datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
		}else if(particle_type == ParticlesType::Electron){
			mass = datum::m_e * datum::c_0 * datum::c_0 / datum::eV;
		}else{
			mass =  2 * datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
		}
		// Hard spheres cross section
		if(particle_type == ParticlesType::H || particle_type == ParticlesType::H_2 || particle_type == ParticlesType::H_2Plus){
			hard_speres_cross_section = datum::a_0 * datum::a_0 * 1e4 / 4;
		}
	}

	const ParticlesType particle_type;
	const CollisionModel collision_model;
	double mass;
	double hard_speres_cross_section;
};
