#pragma once

#include <armadillo>
#include <set>

using namespace std;
using namespace arma;

void RecordPhysParams_LocalCollTest(double temperature, double dencity, double Vmax,
		double init_distr_max_val);

string GetTimeID(const size_t t);

void SaveCubeSlice(const cube& c, const set<size_t>& free_pos, const size_t fixed_pos, const string& file_name);
