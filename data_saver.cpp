#include "data_saver.h"

void RecordPhysParams_LocalCollTest(double temperature, double dencity, double Vmax,
		double init_distr_max_val){
	vec params(5, fill::zeros);
	params(0) = temperature;
	params(1) = dencity;
	params(2) = Vmax;
	params(3) = init_distr_max_val;
	params.save("params.bin", raw_binary);
}

string GetTimeID(const size_t t){
	string result = "";
	size_t a = t / 10;
	size_t b = t % 10;
	for(size_t i = 0; i < a; ++i){
		result += 'a';
	}
	result += to_string(b);
	return result;
}

void SaveCubeSlice(const cube& c, const set<size_t>& free_pos, const size_t fixed_pos, const string& file_name){
	size_t v_size = c.n_slices;
	mat for_saving(v_size, v_size, fill::zeros);
	for(size_t i = 0; i < v_size; ++i){
		for(size_t j = 0; j < v_size; ++j){
			if(free_pos.count(0) == 1 and free_pos.count(1) == 1){
				for_saving(j,i) = c(j,i,fixed_pos);
			}else if(free_pos.count(1) == 1 and free_pos.count(2) == 1){
				for_saving(j,i) = c(fixed_pos,j,i);
			}else if(free_pos.count(0) == 1 and free_pos.count(2) == 1){
				for_saving(j,i) = c(j,fixed_pos,i);
			}
		}
	}
	for_saving.save(file_name, raw_binary);
}
