#include "ReconstructionStructure.hpp"


ReconstructionStructure::ReconstructionStructure(){
	number_voxels_grid = 0;
}

ReconstructionStructure::~ReconstructionStructure(){
}

void ReconstructionStructure::CreateGridMinimal(vector< vector<double> >& boundingbox, double division){

	inter_voxel_distance = division;
	BB = boundingbox;
	initial_offset.resize(3, 0);
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}

	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];

	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);
}



int_type_t ReconstructionStructure::ReturnIndexFromXYZIndices(int_type_t x, int_type_t y, int_type_t z){
	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z >= number_voxels_grid){
		cout << "ERROR ON size " <<x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z << ", n voxesl " << number_voxels_grid << endl;
		cout << "x , y, z " << x << ", " << y << ",  " << z << endl;
		exit(1);
	}

	return x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z;

}

void ReconstructionStructure::ReturnXYZIndicesFromIndex(int_type_t voxel_index, int_type_t& x, int_type_t& y, int_type_t& z){

	int_type_t temp_index = voxel_index;
	x = temp_index/(number_voxels_per_dim[1]*number_voxels_per_dim[2]);

	temp_index -= x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]);
	y = temp_index/(number_voxels_per_dim[2]);

	temp_index -= y*(number_voxels_per_dim[2]);

	z = temp_index;

	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z != voxel_index){
		cout << "ERROR on vox computation! " << endl << endl;
		exit(1);
	}
}


