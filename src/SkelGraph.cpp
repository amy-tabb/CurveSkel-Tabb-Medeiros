#include "SkelGraph.hpp"


double sqrt2 = sqrt(2.0);
double sqrt3 = sqrt(3.0);

SkelGraph::SkelGraph(){
	voxel_id = 0;
	path_id = 0;
	grid_id = 0;
	cc_id = 0;
}

SkelGraph::SkelGraph(int_type_t vi, int_type_t gi){
	voxel_id = vi;
	grid_id = gi;
	path_id = 0;
	cc_id = 0;
}

SkelGraph::~SkelGraph(){}

void SkelGraph::Print(){
	cout << "Voxel id " << voxel_id << endl;
	cout << "Grid id " << grid_id << endl;
	cout << "Path id " << path_id << endl;
	cout << "Number neighbors " << neighbors.size();
}

bool SkelGraph::IsLocalMaxima(vector<int_type_t>& labels){

	int_type_t self_label = labels[voxel_id];

	bool ismax = true;
	for (int_type_t j = 0, jn = neighbors.size(); j < jn && ismax; j++){
		if (labels[neighbors[j]->voxel_id] > self_label){
			ismax = false;
		}
	}
	return ismax;
}



int_type_t SquaredDifferenceOfUints(int_type_t A, int_type_t B){
	return A*A + B*B - 2*A*B;
}

pair<int_type_t, int_type_t> CreateGraphFromGridReturnPair(ReconstructionStructure& RS, vector<SkelGraph>& SG,
		vector<int_type_t>& cc_counts, DISTANCE dist_type){


	int_type_t* dense_to_sparse = new int_type_t[RS.number_voxels_grid]; // will be dealloc'd after this function.
	int_type_t sparse_counts = 0;

	for (int_type_t self_index = 0; self_index < RS.number_voxels_grid; self_index++){


		if (RS.configuration_grid[self_index] == 0){
			dense_to_sparse[self_index] = sparse_counts;
			SG.push_back(SkelGraph(sparse_counts, self_index));
			sparse_counts++;
		}	else {
			dense_to_sparse[self_index] = RS.number_voxels_grid;
		}
	}


	int_type_t arr[26];
	int_type_t distance[26];
	int_type_t number_ns;
	int_type_t n_index;
	int_type_t d;
	double dd;
	int_type_t j;
	int_type_t biggest_cc = 0;
	int_type_t index_biggest_cc = 0;
	vector<int_type_t> triple(3, 0);


	for (int_type_t voxel_id = 0; voxel_id < sparse_counts; voxel_id++){

		number_ns = Return26ConnectedNeighbors(RS, SG[voxel_id].grid_id, &arr[0], &distance[0]);

		for (j = 0; j < number_ns; j++){
			n_index = arr[j];

			d = distance[j];

			if (dense_to_sparse[n_index] != RS.number_voxels_grid){
				// link in
				SG[voxel_id].neighbors.push_back(&SG.at(dense_to_sparse[n_index]));

				// assumes regular voxels ...
				switch (dist_type){
				case L1: {
					SG[voxel_id].neighbor_distance.push_back(d);
					//cout << "d" << d << endl;
				} break;
				case L2_induction:{
					if (d == 1){
						dd = 1;
					}	else {
						if (d == 2){
							dd = sqrt2;
						}	else {
							dd = sqrt3;
						}
					}
					SG[voxel_id].neighbor_distanced.push_back(dd);
				}	break;
				case L2: {
					if (d == 1){
						dd = 1;
					}	else {
						if (d == 2){
							dd = sqrt2;
						}	else {
							dd = sqrt3;
						}
					}
					SG[voxel_id].neighbor_distanced.push_back(dd);
				} break;
				case L2_approx: {

					if (d == 1){
						dd = 1;
					}	else {
						if (d == 2){
							dd = sqrt2;
						}	else {
							dd = sqrt3;
						}
					}
					SG[voxel_id].neighbor_distanced.push_back(dd);
				} break;

				case L_g: {
					SG[voxel_id].neighbor_distance.push_back(1);
				} break;
				default: {
					cout << "ERROR!  faulty dist type" << endl;
				}
				}

			}	else {

			}
		}
	}

	int_type_t cc = 1;
	vector<int_type_t> to_search;
	int_type_t search_index;



	for (int_type_t voxel_id = 0; voxel_id < sparse_counts; voxel_id++){
		if (SG[voxel_id].cc_id == 0){

			to_search.push_back(voxel_id);
			search_index = 0;
			SG[voxel_id].cc_id = cc;

			while (search_index < to_search.size()){

				for (j = 0; j < SG[to_search[search_index]].neighbors.size(); j++){
					if (SG[to_search[search_index]].neighbors[j]->cc_id == 0){
						SG[to_search[search_index]].neighbors[j]->cc_id = cc;
						to_search.push_back(SG[to_search[search_index]].neighbors[j]->voxel_id);
					}
				}

				search_index++;
			}

			if (search_index > biggest_cc){
				biggest_cc = search_index;
				index_biggest_cc = cc;
			}
			cc_counts.push_back(search_index);
			cc++;

			to_search.clear();

		}

	}

	delete [] dense_to_sparse;

	return pair<int_type_t, int_type_t>(index_biggest_cc, cc);

}

void ReworkUsingOnlyBiggestConnectedComponent(ReconstructionStructure& RS, vector<SkelGraph>& SG, vector<SkelGraph>& newSG, int_type_t biggest_connected_component){

	int_type_t* dense_to_sparse = new int_type_t[SG.size()]; // will be dealloc'd after this function.
	int_type_t sparse_counts = 0;
	//char ch;

	for (int_type_t i = 0, in = SG.size(); i < in; i++){
		if (SG[i].cc_id == biggest_connected_component){
			dense_to_sparse[i] = sparse_counts;
			newSG.push_back(SkelGraph(sparse_counts, SG[i].grid_id));
			newSG[sparse_counts].neighbor_distance = SG[i].neighbor_distance;
			newSG[sparse_counts].neighbor_distanced = SG[i].neighbor_distanced;
			newSG[sparse_counts].cc_id = SG[i].cc_id;

			sparse_counts++;
		}	else {
			dense_to_sparse[i] = RS.number_voxels_grid;
			RS.configuration_grid[SG[i].grid_id] = true;
		}
	}

	int_type_t sparse_index;
	// now, copy over ....
	for (int_type_t i = 0, in = SG.size(); i < in; i++){
		if (SG[i].cc_id == biggest_connected_component){
			sparse_index = dense_to_sparse[i];
			for (int_type_t j = 0, jn = SG[i].neighbors.size(); j < jn; j++){
				newSG[sparse_index].neighbors.push_back(&newSG[dense_to_sparse[SG[i].neighbors[j]->voxel_id]]);
			}
		}
	}

	//SG[0].Print();
	delete [] dense_to_sparse;

}

void CreateSurfaceGraphFromGrid(ReconstructionStructure& RS, vector<SkelGraph>& SG){

	int_type_t arr[26];
	int_type_t distance[26];
	int_type_t number_ns;
	int_type_t n_index;
	int_type_t d;
	int_type_t j;
	int_type_t number_ns_2;


	int_type_t* dense_to_sparse = new int_type_t[RS.number_voxels_grid]; // will be dealloc'd after this function.
	int_type_t sparse_counts = 0;
	//char ch;
	for (int_type_t self_index = 0; self_index < RS.number_voxels_grid; self_index++){


		if (RS.configuration_grid[self_index] == 0){
			number_ns = Return26ConnectedNeighbors(RS, self_index, &arr[0], &distance[0]);

			number_ns_2 = 0;

			for (j = 0; j < number_ns; j++){
				n_index = arr[j];
				if (RS.configuration_grid[n_index] == false){
					number_ns_2++;
				}
			}

			if (number_ns_2 < 26)
			{
				dense_to_sparse[self_index] = sparse_counts;
				SG.push_back(SkelGraph(sparse_counts, self_index));
				sparse_counts++;
			}
			else {
				dense_to_sparse[self_index] = RS.number_voxels_grid;
			}
		}	else {
			dense_to_sparse[self_index] = RS.number_voxels_grid;
		}
	}


	for (int_type_t voxel_id = 0; voxel_id < sparse_counts; voxel_id++){

		number_ns = Return26ConnectedNeighbors(RS, SG[voxel_id].grid_id, &arr[0], &distance[0]);

		for (j = 0; j < number_ns; j++){
			n_index = arr[j];
			d = distance[j];

			if (dense_to_sparse[n_index] != RS.number_voxels_grid){
				// link in
				SG[voxel_id].neighbors.push_back(&SG.at(dense_to_sparse[n_index]));
				SG[voxel_id].neighbor_distance.push_back(d);
			}	else {
				// not a surface node .... ignore
			}
		}


	}

	cout << "Number of connected nodes .... " << SG.size() << endl;
	delete [] dense_to_sparse;

}

int Return26ConnectedNeighbors(ReconstructionStructure& RS, int_type_t start_voxel, int_type_t* member_array, int_type_t* distance){

	int number_neighbors = 0;
	int_type_t nindex;

	int_type_t x_index, y_index, z_index;

	RS.ReturnXYZIndicesFromIndex(start_voxel, x_index, y_index, z_index);


	int_type_t minx, maxx, miny, maxy, minz, maxz;

	if (x_index == 0){
		minx = 0;
	}	else {
		minx = x_index - 1;
	}

	if (y_index == 0){
		miny = 0;
	}	else {
		miny = y_index - 1;
	}

	if (z_index == 0){
		minz = 0;
	}	else {
		minz = z_index - 1;
	}

	if (z_index == RS.number_voxels_per_dim[2] - 1){
		maxz = z_index;
	}	else {
		maxz = z_index + 1;
	}

	if (y_index == RS.number_voxels_per_dim[1] - 1){
		maxy = y_index;
	}	else {
		maxy = y_index + 1;
	}

	if (x_index == RS.number_voxels_per_dim[0] - 1){
		maxx = x_index;
	}	else {
		maxx = x_index + 1;
	}

	for (int_type_t x0 = minx; x0 <= maxx; x0++){
		for (int_type_t y0 = miny; y0 <= maxy; y0++){
			for (int_type_t z0 = minz; z0 <= maxz; z0++){
				nindex = RS.ReturnIndexFromXYZIndices(x0, y0, z0);
				if (nindex != start_voxel){
					member_array[number_neighbors] = nindex;
					distance[number_neighbors] = (x0 != x_index) + (y0 != y_index) + (z0 != z_index);
					number_neighbors++;
				}
			}
		}
	}


	return number_neighbors;


}


bool distance_comp_function_skel_distance_double(pair<int_type_t, double> a, pair<int_type_t, double> b){
	return a.second < b.second;
}

