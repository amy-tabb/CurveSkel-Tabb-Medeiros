#include "SubSteps.hpp"
#include "DistanceTransforms.hpp"

void SortSetsWithStops(vector<SkelGraph>& SG, double* bfs_labels,  vector<int_type_t>& F_set, vector<int_type_t>& F1_prime, vector<int_type_t>& F2_prime,
		double l_min_final, bool* frontier,  DISTANCE dist_type){

	double self_label;
	SkelGraph* self_ptr = 0;
	int_type_t self_index;
	int_type_t n_index;

	int_type_t j, jn;

	for (int_type_t i = 0, in = F_set.size(); i < in; i++){
		self_index = F_set[i];
		self_label = bfs_labels[self_index];
		self_ptr = &SG[self_index];

		if (frontier[self_index] == false){
			cout << "Error!  Frontier not correct " << endl;
			char ch; cin >> ch;
		}

		// move to F2
		if (self_label < l_min_final){
			for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

				n_index = self_ptr->neighbors[j]->voxel_id;

				if (frontier[n_index] == false && bfs_labels[n_index] > self_label){

					F2_prime.push_back(n_index);
					frontier[n_index] = true;
				}
			}
		} else {
			F1_prime.push_back(self_index);
		}
	}


	for (int_type_t i = 0, in = F_set.size(); i < in; i++){

		if (bfs_labels[F_set[i]] < l_min_final){
			frontier[F_set[i]] = false;
		}
	}


}



double  BFS_general_outside_to_inside_master_induction1(vector<SkelGraph>& SG, double* bfs_labels, vector<int_type_t>& grid_structure, int_type_t big_number,
		int_type_t connectivity, bool* frontier_map, DISTANCE dist_type, bool parallel){

	int_type_t n = SG.size();

	big_number = pow(grid_structure[0], 2) + pow(grid_structure[1], 2) + pow(grid_structure[2], 2);

	DistanceTransform1Parallel(SG, bfs_labels, grid_structure, big_number, connectivity);

	double max_distance = 0;
	for (int_type_t i = 0; i < n; i++){
		if (bfs_labels[i] > max_distance && bfs_labels[i] < big_number){
			max_distance = bfs_labels[i];
		}
	}

	return max_distance;

}


int_type_t FindALocalMaximaForSeed(vector<SkelGraph>& SG, double* bfs_labels, int_type_t big_number){

	int_type_t n = SG.size();
	int_type_t best_index =0;

	double max_distance = 0;
	for (int_type_t i = 0; i < n; i++){
		if (bfs_labels[i] > max_distance && bfs_labels[i] < big_number){
			max_distance = bfs_labels[i];
			best_index = i;
		}
	}

	return best_index;
}


void BFS_search_from_start_set_detail(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse, double* bfs_labels,
		vector<int_type_t>& F1, vector<int_type_t>& F2, bool* frontier, int_type_t big_number, DISTANCE dist_type){

	vector<int_type_t> F1_prime, F2_prime;
	int_type_t j = 0;
	int_type_t jn = 0;
	int_type_t n_index;
	double d = 0;
	double self_label;


	SkelGraph* self_ptr = 0;
	double l_min_final = 0;

	int_type_t self_index;

	l_min_final= big_number;

	for (int_type_t i = 0, in = F2.size(); i < in; i++){
		self_index = F2[i];

		if (frontier[self_index] == true)
		{

			self_label = bfs_labels[self_index];
			self_ptr = &SG[self_index];
			for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

				n_index = self_ptr->neighbors[j]->voxel_id;
				d = self_ptr->neighbor_distanced[j] + bfs_labels_oTi_inverse[n_index];
				// we stop when we meet a skeleton voxel ....
				if (bfs_labels[n_index] > self_label + d){
					bfs_labels[n_index] = self_label + d;
				}

				if (bfs_labels[n_index] > self_label && frontier[n_index] == false && l_min_final > bfs_labels[n_index]){
					l_min_final = bfs_labels[n_index];
				}
			}
		}
	}


	// Find set N, and l_min for F1
	for (int_type_t i = 0, in = F1.size(); i < in; i++){


		self_index = F1[i];

		if (frontier[self_index] == true)
		{
			self_label = bfs_labels[self_index];
			self_ptr = &SG[self_index];

			for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

				n_index = self_ptr->neighbors[j]->voxel_id;
				if (frontier[n_index] == false && bfs_labels[n_index] > self_label && l_min_final > bfs_labels[n_index]){
					l_min_final = bfs_labels[n_index];
				}
			}
		}
	}


	F1.insert(F1.end(), F2.begin(), F2.end());

	SortSetsWithStops(SG, bfs_labels,  F1, F1_prime, F2_prime,
			l_min_final, frontier, dist_type);




	F1.swap(F1_prime);
	F2.swap(F2_prime);

}


int_type_t  BFS_search_from_start_set(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse,
		double* bfs_labels, int_type_t big_number,
		bool* frontier_map, vector<int_type_t>& start_set, DISTANCE dist_type, bool initialize){

	// We can use this for the intial BFs as well as the update, simply set the initialize flag.

	int_type_t n = SG.size();
	for (int_type_t i = 0; i < n; i++){
		if (initialize){
			bfs_labels[i] = big_number;
		}
		frontier_map[i] = false;
	}

	vector<int_type_t> F1, F2;

	for (int_type_t i = 0, in = start_set.size(); i < in; i++){
		F2.push_back(start_set[i]);
		bfs_labels[start_set[i]] = 1;
		frontier_map[start_set[i]] = true;
	}


	int loop_counter = 0;
	while (F2.size() > 0){

		BFS_search_from_start_set_detail(SG, bfs_labels_oTi_inverse,
				bfs_labels, F1, F2, frontier_map,
				big_number, dist_type);
		loop_counter++;

	}


	return 0;

}



void FindLocalMaximaSurfaceSingle(vector<SkelGraph>& SG, double* bfs_labels, bool* already_explored,
		vector<int_type_t>& local_maxima,
		vector<int_type_t>& grid_structure, int_type_t big_number, int max_connectivity, bool* already_treated){

	int_type_t n = SG.size();
	vector<int_type_t> one_connected_component;


	//already_explored  = already_treated;
	OverwriteFirstArrayWithSecond<bool>(already_explored, already_treated, n);


	double self_label;
	int_type_t j, jn;
	int_type_t n_index;
	SkelGraph* self_ptr = 0;
	int_type_t x, y, z;
	int_type_t av_x, av_y, av_z;
	int_type_t temp_index;

	double best_distance;
	double current_distance;
	int_type_t best_index;
	int_type_t max_distance_index = n;

	double max_distance = 0;
	for (int_type_t i = 0; i < n; i++){
		// late -- revise this -- we just copied already treated to already exlored
		if (bfs_labels[i] > max_distance && bfs_labels[i] != big_number && already_explored[i] == false && already_treated[i] == false){
			//if (bfs_labels[i] > max_distance && bfs_labels[i] != big_number && already_explored[i] == false ){
			max_distance = bfs_labels[i];
			max_distance_index = i;

		}
	}
	//cout << "Max squared distance " << max_distance << endl; //char ch; cin >> ch;

	bool is_max;
	double d = 0;

	//for (int_type_t d = max_distance; d > lower_bound; d--){
	//cout << "d " << d << endl;
	int_type_t i = max_distance_index;

	self_label = bfs_labels[i];

	if (already_explored[i] == false && self_label != big_number){
		//cout << "First level 215 " << endl;
		self_ptr = &SG[i];
		d = self_label;

		if (int(self_ptr->neighbors.size()) < max_connectivity){
			//cout << "Second level 220" << endl;
			is_max = true;
			for (j = 0, jn = self_ptr->neighbors.size(); j < jn && is_max == true; j++){
				if (bfs_labels[self_ptr->neighbors[j]->voxel_id] > self_label){
					//cout << "This one has a greater distance... " << bfs_labels[self_ptr->neighbors[j]->voxel_id] << " already treated? " << already_treated[self_ptr->neighbors[j]->voxel_id] << endl;
					is_max = false;
				}
			}

			already_explored[i] = true;
			// need to explore all of them to make sure that this is a local max
			if (is_max){
				//cout << "Third level 231" << endl;
				one_connected_component.push_back(i);

				for (int_type_t p = 0; p < one_connected_component.size() && is_max == true; p++){

					self_ptr = &SG[one_connected_component[p]];


					for (j = 0, jn = self_ptr->neighbors.size(); j < jn && is_max == true; j++){
						n_index = self_ptr->neighbors[j]->voxel_id;

						if (bfs_labels[n_index] > self_label){
							is_max = false;
						}	else {
							// bfs_labels[n_index] == self_label will never happen b/c we're working with doubles
							if (already_explored[n_index] == false && bfs_labels[n_index] == self_label){
								already_explored[n_index] = true;
								one_connected_component.push_back(n_index);
							}
						}
					}

				}

				if (is_max){
					//cout << "number of components .... " << one_connected_component.size() <<  endl;
					//cout << "Value of max " << self_label << endl;

					av_x = 0;
					av_y = 0;
					av_z = 0;

					for (int_type_t p = 0; p < one_connected_component.size();  p++){
						temp_index = SG[one_connected_component[p]].grid_id;

						x = temp_index/(grid_structure[1]*grid_structure[2]);

						temp_index -= x*(grid_structure[1])*(grid_structure[2]);
						y = temp_index/(grid_structure[2]);

						temp_index -= y*(grid_structure[2]);

						z = temp_index;

						av_x += x;
						av_y += y;
						av_z += z;


					}


					av_x /= one_connected_component.size();
					av_y /= one_connected_component.size();
					av_z /= one_connected_component.size();

					// find closest one ...
					best_distance = max_distance*max_distance;
					best_index = one_connected_component.size();

					for (int_type_t p = 0; p < one_connected_component.size();  p++){
						temp_index = SG[one_connected_component[p]].grid_id;

						x = temp_index/(grid_structure[1]*grid_structure[2]);

						temp_index -= x*(grid_structure[1])*(grid_structure[2]);
						y = temp_index/(grid_structure[2]);

						temp_index -= y*(grid_structure[2]);

						z = temp_index;

						// perhaps there's a better way to do this ....
						current_distance = SquaredDifferenceOfUints(x, av_x) + SquaredDifferenceOfUints(y, av_y)  + SquaredDifferenceOfUints(z, av_z);

						if (current_distance < best_distance){
							best_distance = current_distance;
							best_index = p;
						}
					}

					if (best_index != one_connected_component.size()){
						local_maxima.push_back(one_connected_component[best_index]);
					}	else {
						cout << "Error! best index not found for connected component averaging" << endl;
						exit(1);
					}
				}
			}

			one_connected_component.clear();
		}
	}
}


void SortSetsWithStopsNeighborStop(vector<SkelGraph>& SG, double* bfs_labels,  vector<int_type_t>& F_set, vector<int_type_t>& F1_prime, vector<int_type_t>& F2_prime,
		double l_min_final, bool* frontier, bool* stop_lines, DISTANCE dist_type){

	double self_label;
	SkelGraph* self_ptr = 0;
	int_type_t self_index;
	int_type_t n_index;

	int_type_t j, jn;
	bool found_stop_ancestor;

	for (int_type_t i = 0, in = F_set.size(); i < in; i++){

		self_index = F_set[i];
		self_label = bfs_labels[self_index];
		self_ptr = &SG[self_index];

		if (stop_lines[self_index] == false){
			if (frontier[self_index] == false){
				cout << "Error!  Frontier not correct " << endl;
				char ch; cin >> ch;
			}

			// move to F2
			if (self_label < l_min_final){
				for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

					n_index = self_ptr->neighbors[j]->voxel_id;

					// only add a certain
					if (frontier[n_index] == false && bfs_labels[n_index] > self_label){

						found_stop_ancestor = false;

						if (found_stop_ancestor == false){
							F2_prime.push_back(n_index);
							frontier[n_index] = true;
						}
					}
				}
			} else {
				// move to F1

				found_stop_ancestor = false;
				if (found_stop_ancestor == false){
					F1_prime.push_back(self_index);
				}	else {
					cout << "Did not add b/c ancestor was part of stop set." << endl;
					char df; cin >> df;
				}
			}
		}
	}


	for (int_type_t i = 0, in = F_set.size(); i < in; i++){

		if (stop_lines[F_set[i]] == false){

			if (bfs_labels[F_set[i]] < l_min_final){
				frontier[F_set[i]] = false;
			}
		}

	}

}


void MarkConnectedComponentsFromFrontier_detail(vector<SkelGraph>& SG,  vector< int_type_t >& current_voxels,
		bool* frontier, bool* stop_lines){

	vector< int_type_t > next_layer;
	int_type_t j = 0;
	int_type_t jn = 0;
	int_type_t n_index;

	SkelGraph* self_ptr = 0;
	for (int_type_t i = 0, in = current_voxels.size(); i < in; i++){

		self_ptr = &SG[current_voxels[i]];



		for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

			n_index = self_ptr->neighbors[j]->voxel_id;

			if (frontier[n_index] == true && stop_lines[n_index] == false){
				next_layer.push_back(n_index);
				stop_lines[n_index] = true;
			}
		}
	}


	current_voxels.clear();
	current_voxels.swap(next_layer);

}


void  MarkConnectedComponentsFromFrontier(vector<SkelGraph>& SG, int_type_t start_index,
		bool* frontier, bool* stop_lines){

	// all of the nodes are marked with path_id = 1

	// idea do level sets .... also update this so that the frontier is more accurate.. Now we're pushing a traditional bfs
	// and the values get replace often.
	//int_type_t n = SG.size();

	vector<int_type_t> current_voxels;

	current_voxels.push_back(start_index);
	stop_lines[start_index]= true;

	int loop_counter = 0;
	// && loop_counter < 12
	while (current_voxels.size() > 0 ){

		MarkConnectedComponentsFromFrontier_detail(SG, current_voxels, frontier, stop_lines);
		loop_counter++;

	}

}




void BFS_interior_start_detail_for_skeleton(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse,
		double* bfs_labels, vector<int_type_t>& F1, vector<int_type_t>& F2, double& l_min1, bool* pre_skeleton,
		bool* stop_lines,
		bool* frontier, int_type_t big_number, DISTANCE dist_type){

	vector<int_type_t> F1_prime, F2_prime;
	int_type_t j = 0;
	int_type_t jn = 0;
	int_type_t n_index;
	double d = 0;
	double self_label;
	bool skeleton_found;


	SkelGraph* self_ptr = 0;
	double l_min_final = 0;

	int_type_t self_index;

	l_min_final= big_number;
	for (int_type_t i = 0, in = F2.size(); i < in; i++){

		self_index = F2[i];

		if (frontier[self_index] == true)
		{


			self_label = bfs_labels[self_index];
			self_ptr = &SG[self_index];

			//cout << "self label " <<  self_label << endl;

			// is any of those voxel's neighbors a skeleon voxel?
			skeleton_found = false;

			if (self_ptr->path_id == 0){

				for (j = 0, jn = self_ptr->neighbors.size(); j < jn && skeleton_found == false; j++){
					if (self_ptr->neighbors[j]->path_id > 0){
						skeleton_found = true;
					}
				}
			}

			if (skeleton_found){
				pre_skeleton[self_index] = true;
				stop_lines[self_index] = true;
				//cout << "Marking connected components " << endl;
				// fill in stop lines on the frontier.
				MarkConnectedComponentsFromFrontier(SG, self_index, frontier, stop_lines);
				//cout << "YES a skeleton node is found" << endl;
				//cout << "Return from mark. " << endl;
				//char ch; cin >> ch;
			}
			else
			{
				for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

					n_index = self_ptr->neighbors[j]->voxel_id;
					d = self_ptr->neighbor_distanced[j] + bfs_labels_oTi_inverse[n_index];
					// we stop when we meet a skeleton voxel ....
					if (bfs_labels[n_index] > self_label + d){
						bfs_labels[n_index] = self_label + d;

						//some_changed = true;
						//cout << "Update " << n_index << " by " << self_index << ", " <<  bfs_labels[n_index] << endl;
					}

					if (bfs_labels[n_index] > self_label && frontier[n_index] == false && l_min_final > bfs_labels[n_index]){
						l_min_final = bfs_labels[n_index];
					}

					//					if (frontier[n_index] == false && bfs_labels[n_index] > bfs_labels[self_index] && l_min2 > bfs_labels[n_index]){
					//						l_min2 = bfs_labels[n_index];
					//						//cout << "Updating lmin2, self and neighbor indices " << self_label << ", " << bfs_labels[n_index] << endl;
					//					}

				}
			}
		}
	}


	// Find set N, and l_min for F1
	for (int_type_t i = 0, in = F1.size(); i < in; i++){


		self_index = F1[i];

		if (frontier[self_index] == true)
		{
			self_label = bfs_labels[self_index];
			self_ptr = &SG[self_index];

			for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

				n_index = self_ptr->neighbors[j]->voxel_id;
				if (frontier[n_index] == false && bfs_labels[n_index] > self_label && l_min_final > bfs_labels[n_index]){
					l_min_final = bfs_labels[n_index];
				}
			}
		}
	}


	F1.insert(F1.end(), F2.begin(), F2.end());

	//	SortSetsWithStops(SG, bfs_labels,  F1, F1_prime, F2_prime,
	//			l_min_final, frontier,  stop_lines);

	//	SortSetsWithStops(SG, bfs_labels,  F1, F1_prime, F2_prime,
	//			l_min_final, frontier, stop_lines, dist_type);

	// ONLY CHANGE HERE
	SortSetsWithStopsNeighborStop(SG, bfs_labels,  F1, F1_prime, F2_prime,
			l_min_final, frontier, stop_lines, dist_type);




	F1.swap(F1_prime);
	F2.swap(F2_prime);

}
//
///////// KEEP ICCV last working version OCt 15 -- updated for CVPR
//int_type_t  BFS_search_from_tip_to_existing_skeleton(vector<SkelGraph>& SG, int_type_t* bfs_labels_oTi_inverse,
//		int_type_t* bfs_labels, int_type_t big_number, int_type_t start_index,
//		bool* pre_skeleton, bool* stop_lines, bool* frontier_map, DISTANCE dist_type){
//
//	// all of the nodes are marked with path_id = 1
//
//	// idea do level sets .... also update this so that the frontier is more accurate.. Now we're pushing a traditional bfs
//	// and the values get replace often.
//
//
//	int_type_t n = SG.size();
//	for (int_type_t i = 0; i < n; i++){
//		bfs_labels[i] = big_number;
//		pre_skeleton[i] = false;
//		stop_lines[i] = false;
//		frontier_map[i] = false;
//	}
//
//	vector<int_type_t> F1, F2;
//
//	F2.push_back(start_index);
//	bfs_labels[start_index] = 1;
//	frontier_map[start_index] = true;
//
//	int_type_t l_min1 = big_number;
//
//	int loop_counter = 0;
//	// && loop_counter < 12
//	while (F2.size() > 0){
//		//cout << "Before interior for " << loop_counter << endl;
//
//		BFS_interior_start_detail_for_skeleton(SG, bfs_labels_oTi_inverse,
//				bfs_labels, F1, F2, l_min1, pre_skeleton, stop_lines, frontier_map,
//				big_number, dist_type);
//
//		//		if (loop_counter % 1000 == 0){
//		//			cout << "Current number voxels for next iter " << current_voxels.size() << endl;
//		//		}
//
//		loop_counter++;
//
//	}
//
//	// just for visualization
//	//stop_lines = frontier_map;
//
//
//	return 0;
//
//}

int_type_t  BFS_search_from_tip_to_existing_skeleton(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse,
		double* bfs_labels, int_type_t big_number, int_type_t start_index,
		bool* pre_skeleton, bool* stop_lines, bool* frontier_map, DISTANCE dist_type){

	int_type_t n = SG.size();
	for (int_type_t i = 0; i < n; i++){
		bfs_labels[i] = big_number;
		pre_skeleton[i] = false;
		stop_lines[i] = false;
		frontier_map[i] = false;
	}

	vector<int_type_t> F1, F2;

	F2.push_back(start_index);
	bfs_labels[start_index] = 1;
	frontier_map[start_index] = true;

	double l_min1 = big_number;

	int loop_counter = 0;

	while (F2.size() > 0){

		BFS_interior_start_detail_for_skeleton(SG, bfs_labels_oTi_inverse,
				bfs_labels, F1, F2, l_min1, pre_skeleton, stop_lines, frontier_map,
				big_number, dist_type);
		loop_counter++;

	}


	return 0;

}

void ClearConnectedComponent_detail(vector<SkelGraph>& SG, vector<int_type_t>& F1,
		bool* stop_lines){

	vector<int_type_t> F1_prime;
	int_type_t j = 0;
	int_type_t jn = 0;
	int_type_t n_index;
	SkelGraph* self_ptr = 0;
	int_type_t self_index;



	for (int_type_t i = 0, in = F1.size(); i < in; i++){
		self_index = F1[i];


		self_ptr = &SG[self_index];

		for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

			n_index = self_ptr->neighbors[j]->voxel_id;

			if (stop_lines[n_index] == true){
				stop_lines[n_index] = false;
				F1_prime.push_back(n_index);
			}
		}

	}

	F1.swap(F1_prime);


}


int_type_t  ClearConnectedComponent(vector<SkelGraph>& SG, int_type_t start_index, bool* stop_lines){


	vector<int_type_t> F1;

	F1.push_back(start_index);


	int loop_counter = 0;
	// && loop_counter < 12
	while (F1.size() > 0){

		ClearConnectedComponent_detail(SG, F1, stop_lines);
		loop_counter++;

	}
	return 0;

}
void FillInTreatedMap(vector<SkelGraph>& SG, double* previous_bfs_map, bool* already_treated,
		int_type_t big_number){


	int_type_t n = SG.size();
#pragma omp parallel for
	for (int_type_t i = 0; i < n; i++){
		if (previous_bfs_map[i] < big_number){
			already_treated[i] = true;
		}
	}
}


void IdentifyLocalMinimaFromPreSkeletonsSetInKeepUpdatedComponents(vector<SkelGraph>& SG, bool* in_pre_skeleton, bool* updated_by_skeleton,
		vector<int_type_t>& pre_skeleton,
		double* bfs_map,
		vector< vector<int_type_t> >& keep_components){

	vector<int_type_t> confirmed_pre_skeletons;
	int_type_t self_index;
	//int_type_t self_distance;
	int_type_t j, jn;
	SkelGraph* self_ptr;

	//bool is_local_min;
	// eliminate until the local minimas are alone ....
	int_type_t number_nodes = SG.size();
	double smallest_distance = number_nodes*10;
	int_type_t current_pre_skel_index;

	vector<int_type_t> current_component;

	for (int_type_t i = 0; i < number_nodes; i++){
		if (updated_by_skeleton[i] == true){
			//cout << "starting connected component" << endl;
			smallest_distance = number_nodes*10;

			updated_by_skeleton[i] = false;
			current_component.push_back(i);
			keep_components.push_back(current_component);
			current_pre_skel_index = number_nodes;

			while (current_component.size() > 0){
				self_index = current_component.back();
				current_component.pop_back();

				// does this one in the pre_skeleton?
				self_ptr = &SG[self_index];

				if (in_pre_skeleton[self_index] == true){
					if (bfs_map[self_index] < smallest_distance){
						smallest_distance = bfs_map[self_index];
						current_pre_skel_index = self_index;
					}
				}

				for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

					//if (in_pre_skeleton[self_ptr->neighbors[j]->voxel_id] && updated_by_skeleton[self_ptr->neighbors[j]->voxel_id] && bfs_map[self_ptr->neighbors[j]->voxel_id] < self_distance){
					if (updated_by_skeleton[self_ptr->neighbors[j]->voxel_id] == true){
						current_component.push_back(self_ptr->neighbors[j]->voxel_id);
						keep_components.back().push_back(self_ptr->neighbors[j]->voxel_id);
						updated_by_skeleton[self_ptr->neighbors[j]->voxel_id] = false;
					}
				}
			}

			if (current_pre_skel_index != number_nodes){
				pre_skeleton.push_back(current_pre_skel_index);
			}	else {
				cout << "Error!  Pre skel not found .... " << endl;
				char ch; cin >> ch;
			}
		}
	}
}


void WalkPreSkeletonToTipToFindJunctionsWithDeepIntoCornersModificationConnectivity(vector<SkelGraph>& SG, double* bfs_labels, double* bfs_labels_dis,
		int_type_t tip_index,
		vector<int_type_t>& path_contents){

	SkelGraph* self_ptr = 0;
	SkelGraph* next_ptr =0;
	//SkelGraph* next_ptr =0;
	//int_type_t n_index;
	double best_distance;
	int_type_t j, jn;
	//bool other_skel_found;
	double self_distance;
	// just one for now

	self_ptr = &SG[tip_index];
	double max_d;

	best_distance = bfs_labels[self_ptr->voxel_id];
	self_distance = best_distance;

	path_contents.push_back(self_ptr->voxel_id);

	next_ptr = self_ptr;
	while (next_ptr != 0){

		// initialize loop

		next_ptr = 0;
		self_distance = best_distance;
		max_d = 0;

		// we short circuit to a close skeleton if one is a neighbor to eliminate wierd complexes
		// TODO -- represent this analytically with a path ........
		for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

			if (bfs_labels[self_ptr->neighbors[j]->voxel_id] < self_distance  &&  bfs_labels_dis[self_ptr->neighbors[j]->voxel_id] > max_d){
				max_d = bfs_labels_dis[self_ptr->neighbors[j]->voxel_id];
			}
		}

		for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

			// leave this for now, but the == max_d is terrible for numerical precision reasons.
			if (bfs_labels[self_ptr->neighbors[j]->voxel_id] < best_distance && bfs_labels_dis[self_ptr->neighbors[j]->voxel_id] >= max_d - 0.1){
				best_distance = bfs_labels[self_ptr->neighbors[j]->voxel_id];
				next_ptr = self_ptr->neighbors[j];
			}
		}


		if (next_ptr == 0){
			//cout << "Error!  Next ptr is zero! " << endl; exit(1);
		}	else {
			self_ptr = next_ptr;
		}


		if (next_ptr != 0){
			next_ptr->path_id = 1;
			path_contents.push_back(next_ptr->voxel_id);
		}

	}

	vector<int_type_t> new_path;

	//cout << "Round 2 " << endl;
	//-----------------------------
	// now, go through and pick the fastest descent path without the choices
	//
	self_ptr = &SG[tip_index];
	next_ptr = self_ptr;
	best_distance = bfs_labels[tip_index];
	new_path.push_back(self_ptr->voxel_id);
	while (next_ptr != 0){

		// initialize loop

		next_ptr = 0;
		max_d = 0;

		for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

			// leave this for now, but the == max_d is terrible for numerical precision reasons.
			if (bfs_labels[self_ptr->neighbors[j]->voxel_id] < best_distance && self_ptr->neighbors[j]->path_id == 1){
				best_distance = bfs_labels[self_ptr->neighbors[j]->voxel_id];
				next_ptr = self_ptr->neighbors[j];
			}
		}


		if (next_ptr == 0){
			//cout << "Error!  Next ptr is zero! " << endl; exit(1);
		}	else {
			//next_ptr->path_id = current_path_id;

			self_ptr = next_ptr;
		}


		if (next_ptr != 0){
			//next_ptr->path_id = 1;
			new_path.push_back(next_ptr->voxel_id);
		}

	}

	for (int_type_t i = 0, in = path_contents.size(); i < in; i++){
		SG[path_contents[i]].path_id = 0;
	}
	path_contents.swap(new_path);
}

void FindLinks(vector<SkelGraph>& SG, double* bfs_labels, int_type_t tip_index,
		vector<int_type_t>& path_contents){

	SkelGraph* self_ptr = 0;
	SkelGraph* next_ptr =0;
	//SkelGraph* next_ptr =0;
	//int_type_t n_index;
	double best_distance;
	int_type_t j, jn;
	bool other_skel_found;
	double self_distance;
	// just one for now

	self_ptr = &SG[tip_index];

	//cout << "Starting from distance " << distance_node_map[i].first << endl;

	// make this one have the same path id so ccs are easy ...
	//self_ptr->path_id = current_path_id;
	other_skel_found = false;
	// distance will always decrease from self
	best_distance = bfs_labels[self_ptr->voxel_id];
	self_distance = best_distance;

	path_contents.push_back(self_ptr->voxel_id);

	while (other_skel_found == false){

		// initialize loop

		next_ptr = 0;
		self_distance = best_distance;

		// we short circuit to a close skeleton if one is a neighbor to eliminate wierd complexes
		// TODO -- represent this analytically with a path ........
		for (j = 0, jn = self_ptr->neighbors.size(); j < jn && other_skel_found == false; j++){
			//for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){
			//cout << "Next distance " << bfs_labels[self_ptr->neighbors[j]->voxel_id] << ", self distance " << self_distance << " path id " << self_ptr->neighbors[j]->path_id << endl;


			if (bfs_labels[self_ptr->neighbors[j]->voxel_id] < best_distance){
				best_distance = bfs_labels[self_ptr->neighbors[j]->voxel_id];
				next_ptr = self_ptr->neighbors[j];
			}

			// path_id > 0 means part of an existing path.
			if (self_ptr->neighbors[j]->path_id > 0 && bfs_labels[self_ptr->neighbors[j]->voxel_id] < self_distance){
				//cout << bfs_labels[self_ptr->neighbors[j]->voxel_id] << " neighbor that is skeleton distance " << endl;
				other_skel_found = true;
			}
		}

		//cout << "Next distance " << best_distance << endl;




		if (other_skel_found == false){
			if (next_ptr == 0){
				cout << "in Find links " << endl;
				cout << "Error!  Next ptr is zero! " << endl; exit(1);
			}
			//next_ptr->path_id = current_path_id;
			self_ptr = next_ptr;
		}	else {
			// there may be several options of the other skel found .. we choose the best one .

			next_ptr = 0;
			self_distance = bfs_labels[self_ptr->voxel_id];
			best_distance = self_distance;
			//cout << "Best distance, self distance " << best_distance << ", " << self_distance << endl;
			for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){
				//for (j = 0, jn = self_ptr->neighbors.size(); j < jn; j++){

				if (self_ptr->neighbors[j]->path_id > 0 && bfs_labels[self_ptr->neighbors[j]->voxel_id] < self_distance){
					if (bfs_labels[self_ptr->neighbors[j]->voxel_id] < best_distance){
						best_distance = bfs_labels[self_ptr->neighbors[j]->voxel_id];
						next_ptr = self_ptr->neighbors[j];
						//cout << "Assigned next " << best_distance << endl;
					}
				}
			}

			if (next_ptr == 0){
				cout << "Error!  Next ptr is zero 906! " << endl; exit(1);
			}

		}


		path_contents.push_back(next_ptr->voxel_id);

	}
	//cout << "Skel found? " << other_skel_found << endl;
}



double FindMaxInBFSMap(vector<int_type_t>& set_to_consider,	double* bfs_map){

	double max = 0;

	for (int_type_t i = 0, in = set_to_consider.size(); i < in; i++){

		if (bfs_map[set_to_consider[i]] > max){
			max = bfs_map[set_to_consider[i]];
		}
	}

	return max;
}


void ReturnXYZFromIndicesSkelGraph(SkelGraph& S, int_type_t& x, int_type_t& y, int_type_t& z, vector<int_type_t>& grid_structure){
	int_type_t temp_index = S.grid_id;

	x = temp_index/(grid_structure[1]*grid_structure[2]);

	temp_index -= x*(grid_structure[1])*(grid_structure[2]);
	y = temp_index/(grid_structure[2]);

	temp_index -= y*(grid_structure[2]);

	z = temp_index;
}



bool ClassifyUsingChiSquaredAssumption(vector<SkelGraph>& SG, int_type_t skeleton_index, bool* stop_lines,
		vector<vector<int_type_t> >& paths, vector<int_type_t>& grid_structure,
		int_type_t big_big_number, int_type_t connectivity, vector<int_type_t>& proposed_path, double probability_threshold){

	vector<double> v_diff;
	vector<double> v(3,0);
	vector<double> v_sum(3,0);
	vector<double> summed_for_all(3,0);
	vector<double> composite_center(3, 0);
	vector<double> one(3, 0);
	vector<double> point(4, 1);
	int_type_t x, y, z;
	int_type_t x0, y0, z0;

	int_type_t number_nodes = SG.size();

	int_type_t skel_index = proposed_path[0];

	int_type_t number_in_stop_set = 0;

	ReturnXYZFromIndicesSkelGraph(SG[skel_index], x0, y0, z0, grid_structure);

	for (int_type_t k = 0; k < number_nodes; k++){

		if (stop_lines[k] == true && SG[k].neighbors.size() < connectivity){
			v[0] = 0;  v[1] = 0; v[2] = 0;
			number_in_stop_set++;
			// find the closest skeletal node ...
			ReturnXYZFromIndicesSkelGraph(SG[k], x, y, z, grid_structure);
			v[0] = double(x) - double(x0);
			v[1] = double(y) - double(y0);
			v[2] = double(z) - double(z0);

			v_sum[0] += v[0];
			v_sum[1] += v[1];
			v_sum[2] += v[2];
			v_diff.push_back(v[0]);
			v_diff.push_back(v[1]);
			v_diff.push_back(v[2]);

		}
	}

	vector<double> mu0(3,0);
	mu0[0] = v_sum[0]/double(number_in_stop_set);
	mu0[1] = v_sum[1]/double(number_in_stop_set);
	mu0[2] = v_sum[2]/double(number_in_stop_set);

	vector<vector<double> > Sigma(3, vector<double>(3, 0));
	vector<vector<double> > Sv(3, vector<double>(3, 0));
	vector<vector<double> > temp_product(3, vector<double>(3, 0));
	vector< vector< double> > diffT(3, vector<double>(1, 0));  vector< vector< double> > diff(1, vector<double>(3, 0));

	for (int_type_t i = 0; i < number_in_stop_set; i++){

		for (int j = 0; j < 3; j++){
			diff[0][j] = mu0[j] - v_diff[3*i + j];
			diffT[j][0] = mu0[j] - v_diff[3*i + j];
		}

		MultiplyMatricesWithSizes(diffT, diff, 3, 1, 3, temp_product);

		AddMatrixToMatrix(Sigma, temp_product, Sigma, 3, 3);
	}

	bool well_behaved = true;

	for (int i = 0; i < 3; i++){
		if (Sigma[i][i] < (1e-8)){
			well_behaved = false;
		}
	}

	if (well_behaved == false){
		for (int i = 0; i < 3; i++){
			Sigma[i][i] += 0.1;
		}

	}

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			Sigma[i][j] /= double(number_in_stop_set);
		}
	}


	cv::Mat sigma(3, 3, cv::DataType<double>::type);
	cv::Mat S(3, 3, cv::DataType<double>::type);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			sigma.at<double>(i, j) = Sigma[i][j];
		}
	}

	S = sigma.inv(cv::DECOMP_SVD);

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			Sv[i][j] = S.at<double>(i, j);
		}
	}


	vector<vector<double> > x1T(3,vector<double>(1, 0));
	vector< vector< double> > x1(1, vector<double>(3, 0));
	vector<vector<double> > d_var(1, vector<double>(1, 0));
	vector<vector<double> > intermediate(1, vector<double>(3, 0));
	double pdf = 0;
	bool passed = false;
	double denominator = pow(2, 1.5)*0.88622692545;

	double prod = 1;
	for (int_type_t pcounter = proposed_path.size(); pcounter > proposed_path.size() - 1; pcounter--){
		ReturnXYZFromIndicesSkelGraph(SG[proposed_path[pcounter-1]], x, y, z, grid_structure);


		v[0] = double(x) - double(x0);
		v[1] = double(y) - double(y0);
		v[2] = double(z) - double(z0);

		v[0] = mu0[0] - v[0];
		v[1] = mu0[1] - v[1];
		v[2] = mu0[2] - v[2];


		x1[0][0] = v[0];
		x1[0][1] = v[1];
		x1[0][2] = v[2];

		x1T[0][0] = x1[0][0];
		x1T[1][0] = x1[0][1];
		x1T[2][0] = x1[0][2];


		MultiplyMatricesWithSizes( x1, Sv, 1, 3, 3, intermediate);
		MultiplyMatricesWithSizes( intermediate, x1T, 1, 3, 1, d_var);


		{

			pdf = pow(d_var[0][0], 0.5)*exp(-d_var[0][0]/2.0)/denominator;
			prod = prod * pdf;

			if (pdf < probability_threshold && pcounter == proposed_path.size() ){
				passed = true;
			}
		}
	}

	return passed;
}





