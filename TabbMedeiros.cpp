#include "TabbMedeiros.hpp"

#include <sys/stat.h>
#include <sys/time.h>
#include <chrono>
#include <opencv2/core/core.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <chrono>

#include "SubSteps.hpp"
#include "Includes.hpp"
#include "ReconstructionStructure.hpp"
#include "SkelGraph.hpp"
#include "ConversionFunctions.hpp"

using std::list;
using std::cout;
using std::endl;
using std::string;


int main(int argc, char **argv)
{
	/// program has 1 required and 2 optional argument
	// mandatory -- directory with the BB.txt and 0.txt files
	// optional 1 --  threshold for spurious curve segment classification.  if not specified, it is set to 1e-12.
	// optional 2 -- flag to indicate that the model needs to be converted from a directory of images to the format in 0.txt and BB.txt.  If conversion is needed, then threshold needs to be specified.

	double threshold = 1e-12;
	bool need_to_convert = false;
	string directory = "";

	if (argc >= 2){
		directory = argv[1];
	}	else {
		cout << "Not enough arguments.  prog_name directory threshold." << endl;
		exit(1);
	}

	if (argc >= 3){
		threshold = FromString<double>(argv[2]);
	}

	if (directory[directory.size() - 1] != '/'){
		directory = directory + "/";
	}

	if (argc == 4){
		need_to_convert = true;
	}

	if (need_to_convert){
		ReadIndividualImages(directory);
	}

	TabbMedeiros( L2_induction, directory, threshold);

	if (need_to_convert){
		WriteIndividualImages(directory);
	}

	cout << "Before final return" << endl;


	return(0);
}

void TabbMedeiros(DISTANCE dist_type, string supplied_write_directory, double threshold){

	string supplied_BBfile = supplied_write_directory + "BB.txt";
	string supplied_object_file = supplied_write_directory + "0.txt";

	ifstream in_test;
	in_test.open(supplied_BBfile.c_str());

	if (!in_test){
		cout << "You forgot the BB file or the path is wrong" << endl << supplied_BBfile << endl;;
		exit(1);
	}
	in_test.close();

	in_test.open(supplied_write_directory.c_str());
	if (!in_test){
		cout << "You forgot the object file or the path is wrong" << endl << supplied_object_file << endl;;
		exit(1);
	}
	in_test.close();


	string source_file;
	float division = 100.0;

	string write_directory ="";

	vector<double> pA(3, 0);
	vector<double> pB(3, 1000);

	ifstream in;
	vector<vector<double> > BB;
	string filename;

	string object_file;

	bool* frontier = 0;
	double* bfs_label_outside_to_inside_double = 0;
	double* bfs_label_outside_to_inside_inverse_double = 0;
	int_type_t big_big_number = 0;
	vector<int_type_t> local_maxima_sparse_rep_int;
	vector<double> local_maxima_sparse_rep_double;
	int_type_t greatest_OtI_distance_double = 0;
	int_type_t seed_node_index = 0;
	int_type_t best_value_int = 0;
	double best_value_double = 0;
	vector<int_type_t> temp_vector;
	double* bfs_label_node_to_node_double = 0;
	double greatest_mod_distance_double = 0;
	bool* already_explored = 0;
	bool* this_voxel_already_treated = 0;
	int loop_counter = 0;
	int_type_t local_tip;
	vector<int_type_t> pre_skeleton_nodes;
	vector<vector<int_type_t> > multiple_paths_from_same_start;
	double local_max_distance_double = 0;
	double* bfs_from_tip_double =0 ;
	bool* pre_skeleton_bool = 0;
	bool* stop_lines = 0;
	vector<pair<int_type_t, double> > pre_skel_distance_double;
	bool passed_spurious_tests;
	bool* stop_lines_copy = 0;
	vector<vector<int_type_t> > pending_paths;
	bool walks_through_another_stop_set0;
	bool walks_through_another_stop_set1;


	vector<int_type_t> linked_path;
	vector< vector<int_type_t> > keep_components;

	vector<int_type_t> connected_components_info;
	int_type_t biggest_cc;
	int_type_t path_counter;
	int_type_t last_path_added;
	pair<int_type_t, int_type_t> biggest_cc_total_cc_pair;
	int_type_t current_path_id;

	vector<vector<int_type_t> > paths;
	vector<int_type_t> current_path;
	vector<vector<double> > path_diameters_single;
	vector<vector<double> > path_diameters;
	vector<pair<int_type_t, int_type_t> > distance_node_map;
	vector<int_type_t> further_starts;
	int_type_t index_first_path;
	vector<SkelGraph> SkeletonGraph;
	int_type_t number_nodes;
	vector<int_type_t> local_maxes_skeleton_hypotheses;

	in.open(supplied_BBfile.c_str());

	in >> division;
	in >> pA[0] >> pA[1] >> pA[2];
	in >> pB[0] >> pB[1] >> pB[2];

	object_file = supplied_object_file;
	write_directory = supplied_write_directory;

	BB.push_back(pA);
	BB.push_back(pB);

	// this version only does L2 induction
	switch (dist_type){
	case L2_induction: {} break;
	default: {
		cout << "Bad cost function selection. " << endl;
		exit(1);
	}
	}

	string current_write_file = write_directory + "trace.txt";

	ReconstructionStructure R_model;
	int max_connectivity = 26;

	/////////////////////////////// Read Model //////////////////////////////////////////////////
	ifstream oin;
	oin.open(object_file.c_str());

	std::ofstream out;
	out.open(current_write_file.c_str());

	if (oin.fail()){
		oin.close();

		out << "file " << object_file << " returns error " << endl;
		exit(1);

	}	else {
		R_model.CreateGridMinimal(BB, division);

		string val;
		oin >> val;

		int_type_t index_i;
		while (val != "-1"){

			index_i = FromString<int_type_t>(val);
			R_model.configuration_grid[index_i] = false;

			oin >> val;
		}

		out << "Read from " << object_file << endl;
	}

	//////////////////////////////////////////////////////////  AFTER READING //////////////////////////
	// timer starts after loaded. ...
	auto t0 = std::chrono::high_resolution_clock::now();

	vector<bool> skeleton(R_model.number_voxels_grid, true);


	// 1.  create nodes from the grid -- with connections.  Assume that there's only one connected component.
	biggest_cc_total_cc_pair = CreateGraphFromGridReturnPair(R_model, SkeletonGraph, connected_components_info, dist_type);
	// We select the largest connected component
	biggest_cc = biggest_cc_total_cc_pair.first;

	if (biggest_cc_total_cc_pair.second > 2){
		cout << "Alert!  This file had more than one component; we're only dealing with the largest one in this implementation." << endl;
		vector<SkelGraph> SGswap;

		ReworkUsingOnlyBiggestConnectedComponent(R_model, SkeletonGraph, SGswap, biggest_cc);

		SGswap.swap(SkeletonGraph);
	}
	number_nodes = SkeletonGraph.size();

	// allocations of storage locations ... ///////////////////////////////

	bfs_label_outside_to_inside_double = new double[number_nodes];
	bfs_label_outside_to_inside_inverse_double = new double[number_nodes];
	bfs_label_node_to_node_double = new double[number_nodes];
	bfs_from_tip_double = new double[number_nodes];


	frontier = new bool[number_nodes];
	already_explored = new bool[number_nodes];
	this_voxel_already_treated = new bool[number_nodes];
	pre_skeleton_bool = new bool[number_nodes];
	stop_lines = new bool[number_nodes];
	stop_lines_copy = new bool[number_nodes];


	// Initialization
	for (int_type_t i = 0; i < number_nodes; i++){
		already_explored[i] =false;
		this_voxel_already_treated[i] = false;
	}
	for (int_type_t i = 0; i < number_nodes; i++){
		bfs_label_outside_to_inside_double[i] = 0;
		bfs_from_tip_double[i] = 0;
		frontier[i] = false;
	}

	big_big_number = R_model.number_voxels_grid*2;
	//////////////////////////////// End reading/initialization /////////////////////////////////

	// Step 1
	auto t_timer_distance0 = std::chrono::high_resolution_clock::now();
	// Step 1.1, compute d_i
	greatest_OtI_distance_double = BFS_general_outside_to_inside_master_induction1(SkeletonGraph, bfs_label_outside_to_inside_double, R_model.number_voxels_per_dim,
			big_big_number,
			max_connectivity, frontier, dist_type, true);
	auto t_timer_distance1 = std::chrono::high_resolution_clock::now();
	cout << "Time for distance computation without a grid ... " << std::chrono::duration_cast<std::chrono::milliseconds>(t_timer_distance1 - t_timer_distance0).count()	<< " milliseconds "<< endl;
	// Step 1.2, find v^*
	seed_node_index = FindALocalMaximaForSeed(SkeletonGraph, bfs_label_outside_to_inside_double, big_big_number);
	best_value_double = bfs_label_outside_to_inside_double[seed_node_index];
	local_maxima_sparse_rep_double.clear();
	local_maxima_sparse_rep_double.push_back(seed_node_index);

	// create inverse of the bfs_map to get w_i's
	for (int_type_t i = 0; i < number_nodes; i++){
		bfs_label_outside_to_inside_inverse_double[i] = greatest_OtI_distance_double - bfs_label_outside_to_inside_double[i];
	}


	// Step 2.1 -- compute BFS1 map
	temp_vector.resize(1, seed_node_index);
	greatest_mod_distance_double = BFS_search_from_start_set(SkeletonGraph,
			bfs_label_outside_to_inside_inverse_double, bfs_label_node_to_node_double,
			big_big_number, frontier, temp_vector, dist_type, true);


	// setting initial conditions for the loop
	current_path_id = 1;
	SkeletonGraph[seed_node_index].path_id = current_path_id;

	// Step 2.2 -- find local maxima in the labels.
	local_maxes_skeleton_hypotheses.clear();
	FindLocalMaximaSurfaceSingle(SkeletonGraph, bfs_label_node_to_node_double, already_explored,
			local_maxes_skeleton_hypotheses,R_model.number_voxels_per_dim,
			big_big_number, max_connectivity, this_voxel_already_treated);




	while (local_maxes_skeleton_hypotheses.size() > 0){
		loop_counter++;


		multiple_paths_from_same_start.clear();
		local_tip = local_maxes_skeleton_hypotheses[0];


		if (this_voxel_already_treated[local_tip] == false){ // determine if this voxel has been tested already.

			pre_skeleton_nodes.clear();

			local_max_distance_double = BFS_search_from_tip_to_existing_skeleton(SkeletonGraph, bfs_label_outside_to_inside_inverse_double,
					bfs_from_tip_double, big_big_number, local_tip, pre_skeleton_bool, stop_lines, frontier, dist_type);


			// make a copy ....this gets reset in a later step.
			for (int_type_t k = 0; k < number_nodes; k++){
				stop_lines_copy[k] = stop_lines[k];
			}

			pre_skel_distance_double.clear();


			keep_components.clear();
			IdentifyLocalMinimaFromPreSkeletonsSetInKeepUpdatedComponents(SkeletonGraph, pre_skeleton_bool, stop_lines, pre_skeleton_nodes, bfs_from_tip_double, keep_components);


			for (int_type_t k = 0; k < pre_skeleton_nodes.size(); k++){
				pre_skel_distance_double.push_back(pair<int_type_t, double>(pre_skeleton_nodes[k], bfs_from_tip_double[pre_skeleton_nodes[k]]));
			}

			// sort by arrival ....

			std::sort(pre_skel_distance_double.begin(), pre_skel_distance_double.end(), distance_comp_function_skel_distance_double);

			path_counter = 0;
			last_path_added = paths.size();

			pending_paths.clear();
			for (int_type_t k = 0; k < pre_skeleton_nodes.size(); k++){
				current_path.clear();
				linked_path.clear();

				WalkPreSkeletonToTipToFindJunctionsWithDeepIntoCornersModificationConnectivity(SkeletonGraph, bfs_from_tip_double, bfs_label_outside_to_inside_double,
						pre_skel_distance_double[k].first, current_path);

				FindLinks(SkeletonGraph, bfs_label_node_to_node_double, current_path[0], linked_path);

				linked_path[0] = linked_path[1];  // housekeeping, linked path is only two elements long, we only need the second one.
				linked_path.pop_back();
				linked_path.insert(linked_path.end(), current_path.begin(), current_path.end());

				pending_paths.push_back(linked_path);
			}

			// pre Loop handling -- step 4.2 -- deal with some noisy regions
			if (pre_skeleton_nodes.size() > 1){
				vector< vector<int_type_t> > temp_pending_paths;

				for (int_type_t k = 0; k < pre_skeleton_nodes.size(); k++){
					walks_through_another_stop_set0 = true;
					walks_through_another_stop_set1 = false;
					int_type_t out_of_first_stop_set_index = pending_paths[k].size();

					int_type_t limit = pending_paths[k].size();
					if (limit > 10){ limit = 10;}


					// first one is the skeleton.
					for (int_type_t i = 1, in  = pending_paths[k].size(); i < in && walks_through_another_stop_set0 == true; i++){
						if (stop_lines_copy[pending_paths[k][i]] == false){
							walks_through_another_stop_set0 = false;
							out_of_first_stop_set_index = i;
						}
					}

					for (int_type_t i = out_of_first_stop_set_index, in  = pending_paths[k].size(); i < in && walks_through_another_stop_set1 == false; i++){
						if (stop_lines_copy[pending_paths[k][i]] == true){
							walks_through_another_stop_set1 = true;
						}
					}

					if (walks_through_another_stop_set1 == false){
						temp_pending_paths.push_back(pending_paths[k]);
					}	else {
						ClearConnectedComponent(SkeletonGraph, pending_paths[k][out_of_first_stop_set_index - 1],  stop_lines_copy);
					}
				}

				temp_pending_paths.swap(pending_paths);

				// need to go through and zero out the extra stop set components as well.
			}

			if (pending_paths.size() == 0){
				cout << "We have no paths" << endl; exit(1);
			}


			// the size of pending paths may have changed, here we have the true loops.
			// Step 4.2 -- loop handling
			if (pending_paths.size() > 1){ ///loop case
				int_type_t first_intersection = 0;
				int_type_t first_path_max_index_of_intersection = 0;

				int_type_t first_path_id = current_path_id;
				int_type_t last_path_id = 0;

				// this changes the labels to be the maximum numbered-branch that is coicident with the current branch.
				// the last branch has all the same label.
				for (int_type_t k = 0; k < pending_paths.size(); k++){
					// don't overwrite the ends ....
					first_intersection = 0;
					for (int_type_t j = 1, jn = pending_paths[k].size(); j < jn && first_intersection == 0; j++){
						SkeletonGraph[pending_paths[k][j]].path_id = current_path_id;
					}

					if (SkeletonGraph[pending_paths[k][0]].path_id == 0){
						SkeletonGraph[pending_paths[k][0]].path_id = current_path_id;
					}

					path_counter++;
					current_path_id++;
				}

				last_path_id = current_path_id - 1;
				bool stop_descent = false;
				int_type_t path_id_for_this_one;

				for (int_type_t k = 0; k < pending_paths.size() - 1; k++){
					path_id_for_this_one = k + first_path_id;

					stop_descent = false;
					for (int_type_t i = pending_paths[k].size(); i > 2 && stop_descent == false; i--){
						if (SkeletonGraph[pending_paths[k][i-2]].path_id == path_id_for_this_one){
							stop_descent = true;
						}	else {
							pending_paths[k].pop_back();
						}
					}
					paths.push_back(pending_paths[k]);
				}

				path_id_for_this_one = pending_paths.size() - 1 + first_path_id;
				for (int_type_t k = 0; k < pending_paths.size() - 1; k++){
					if (SkeletonGraph[pending_paths[k].back()].path_id == path_id_for_this_one){
						SkeletonGraph[pending_paths[k].back()].path_id = 0;
					}
				}

				// pop off the 'tail' for this local max
				stop_descent = false;
				for (int_type_t i = pending_paths.back().size() + 1; i > 1 && stop_descent == false; i--){
					if (SkeletonGraph[pending_paths[pending_paths.size() - 1][i-2]].path_id == 0){
						stop_descent = true;
					}	else {
						SkeletonGraph[pending_paths[pending_paths.size() - 1].back()].path_id = 0;
						pending_paths[pending_paths.size() - 1].pop_back();
					}
				}

				// then rewrite
				for (int_type_t k = 0; k < pending_paths.size() - 1; k++){
					if (SkeletonGraph[pending_paths[k].back()].path_id == 0){
						SkeletonGraph[pending_paths[k].back()].path_id = path_id_for_this_one;
					}
				}

				//cout << "Line 910 " << endl;
				for (int_type_t k = 0; k < pending_paths.size(); k++){
					paths.push_back(pending_paths[k]);
				}
			}	else {
				// mark the tip
				this_voxel_already_treated[local_tip] = true;
				// non loop case

				if (current_path_id > 1){


					passed_spurious_tests = ClassifyUsingChiSquaredAssumption(SkeletonGraph, pending_paths[0][0], stop_lines_copy, paths, R_model.number_voxels_per_dim,
							big_big_number, max_connectivity, pending_paths[0], threshold);

				}	else {
					passed_spurious_tests = true;
				}

				if (passed_spurious_tests){
					//	cout << "Entering 1394" << pending_paths.size() << endl;
					// this takes care of the case with more than one segment per tip

					for (int_type_t j = 1, jn = pending_paths[0].size(); j < jn; j++){
						SkeletonGraph[pending_paths[0][j]].path_id = current_path_id;
					}
					if (SkeletonGraph[pending_paths[0][0]].path_id == 0){
						SkeletonGraph[pending_paths[0][0]].path_id = current_path_id;
					}

					path_counter++;

					current_path_id++;
					paths.push_back(pending_paths[0]);
				}	else {

					FillInTreatedMap(SkeletonGraph, bfs_from_tip_double, this_voxel_already_treated, big_big_number);
				}
			}
			current_path.clear();

			// Step 2.1 -- update BFS1
			if (paths.size() != last_path_added){
				vector<int_type_t> all_paths = paths[last_path_added];

				for (int_type_t new_path = last_path_added + 1; new_path < paths.size(); new_path++){
					all_paths.insert(all_paths.end(), paths[new_path].begin(), paths[new_path].end());
				}
				greatest_mod_distance_double = BFS_search_from_start_set(SkeletonGraph,
						bfs_label_outside_to_inside_inverse_double, bfs_label_node_to_node_double,
						big_big_number, frontier, all_paths, dist_type, false);
			}
		}

		// Step 2.2 -- locate potential endpoints in BFS1 map
		// find a new local maxima for the next round .....
		local_maxes_skeleton_hypotheses.clear();

		FindLocalMaximaSurfaceSingle(SkeletonGraph, bfs_label_node_to_node_double, already_explored,
				local_maxes_skeleton_hypotheses,R_model.number_voxels_per_dim,
				big_big_number, max_connectivity, this_voxel_already_treated);
	}
	out.close();

	// stop timer, start writing.
	auto t1 = std::chrono::high_resolution_clock::now();
	cout << "Time for skeletonization ... " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()<< " milliseconds "<< endl;

	filename = write_directory + "details.txt";
	cout << "Filename " << filename << endl;
	out.open(filename.c_str());
	out << "time " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1000.0 << " seconds "<< endl;
	out << "Number of occupied voxels " << number_nodes << endl;
	out << "Number possible voxels " << R_model.number_voxels_grid << endl;
	out << "Number of proposed tips " << loop_counter << endl;
	out << "Threshold for spurious segment classification " << threshold << endl;
	out.close();

	for (int_type_t i = 0, in = paths.size(); i < in; i++){
		vector<double> temp(paths[i].size(), R_model.inter_voxel_distance);
		path_diameters_single.push_back(temp);
	}

	for (int_type_t i = 0, in = paths.size(); i < in; i++){
		vector<double> temp(paths[i].size(), R_model.inter_voxel_distance);
		path_diameters.push_back(temp);
	}



	filename = write_directory + "paths_ECCV257.ply";

	WritePaths(R_model, SkeletonGraph, paths, path_diameters_single, filename );

	for (int_type_t i = 0, in = paths.size(); i < in; i++){
		for (int_type_t j = 0, jn= paths[i].size(); j < jn; j++){
			path_diameters[i][j] = 1 + bfs_label_outside_to_inside_double[paths[i][j]]*R_model.inter_voxel_distance;
		}
	}

	for (int_type_t k = 0; k < skeleton.size(); k++){
		skeleton[k] = true;
	}
	for (int_type_t i = 0, in = paths.size(); i < in; i++){
		for (int_type_t j = 0, jn= paths[i].size(); j < jn; j++){
			skeleton[SkeletonGraph[paths[i][j]].grid_id] = false;
		}
	}


	R_model.GenerateAndWriteSurfaceInPlyFormat(write_directory, 0, "initial", NULL, false);
	R_model.configuration_grid = skeleton;

	R_model.GenerateAndWriteSurfaceInPlyFormat(write_directory, 0, "skel_ECCV257_", NULL, false);

	string write_file = write_directory + "result.txt";
	out.open(write_file.c_str());

	for (int_type_t i = 0; i < R_model.number_voxels_grid; i++){
		if (skeleton[i] == false){
			out << i << " ";
		}

	}
	out << "-1" << endl;

	out.close();


	delete [] bfs_label_outside_to_inside_double;
	delete [] bfs_label_outside_to_inside_inverse_double;
	delete [] bfs_label_node_to_node_double;
	delete [] bfs_from_tip_double;

	delete [] frontier;
	delete [] already_explored;
	delete [] this_voxel_already_treated;
	delete [] pre_skeleton_bool;
	delete [] stop_lines;
}


