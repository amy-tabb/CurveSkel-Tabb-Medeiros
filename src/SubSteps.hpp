#ifndef SUBSTEPS_HPP_
#define SUBSTEPS_HPP_


#include "SkelGraph.hpp"


void SortSetsWithStops(vector<SkelGraph>& SG, double* bfs_labels,  vector<int_type_t>& F_set, vector<int_type_t>& F1_prime, vector<int_type_t>& F2_prime,
		double l_min_final, bool* frontier,  DISTANCE dist_type);

int_type_t  FindALocalMaximaForSeed(vector<SkelGraph>& SG, double* bfs_labels, int_type_t big_number);


double  BFS_general_outside_to_inside_master_induction1(vector<SkelGraph>& SG, double* bfs_labels, vector<int_type_t>& grid_structure, int_type_t big_number,
		int_type_t connectivity, bool* frontier_map, DISTANCE dist_type, bool parallel = false);


int_type_t  BFS_search_from_start_set(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse,
		double* bfs_labels, int_type_t big_number,
		bool* frontier_map, vector<int_type_t>& start_set, DISTANCE dist_type, bool initialize);


template <class T>
void OverwriteFirstArrayWithSecondParallel(T* arr0, T* arr1, int_type_t N){

#pragma omp parallel for schedule(static, 10)
	for (int_type_t i = 0; i < N; i++){
		arr0[i] = arr1[i];
	}

}

template <class T>
void OverwriteFirstArrayWithSecond(T* arr0, T* arr1, int_type_t N){

	for (int_type_t i = 0; i < N; i++){
		arr0[i] = arr1[i];
	}

}

void FindLocalMaximaSurfaceSingle(vector<SkelGraph>& SG, double* bfs_labels, bool* already_explored,
		vector<int_type_t>& local_maxima,
		vector<int_type_t>& grid_structure, int_type_t big_number, int max_connectivity, bool* already_treated);


int_type_t  BFS_search_from_tip_to_existing_skeleton(vector<SkelGraph>& SG, double* bfs_labels_oTi_inverse,
		double* bfs_labels, int_type_t big_number, int_type_t start_index,
		bool* pre_skeleton, bool* stop_lines, bool* frontier_map, DISTANCE dist_type);

void FillInTreatedMap(vector<SkelGraph>& SG, double* previous_bfs_map, bool* already_treated,
		int_type_t big_number);

int_type_t  ClearConnectedComponent(vector<SkelGraph>& SG, int_type_t start_index, bool* stop_lines);



void IdentifyLocalMinimaFromPreSkeletonsSetInKeepUpdatedComponents(vector<SkelGraph>& SG, bool* in_pre_skeleton, bool* updated_by_skeleton,
		vector<int_type_t>& pre_skeleton,
		double* bfs_map,
		vector< vector<int_type_t> >& keep_components);

void FindLinks(vector<SkelGraph>& SG, double* bfs_labels, int_type_t tip_index,
		vector<int_type_t>& path_contents);

double FindMaxInBFSMap(vector<int_type_t>& set_to_consider,	double* bfs_map);

void ReturnXYZFromIndicesSkelGraph(SkelGraph& S, int_type_t& x, int_type_t& y, int_type_t& z, vector<int_type_t>& grid_structure);

void WalkPreSkeletonToTipToFindJunctionsWithDeepIntoCornersModificationConnectivity(vector<SkelGraph>& SG, double* bfs_labels, double* bfs_labels_dis,
		int_type_t tip_index,
		vector<int_type_t>& path_contents);

bool ClassifyUsingChiSquaredAssumption(vector<SkelGraph>& SG, int_type_t skeleton_index, bool* stop_lines,
		vector<vector<int_type_t> >& paths, vector<int_type_t>& grid_structure,
		int_type_t big_big_number, int_type_t connectivity, vector<int_type_t>& proposed_path, double probability_threshold);


#endif /* SUBSTEPS_HPP_ */
