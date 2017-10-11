#ifndef SKELGRAPH_HPP_
#define SKELGRAPH_HPP_


#include "Includes.hpp"
#include "ReconstructionStructure.hpp"


#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <list>
#include <ctime>
#include <time.h>

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <inttypes.h>

bool IsIntType(DISTANCE d);

class SkelGraph{
public:
	int_type_t voxel_id; // sparse numbering
	int_type_t path_id; // for skeleton and path
	int_type_t grid_id; // in the voxel grid
	int_type_t cc_id;
	vector<SkelGraph*> neighbors;
	vector<int8_t> neighbor_distance; // this can be much smaller ....
	vector<double> neighbor_distanced; // this can be much smaller ....

	SkelGraph();

	SkelGraph(int_type_t vi, int_type_t gi);

	~SkelGraph();

	bool IsLocalMaxima(vector<int_type_t>& labels);

	void Print();
};

void CreateSurfaceGraphFromGrid(ReconstructionStructure& RS, vector<SkelGraph>& SG);

void WritePaths(ReconstructionStructure& RS, vector<SkelGraph>& SG, vector< vector<int_type_t> >& paths, vector< vector<double> >& path_diameters,
		string outfile);


int_type_t SquaredDifferenceOfUints(int_type_t A, int_type_t B);

void ReworkUsingOnlyBiggestConnectedComponent(ReconstructionStructure& RS, vector<SkelGraph>& SG, vector<SkelGraph>& newSG, int_type_t biggest_connected_component);

pair<int_type_t, int_type_t> CreateGraphFromGridReturnPair(ReconstructionStructure& RS, vector<SkelGraph>& SG, vector<int_type_t>& cc_counts, DISTANCE dist_type = L1);

bool distance_comp_function_skel_distance_double(pair<int_type_t, double> a, pair<int_type_t, double> b);

int Return26ConnectedNeighbors(ReconstructionStructure& RS, int_type_t start_voxel, int_type_t* member_array, int_type_t* distance);

#endif /* SKELGRAPH_HPP_ */
