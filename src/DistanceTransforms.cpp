#include "SkelGraph.hpp"
#include "DistanceTransforms.hpp"
#include "SubSteps.hpp"

double f(int_type_t x, int_type_t i, int_type_t y, vector<vector<double> >& g){
	return pow((double(x) - double(i)), 2) + pow(g[i][y], 2);
}

double fwith_square(int_type_t x, int_type_t i, int_type_t y, vector<vector<double> >& g){

	return pow((double(x) - double(i)), 2) + pow(g[i][y], 2);
}

// for the non-grid version, scan
double f(int_type_t x, int_type_t i, double* g, int_type_t* line){

	return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

// for the non-grid version, scan
double f_with_square(int_type_t x, int_type_t i, double* g, int_type_t* line){

	return pow((double(x) - double(i)), 2) + g[line[i]];
}

// for the non-grid version, scan
double f(int_type_t x, int_type_t i, double* g, int_type_t* line, int_type_t big_number){

	if (line[i] == big_number){
		// this value is going to be zero.
		return pow((double(x) - double(i)), 2);
	}	else {
		return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
	}
}



// NOte: here g is the squared distance
double f(int_type_t z, int_type_t i, int_type_t x, int_type_t y, vector<vector<vector<double> > >& g){
	return pow((double(z) - double(i)), 2) + g[x][y][i];
}

double Sep(int_type_t i, int_type_t u, int_type_t y, vector<vector<double> >& g){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// scan 4, non-grid version
double Sep(int_type_t i, int_type_t u, double* g, int_type_t* line){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + pow(g[line[u]], 2) - pow(g[line[i]], 2))/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// scan 4, non-grid version
double Sep_with_square(int_type_t i, int_type_t u, double* g, int_type_t* line){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + g[line[u]] - g[line[i]])/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// scan 4, non-grid version
double Sep(int_type_t i, int_type_t u, double* g, int_type_t* line, int_type_t big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double a= 0;
	double b = 0;
	if (line[u] != big_number){
		a = pow(g[line[u]], 2);
	}

	if (line[i] != big_number){
		b = pow(g[line[i]], 2);
	}

	double s = int(double(u*u) - double(i*i) + a - b)/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// NOte: here g is the squared distance
double Sep(int_type_t i, int_type_t u, int_type_t x, int_type_t y, vector<vector<vector<double> > >& g){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + g[x][y][u] - g[x][y][i])/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}




SkelGraph* FindNextVoxelWithOnAxisWithDirection(SkelGraph& S, int axis, int direction, vector<int_type_t>& grid_structure){
	//bool has_empty = false;
	//int_type_t temp_index,
	int_type_t x, y, z, x0, y0, z0;

	//if (S.neighbors.size() != connectivity)
	//{


	ReturnXYZFromIndicesSkelGraph(S, x, y, z, grid_structure);

	for (int_type_t i = 0; i < S.neighbors.size(); i++){
		ReturnXYZFromIndicesSkelGraph(*S.neighbors[i], x0, y0, z0, grid_structure);

		switch (axis) {
		case 0: {
			if (direction == 1 && x0 == x + 1 && y0 == y && z0 == z){
				return S.neighbors[i];
			}

			if (direction == -1 && x0 + 1 == x && y0 == y && z0 == z){
				return S.neighbors[i];;
			}
		} break;
		case 1: {
			if (direction == 1 && y0 == y + 1 && x0 == x && z0 == z){
				return S.neighbors[i];
			}

			if (direction == -1 && y0 + 1 == y && x0 == x && z0 == z){
				return S.neighbors[i];
			}
		} break;
		case 2: {
			if (direction == 1 && z0 == z + 1 && x0 == x && y0 == y){
				return S.neighbors[i];
			}

			if (direction == -1 && z0 + 1 == z && x0 == x && y0 == y){
				return S.neighbors[i];
			}
		} break;
		default: {
			cout << "Bad axis option, Find Next voxel " << endl; exit(1);
		} break;
		}
	}

	return 0; /// no match found, so this one has an empty in the direction indicated.
	//}


}

int_type_t FindNextVoxelWithOnAxisWithDirectionReturnInt(SkelGraph& S, int axis, int direction, vector<int_type_t>& grid_structure){
	//bool has_empty = false;
	int_type_t temp_index, x, y, z, x0, y0, z0;

	//if (S.neighbors.size() != connectivity)
	//{


	ReturnXYZFromIndicesSkelGraph(S, x, y, z, grid_structure);

	for (int_type_t i = 0; i < S.neighbors.size(); i++){
		ReturnXYZFromIndicesSkelGraph(*S.neighbors[i], x0, y0, z0, grid_structure);

		switch (axis) {
		case 0: {
			if (direction == 1 && x0 == x + 1 && y0 == y && z0 == z){
				return i;
			}

			if (direction == -1 && x0 + 1 == x && y0 == y && z0 == z){
				return i;
			}
		} break;
		case 1: {
			if (direction == 1 && y0 == y + 1 && x0 == x && z0 == z){
				return i;
			}

			if (direction == -1 && y0 + 1 == y && x0 == x && z0 == z){
				return i;
			}
		} break;
		case 2: {
			if (direction == 1 && z0 == z + 1 && x0 == x && y0 == y){
				return i;
			}

			if (direction == -1 && z0 + 1 == z && x0 == x && y0 == y){
				return i;
			}
		} break;
		default: {
			cout << "Bad axis option, Find Next voxel " << endl; exit(1);
		} break;
		}
	}

	// max number of neighbors in almost any situation
	return 26; /// no match found, so this one has an empty in the direction indicated.
	//}


}


void DistanceTransform1Parallel(vector<SkelGraph>& SG, double* dt_xyz, vector<int_type_t>& grid_structure, double big_number, int_type_t connectivity){


	// This is an implementation of the Meijster linear-time distance transform method, but for connected components instead of a grid
	// within each segment, the distances are going to be constrained to be within the object.  Need some scratch surfaces ....

	int_type_t big_number_int = big_number;
	int_type_t xsize, ysize, zsize;
	int_type_t n = SG.size();

	xsize= grid_structure[0];
	ysize= grid_structure[1];
	zsize= grid_structure[2];

	// we need s and t for each thread ....arrays are faster
	int_type_t number_threads = omp_get_max_threads();
	cout << "Number threads " << number_threads << endl;

//	// AKA don't do hyperthreading.
//	if (number_threads > 24){
//		omp_set_num_threads(24);
//		number_threads = 24;
//	}


	vector< int_type_t* > s(number_threads, 0);
	vector< int_type_t* > t(number_threads, 0);
	vector< int_type_t* > current_line(number_threads, 0);


	for (int_type_t i = 0; i < number_threads; i++){
		s[i] = new int_type_t[max(xsize, zsize)];
		t[i] = new int_type_t[max(xsize, zsize)];
		current_line[i] = new int_type_t[n + 2];
	}
	int q; int w;


	// these can really be ints, b/c we're not doin g sqrts until the end ....
	double* g_y = new double[n + 1];
	double* g_xy = new double[n + 1];

	int_type_t thread_id = 0;

	for (int_type_t i = 0; i < n; i++){
		g_y[i] = big_number;

	}

	g_y[n] = 0;
	g_xy[n] = 0;

	// first, go fishing.  Do scan 1 and 2 where x is static, z static, y increases and decreases.
	// scan 1
	SkelGraph* Sptr = 0;
	int current_distance;

#pragma omp parallel for private(Sptr, current_distance)
	for (int_type_t i = 0; i < n; i++){
		if (SG[i].neighbors.size() != connectivity){// any candidates will have less than the number of neighbors ....

			// is this one have an empty y neighbor?
			Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 1, 1, grid_structure);

			if (Sptr == 0){
				Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 1, -1, grid_structure);

				current_distance = 1;
				g_y[i] = current_distance;
				while (Sptr != 0){// keep on walking ....
					current_distance++;
					g_y[Sptr->voxel_id] = current_distance;
					Sptr = FindNextVoxelWithOnAxisWithDirection(*Sptr, 1, -1, grid_structure);
				}
			}


		}
	}

	// scan 2
#pragma omp parallel for private(Sptr, current_distance)
	for (int_type_t i = 0; i < n; i++){
		if (SG[i].neighbors.size() != connectivity){// any candidates will have less than the number of neighbors ....

			Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 1, -1, grid_structure);

			if (Sptr == 0){
				Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 1, 1, grid_structure);

				current_distance = 1;
				g_y[i] = current_distance;
				while (Sptr != 0){// keep on walking ....
					current_distance++;
					if (g_y[Sptr->voxel_id] > current_distance){
						g_y[Sptr->voxel_id] = current_distance;
						Sptr = FindNextVoxelWithOnAxisWithDirection(*Sptr, 1, 1, grid_structure);
					}	else {
						Sptr = 0;
					}
				}
			}
		}
	}

	int_type_t line_counter = 0;
	// scan 3 and 4 -- y and z are static, x is moving.


#pragma omp parallel for private (Sptr, line_counter, w, q, thread_id)
	for (int_type_t i = 0; i < n; i++){
		if (SG[i].neighbors.size() != connectivity){// any candidates will have less than the number of neighbors .... and not currently be marked

			thread_id = omp_get_thread_num();
			//cout << "thread id " << thread_id << endl;
			Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 0, -1, grid_structure); // only do this is this voxel is on the edge and nothing to the negative side
			if (Sptr == 0)
			{

				current_line[thread_id][0] = n;
				current_line[thread_id][1] = i;
				line_counter = 2;

				Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 0, 1, grid_structure);
				while (Sptr != 0){
					current_line[thread_id][line_counter] = Sptr->voxel_id;
					line_counter++;

					Sptr = FindNextVoxelWithOnAxisWithDirection(*Sptr, 0, 1, grid_structure);
				}
				current_line[thread_id][0] = n;
				current_line[thread_id][line_counter] = n;
				line_counter++;

				q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;

				// scan 3
				for (int_type_t u = 1; u < line_counter; u++){
					while (q >= 0 && f(t[thread_id][q], s[thread_id][q], g_y, current_line[thread_id]) > f(t[thread_id][q], u, g_y, current_line[thread_id])){
						q--;
					}

					if (q < 0){
						//cout << "q is zero " << endl;
						q = 0;  s[thread_id][0] = u;
					} else {
						w = 1 + Sep(s[thread_id][q], u, g_y, current_line[thread_id]);
						//cout << "w " << w << "   " << u << endl;

						if (w < 0){
							cout << "EEEO! w less than zero " << endl;
							exit(1);
						}
						if (w < int(line_counter)){
							q++; s[thread_id][q] = u; t[thread_id][q] = w;
						}
					}
				}

				// scan 4
				for (int_type_t u = line_counter; u > 0; u--){

					if (current_line[thread_id][u-1] != big_number_int){
						g_xy[current_line[thread_id][u-1]] = f(u - 1, s[thread_id][q], g_y, current_line[thread_id]);
						//dt_xyz[current_line[u-1]] = sqrt(f(u - 1, s[q], g_y, current_line));
					}

					if ((u-1) == t[thread_id][q]){
						q--;
					}
				}

				//i = n;
			}
		}
	}


	// scan 5 and 6 -- x and y are static, z is moving.
#pragma omp parallel for private (Sptr, line_counter, w, q, thread_id)
	for (int_type_t i = 0; i < n; i++){
		if (SG[i].neighbors.size() != connectivity){// any candidates will have less than the number of neighbors ....

			thread_id = omp_get_thread_num();
			Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 2, -1, grid_structure); // only do this is this voxel is on the edge and nothing to the negative side

			if (Sptr == 0)
			{
				current_line[thread_id][0] = n;
				current_line[thread_id][1] = i;
				line_counter = 2;

				Sptr = FindNextVoxelWithOnAxisWithDirection(SG[i], 2, 1, grid_structure);
				while (Sptr != 0){
					current_line[thread_id][line_counter] = Sptr->voxel_id;
					line_counter++;

					Sptr = FindNextVoxelWithOnAxisWithDirection(*Sptr, 2, 1, grid_structure);
				}
				current_line[thread_id][line_counter] = n;
				line_counter++;



				q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;

				// debug -- do only one line and the quit.

				// scan 3
				for (int_type_t u = 1; u < line_counter; u++){

					while (q >= 0 && f_with_square(t[thread_id][q], s[thread_id][q], g_xy, current_line[thread_id]) > f_with_square(t[thread_id][q], u, g_xy, current_line[thread_id])){
						q--;
					}

					if (q < 0){
						//cout << "q is zero " << endl;
						q = 0;  s[thread_id][0] = u;
					} else {
						w = 1 + Sep_with_square(s[thread_id][q], u, g_xy, current_line[thread_id]);
						//cout << "w " << w << "   " << u << endl;

						if (w < 0){
							cout << "EEEO! w less than zero " << endl;
							exit(1);
						}
						if (w < line_counter){
							q++; s[thread_id][q] = u; t[thread_id][q] = w;
						}
					}
				}

				// scan 4 -- both ends of current line are n and therefore should not be used
				for (int_type_t u = line_counter; u > 0; u--){

					if (current_line[thread_id][u-1] != n){
						dt_xyz[current_line[thread_id][u-1]] = sqrt(f_with_square(u - 1, s[thread_id][q], g_xy, current_line[thread_id]));
					}

					// still need to do these updates within the context of the Meijster algo.
					if ((u-1) == t[thread_id][q]){
						q--;
					}
				}

			}
		}
	}

	delete [] g_y;
	delete [] g_xy;

	for (int_type_t i = 0; i < number_threads; i++){
		delete [] s[i];
		delete [] t[i];
		delete [] current_line[i];
	}

}


