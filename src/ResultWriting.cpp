#include "ResultWriting.hpp"
#include "ReconstructionStructure.hpp"


void ReconstructionStructure::WritePlyFile(string outfile,
		vector<int_type_t>& subdims,
		vector<double>& points,
		vector<int_type_t>& faces,
		vector<int_type_t>& color){

	// each vertex needs a color ....

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << points.size()/3 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << faces.size()/2 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	out << "end_header" << endl;


	for (int_type_t i = 0; i < (points.size()/3); i++){
		out << points[3*i] << " " <<  points[3*i + 1]  << " " <<  points[3*i + 2]  << " ";
		out << color[0] << " " <<  color[1] << " " << color[2] << " 255" << endl;
	}

	for (int_type_t i = 0; i < (faces.size()/4); i++){
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 1]   << " " <<  faces[4*i + 2]  << endl; //<< " " << faces[4*i + 3] << endl;;
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 2]   << " " <<  faces[4*i + 3]  << endl;
	}

	out << endl;

	out.close();
}


void AxisAngleToDCM(vector<vector<double> >& R, vector<double>& axis, float theta){
	// assume axis normalized


	double xs, ys, zs, xC, yC, zC, xyC, yzC, zxC;
	double x, y, z;
	double c, s, C;
	x = axis[0];  y = axis[1];  z = axis[2];
	c = cos(theta);  s = sin(theta); C = 1-c;
	xs = x*s;  ys = y*s; zs = z*s;
	xC = x*C; yC = y*C; zC = z*C;
	xyC = x*yC; yzC = y*zC; zxC = z*xC;
	R[0][0] = x*xC+c;
	R[0][1] = xyC-zs;
	R[0][2] = zxC+ys;

	R[1][0] = xyC+zs;
	R[1][1] = y*yC+c;
	R[1][2] = yzC-xs;

	R[2][0]=  zxC-ys;
	R[2][1] = yzC+xs;
	R[2][2] = z*zC+c;
}

void PrintVector(vector<double>& p){

	for (int i = 0; i < 3; i++){
		cout << p[i] << " ";
	}
	cout << endl;
}

void PrintMatrix(vector< vector<double> >& p){

	for (int i = 0; i < int(p.size()); i++){
		for (int j = 0; j < int(p[i].size()); j++){
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void CreatePlaneWithTheseThreePoints(vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p){
	vector<double> v0(3, 0);
	vector<double> v1(3, 0);
	vector<double> n(3, 0);

	SubtractVectorFromVector(p1, p0, v0, 3);
	NormalizeVector(v0);

	SubtractVectorFromVector(p1, p2, v1, 3);
	NormalizeVector(v1);

	CrossProduct(v1, v0, n);
	NormalizeVector(n);

	p[0]= n[0];
	p[1] = n[1];
	p[2] = n[2];
	p[3] = -(DotProduct(p1, n, 3));


}

bool IsOnBackSideOfPlane(vector<double>& plane, vector<double>& point){

	double d = DotProduct(plane, point, 3) + plane[3];

	return d < 0;
}

void WritePaths(ReconstructionStructure& RS, vector<SkelGraph>& SG, vector< vector<int_type_t> >& paths, vector< vector<double> >& path_diameters,
		string outfile){

	char ch;

	vector<vector<int> > colors;

	vector<int> c(3);
	// lowest value is grey
	c[0] = 100;
	c[1] = 100;
	c[2] = 100;

	colors.push_back(c);

	// next lowest value is purple
	c[0] = 128;
	c[1] = 0;
	c[2] = 128;
	colors.push_back(c);

	// next lowest value is blue
	c[0] = 0;
	c[1] = 0;
	c[2] = 200;
	colors.push_back(c);

	// next lowest value is cyan
	c[0] = 0;
	c[1] = 255;
	c[2] = 255;
	colors.push_back(c);

	// next lowest value is green
	c[0] = 0;
	c[1] = 255;
	c[2] = 0;
	colors.push_back(c);

	// next lowest value is yellow
	c[0] = 255;
	c[1] = 255;
	c[2] = 0;
	colors.push_back(c);

	// next lowest value is red
	c[0] = 255;
	c[1] = 0;
	c[2] = 0;
	colors.push_back(c);



	vector<int_type_t> faces;
	vector<double> points;
	vector<int_type_t> current_face(4, 0);
	int number_divisions = 16;
	vector<int_type_t> start_circle; // these are point indices
	vector<int_type_t> middle_circle;
	vector<double> direction_vector(3, 0);
	vector<double> start_point(3, 0);
	vector<double> middle_point(3, 0);
	vector<double> next_point(3, 0);
	int_type_t number_points = 0;
	int_type_t x_index, y_index, z_index;
	int_type_t color_group;
	vector<double> v0(3, 0);
	vector<double> v1(3, 0);
	vector<double> plane_normal(3, 0);
	vector<double> neg_plane_normal(3, 0);
	vector<double> rotation_axis(3, 0);
	vector<double> z_axis(3, 0);
	vector< vector<double> > R0(3, vector<double>(3));
	vector< vector<double> > R1(3, vector<double>(3));
	vector< vector<double> > R(3, vector<double>(3));

	vector<double> result0(3, 0);
	vector<double> result1(3, 0);
	vector<double> circle_point(3, 0);
	vector<double> circle_point_placed(3, 0);
	vector<double> radius_vector(3, 0);
	vector<double> diff_vector(3, 0);
	vector<int_type_t> temp_circle;
	vector<double> p0(3, 0);
	vector<double> p1(3, 0);
	vector<double> p2(3, 0);
	vector<double> facet_plane(4, 0);

	z_axis[2] = 1;
	double d0, d1;
	double radius;
	double theta;
	double sd;


	vector< int_type_t > path_by_point;


	for (int_type_t i = 0, in = paths.size(); i < in; i++){

		path_by_point.push_back(number_points);

		start_circle.clear();
		middle_circle.clear();
		if (paths[i].size() > 0){
			for (int_type_t j = 1, jn= paths[i].size() - 1; j < jn; j++){

				if (j == 1){
					RS.ReturnXYZIndicesFromIndex(SG[paths[i][j-1]].grid_id, x_index, y_index, z_index);
					start_point[0] = RS.BB[0][0] + (x_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
					start_point[1] = RS.BB[0][1] + (y_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
					start_point[2] = RS.BB[0][2] + (z_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;

					// first point
					points.insert(points.end(), start_point.begin(), start_point.end());
					start_circle.push_back(number_points);
					number_points++;

					RS.ReturnXYZIndicesFromIndex(SG[paths[i][j]].grid_id, x_index, y_index, z_index);
					middle_point[0] = RS.BB[0][0] + (x_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
					middle_point[1] = RS.BB[0][1] + (y_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
					middle_point[2] = RS.BB[0][2] + (z_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;

					SubtractVectorFromVector(middle_point, start_point, v0, 3);
					NormalizeVector(v0);
				}


				RS.ReturnXYZIndicesFromIndex(SG[paths[i][j + 1]].grid_id, x_index, y_index, z_index);
				next_point[0] = RS.BB[0][0] + (x_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
				next_point[1] = RS.BB[0][1] + (y_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;
				next_point[2] = RS.BB[0][2] + (z_index)*RS.inter_voxel_distance + RS.inter_voxel_distance/2.0;

				// we need all three to create the middle circle ... first we need to determine the plane's normal.
				SubtractVectorFromVector(next_point, middle_point,v1, 3);
				NormalizeVector(v1);

				AddVectorToVector(v0, v1, plane_normal, 3);
				NormalizeVector(plane_normal);

				middle_circle.clear();
				radius = path_diameters[i][j];

				neg_plane_normal[0] = -plane_normal[0];
				neg_plane_normal[1] = -plane_normal[1];
				neg_plane_normal[2] = -plane_normal[2];

				if (SquaredDistance(z_axis, plane_normal, 3) < 0.1 || SquaredDistance(z_axis, neg_plane_normal, 3) < 0.1){


					R[0][0] = 1;  R[0][1] = 0;  R[0][2] = 0;
					R[1][0] = 0;  R[1][1] = 1;  R[1][2] = 0;
					R[2][0] = 0;  R[2][1] = 0;  R[2][2] = 1;

					if (SquaredDistance(z_axis, neg_plane_normal, 3) < 0.1){
						R[0][0] = 1;  R[0][1] = 0;
						R[1][0] = 0; R[1][1] = -1;
						R[2][2] = 1;

					}

				}	else {

					CrossProduct(z_axis, plane_normal, rotation_axis);
					NormalizeVector(rotation_axis);

					theta = acos(DotProduct(z_axis, plane_normal, 3));

					AxisAngleToDCM(R0, rotation_axis, theta);
					AxisAngleToDCM(R1, rotation_axis, -theta);

					MultiplyMatrixVector(R0, 3, 3, z_axis, result0);
					MultiplyMatrixVector(R1, 3, 3, z_axis, result1);

					d0 = SquaredDistance(plane_normal, result0, 3);
					d1 = SquaredDistance(plane_normal, result1, 3);
					if (d0 < d1){
						R = R0;
					}	else {
						R = R1;
					}
				}

				radius_vector[0] = radius; radius_vector[1] = 0; radius_vector[2] = 0;
				// now compute points for this current circle using r ....
				for (double k = 0; k < number_divisions; k++){

					theta = (k/double(number_divisions)) * 3.14 * 2.0;

					radius_vector[0] = radius*cos(theta);
					radius_vector[1] = radius*sin(theta);

					// rotate the radius vector to a point ....
					MultiplyMatrixVector(R, 3, 3, radius_vector, circle_point);

					// now, center this on the middle point ...
					AddVectorToVector(circle_point, middle_point, circle_point_placed, 3);

					points.insert(points.end(), circle_point_placed.begin(), circle_point_placed.end());
					middle_circle.push_back(number_points);
					number_points++;



				}

				if (std::isnan(points.back())){
					cout << "a point " << middle_point[0] << ", " << middle_point[1] << ", " << middle_point[2] << endl;

					cout << "theta  " << theta << endl;
					cout << "plane normal ";  PrintVector(plane_normal);
					cout << "rotation axis ";  PrintVector(rotation_axis);
					cout << "result0 ";  PrintVector(result0);
					cout << "result1 ";  PrintVector(result1);
					cout << "d0 " << d0 << endl;
					cout << "d1 " << d1 << endl;

					cin >> ch;

				}


				// now make facets ....

				if (j == 1){
					// special case .... only one point in the start
					for (int_type_t k = 0, kn = middle_circle.size() - 1; k < kn; k++){
						// make one facet for each pair
						faces.push_back(middle_circle[k]);
						faces.push_back(start_circle[0]);
						faces.push_back(middle_circle[k + 1]);
					}

					faces.push_back(middle_circle[middle_circle.size() - 1]);
					faces.push_back(start_circle[0]);
					faces.push_back(middle_circle[0]);
				}	else {


					bool reversed = false;
					// first find the starting pair ....
					int best_match = start_circle.size();
					double best_distance = SquaredDistance(start_point, middle_point, 3) * 10.0;

					for (int_type_t k = 0, kn = middle_circle.size(); k < kn; k++){

						sd = pow(points[3*middle_circle[0]] - points[3*start_circle[k]], 2) +
								pow(points[3*middle_circle[0] + 1] - points[3*start_circle[k] + 1], 2) +
								pow(points[3*middle_circle[0] + 2] - points[3*start_circle[k] + 2], 2);

						p0[0] = points[3*middle_circle[0]];  p0[1] = points[3*middle_circle[0] + 1]; p0[2] = points[3*middle_circle[0] + 2];
						p1[0] = points[3*start_circle[k]];  p1[1] = points[3*start_circle[k] + 1]; p1[2] = points[3*start_circle[k] + 2];
						p2[0] = points[3*middle_circle[1]];  p2[1] = points[3*middle_circle[1] + 1]; p2[2] = points[3*middle_circle[1] + 2];

						if (sd < best_distance){

							// now, we need to test whether for any configuration both the start point and middle point and BEHIND the plane
							CreatePlaneWithTheseThreePoints(p0, p1, p2, facet_plane);

							bool back0 = IsOnBackSideOfPlane(facet_plane, start_point);
							bool back1 = IsOnBackSideOfPlane(facet_plane, middle_point);

							// id xOR back0 back1, this facet splits the ray between the two center points.
							if (!(back0 ^ back1)){

								best_distance = sd;
								best_match = k;

								if (!back0){
									reversed = true;
								}	else {
									reversed = false;
								}

							}
						}
					}


					{

						// reversed doesn't really make a difference, as much as splitting the place does, b/c we have two-sided faces.
						for (int_type_t k = 0, kn = middle_circle.size(); k < kn; k++){

							// make one facet for each pair
							faces.push_back(middle_circle[k % number_divisions]);
							faces.push_back(start_circle[(best_match + k) % number_divisions]);
							faces.push_back(middle_circle[(k + 1) % number_divisions]);

							// then for the start circle ...
							faces.push_back(start_circle[(best_match + k + 1) % number_divisions]);
							faces.push_back(middle_circle[(k + 1) % number_divisions]);
							faces.push_back(start_circle[(best_match + k) % number_divisions]);

						}

					}


					if (j == paths[i].size() - 2) {
						//cout << "Tidy up at the end" << endl;

						points.insert(points.end(), next_point.begin(), next_point.end());
						start_circle.clear();
						start_circle.push_back(number_points);
						number_points++;

						for (int_type_t k = 0, kn = middle_circle.size() - 1; k < kn; k++){
							// make one facet for each pair
							faces.push_back(middle_circle[k + 1]);
							faces.push_back(start_circle[0]);
							faces.push_back(middle_circle[k]);


						}

						faces.push_back(middle_circle[0]);
						faces.push_back(start_circle[0]);
						faces.push_back(middle_circle[middle_circle.size() - 1]);
					}


				}


				// copy over
				// start and middle points
				// vectors
				// circles
				start_point = middle_point;
				middle_point = next_point;
				v0 = v1;
				start_circle = middle_circle;
				// later -- tidy up? paths are from tip to source
			}
		}
	}
	path_by_point.push_back(number_points);

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << points.size()/3 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << faces.size()/3 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	out << "end_header" << endl;
	color_group = 0;

	for (int_type_t i = 0; i < (points.size()/3); i++){
		if (i >= path_by_point[color_group + 1]){
			color_group++;
		}
		out << points[3*i] << " " <<  points[3*i + 1]  << " " <<  points[3*i + 2]  << " ";
		out << colors[color_group % colors.size()][0] << " "
				<< colors[color_group % colors.size()][1] << " "
				<< colors[color_group % colors.size()][2] << "  255" << endl;
	}

	for (int_type_t i = 0; i < (faces.size()/3); i++){
		out << "3 " << faces[3*i]  << " " <<  faces[3*i + 1]   << " " <<  faces[3*i + 2]  << endl; //<< " " << faces[4*i + 3] << endl;;
	}

	out << endl;

	out.close();
}



int_type_t LookUpAndPlaceAndReturn(int_type_t self_index, vector<SkelGraph>& SG, ReconstructionStructure* RS, int_type_t* voxel_points_map, vector<double>& points,
		int_type_t x_index, int_type_t y_index, int_type_t z_index, int_type_t sub_x, int_type_t sub_y, int_type_t sub_z, int_type_t& point_counter,
		int_type_t big_number){

	// voxel points map  ....
	// relative to self ...
	// x = [0, 1], y =[0, 1], z = [0, 1]
	// then index is x*4 + y*2 + z
	int_type_t nx, ny, nz;
	int_type_t bx, by, bz; // bin x y z
	int_type_t point_index;
	int_type_t point_index_neighbor;
	int_type_t neighbor_index;
	bool contained;

	bx = (sub_x == x_index + 1); // 0 if sub_x = x_index*2 b/c the sub_x%x_index =  0;
	by = (sub_y == y_index + 1);
	bz = (sub_z == z_index + 1);

	point_index = self_index*8 + 4*bx + 2*by + bz;

	if (voxel_points_map[point_index] < big_number){


	}	else {
		voxel_points_map[point_index] = point_counter;

		points.push_back(RS->BB[0][0] + (double(2*sub_x))*RS->inter_voxel_distance/2.0);
		points.push_back(RS->BB[0][1] + (double(2*sub_y))*RS->inter_voxel_distance/2.0);
		points.push_back(RS->BB[0][2] + (double(2*sub_z))*RS->inter_voxel_distance/2.0);

		// need to go through the neighbors and see if they have this as a point.
		for (int j =0, jn = SG[self_index].neighbors.size(); j < jn; j++){

			neighbor_index = SG[self_index].neighbors[j]->voxel_id;
			RS->ReturnXYZIndicesFromIndex(SG[self_index].neighbors[j]->grid_id, nx, ny, nz);

			contained = ((sub_x == nx) || (sub_x == nx + 1)) && ((sub_y == ny) || (sub_y == ny + 1)) && ((sub_z == nz) || (sub_z == nz + 1));


			if (contained){
				//cout << "Contained ! " << contained << endl;
				bx = (sub_x == nx + 1); // 0 if sub_x = x_index*2 b/c the sub_x%x_index =  0;
				by = (sub_y == ny + 1);
				bz = (sub_z == nz + 1);

				point_index_neighbor = neighbor_index*8 + 4*bx + 2*by + bz;
				voxel_points_map[point_index_neighbor] = point_counter;
			}
		}




		point_counter++;
	}

	return voxel_points_map[point_index] ;
}

void ReconstructionStructure::GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration,
		string prefix, vector<int>& cs, bool default_string){

	int newcs[3];
	newcs[0] = cs[0];
	newcs[1] = cs[1];
	newcs[2] = cs[2];

	GenerateAndWriteSurfaceInPlyFormat(outdir, iteration, prefix, newcs, default_string);
}


void ReconstructionStructure::GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration,
		string prefix, int* cs, bool default_string){


	string filename;
	if (default_string){
		filename = outdir + "smoothed_files/" + prefix + ToString<int>(iteration) + ".ply";
	}	else {
		filename = outdir + prefix + ToString<int>(iteration) + ".ply";
	}
	// we'll refer to the points on the mesh by a subsampled version of the
	vector<double> Xhm(3,0);
	vector<int_type_t> sub_number_per_dim(3, 0);

	vector<int_type_t> color(3, 100);
	color[1] = 200;  // was green -- too dark?
	color[2] = 200;
	color[0] = 200;

	if (prefix == "model"){
		color[0] = 255;
		color[1] = 0;
		color[2] = 0;
	}	else {
		color[0] = 0;
		color[1] = 0;
		color[2] = 255;
	}


	if (cs != 0){
		color[0] = cs[0];
		color[1] = cs[1];
		color[2] = cs[2];
	}


	for (int i = 0; i < 3; i++){
		sub_number_per_dim[i] = number_voxels_per_dim[i]*2 + 1;
	}

	vector<SkelGraph> SurfaceGraph;
	CreateSurfaceGraphFromGrid(*this, SurfaceGraph);


	cout << "After SG ... " << SurfaceGraph.size() << endl;


	// now process the points to create unique indices
	int_type_t big_number = sub_number_per_dim[0]* sub_number_per_dim[1]*sub_number_per_dim[2];
	//vector<int_type_t> voxel_points_map;

	int_type_t* voxel_points_map = new int_type_t[SurfaceGraph.size()*8];
	for (int_type_t i = 0, in = SurfaceGraph.size()*8; i < in; i++){
		voxel_points_map[i] = big_number;
	}

	int_type_t point_counter = 0;




	int_type_t above_x, below_x, above_y, below_y, above_z, below_z;
	int_type_t x_index, y_index, z_index;
	bool outer_shell;
	int_type_t voxel_index = 0;


	//	vector<int_type_t> start_index_for_face;
	// faces are all 4 sided
	vector<int_type_t> faces;
	vector<double> points;
	vector<int_type_t> current_face(4, 0);
	outer_shell = false;
	int_type_t current_point_index;
	int_type_t sub_x, sub_y, sub_z;

	//cout << "Line 273 " << endl;cin >> ch;
	for (int_type_t i = 0, in = SurfaceGraph.size(); i < in; i++){
		// Setting up
		voxel_index = SurfaceGraph[i].grid_id;
		outer_shell = false;

		ReturnXYZIndicesFromIndex(voxel_index, x_index, y_index, z_index);

		// TESTING X DIRECTION ..... ////

		if (x_index == 0){
			outer_shell = true;
		}	else {
			below_x = ReturnIndexFromXYZIndices(x_index - 1, y_index, z_index);

			if (configuration_grid[below_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			// below x is a face ....


			// Do I share points on this face with any of the neighbors?
			// come up with a list using the subsample indices.  Look up -- increment and place if so.

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			///////////////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			///////////////////////////




			for (int f = 3; f >= 0; f--){
				faces.push_back(current_face[f]);
			}


		}

		// + x
		outer_shell = false;
		if (x_index == number_voxels_per_dim[0] - 1){
			outer_shell =true;
		}	else {
			above_x = ReturnIndexFromXYZIndices(x_index + 1, y_index, z_index);

			if (configuration_grid[above_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}



		}


		/// -y

		outer_shell = false;
		if (y_index == 0){
			outer_shell =true;
		}	else {
			below_y = ReturnIndexFromXYZIndices(x_index, y_index - 1, z_index);

			if (configuration_grid[below_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////


			for (int f = 3; f >= 0; f--){
				faces.push_back(current_face[f]);
			}

		}

		// + y
		outer_shell = false;
		if (y_index == number_voxels_per_dim[1] - 1){
			outer_shell =true;
		}	else {
			above_y = ReturnIndexFromXYZIndices(x_index, y_index + 1, z_index);

			if (configuration_grid[above_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

		}
		// -z
		outer_shell = false;

		if (z_index == 0){
			outer_shell =true;
		}	else {
			below_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index - 1);

			if (configuration_grid[below_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			//			for (int f = 0; f < 4; f++){
			//				faces.push_back(current_face[f]);
			//			}

			for (int f = 3; f >= 0; f--){
				faces.push_back(current_face[f]);
			}

		}

		// + z
		outer_shell = false;

		if (z_index == number_voxels_per_dim[2] - 1){
			outer_shell =true;
		}	else {
			above_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index + 1);

			if (configuration_grid[above_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}
		}

	}

	cout << "After processing points " << endl;
	cout << "Number points " << point_counter << endl;
	//cin >> ch;

	delete [] voxel_points_map;
	//
	//	//cout << "Line 896 " << endl;cin >> ch;
	//	delete [] dense_to_sparse_map;
	//
	cout << "Before actual write ... " << endl;
	//	char fg; cin >> fg;
	if (points.size() > 0){
		WritePlyFile(filename, sub_number_per_dim, points, faces, color);
	}
	cout << "After actual write" << endl;

	//exit(1);
	//	//	cin >> fg;


}



void WritePlyFileTriangles(string outfile,
		vector<double>& points,
		vector<int_type_t>& faces,
		vector<int>& color){

	// each vertex needs a color ....

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << points.size()/3 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << faces.size()/3 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	out << "end_header" << endl;


	for (int_type_t i = 0; i < (points.size()/3); i++){
		out << points[3*i] << " " <<  points[3*i + 1]  << " " <<  points[3*i + 2] << " " << color[0] << " " <<  color[1] << " " << color[2] << " 255" << endl;
	}

	for (int_type_t i = 0; i < (faces.size()/3); i++){
		out << "3 " << faces[3*i]  << " " <<  faces[3*i + 1]   << " " <<  faces[3*i + 2]  << endl; //" " << color[0] << " " << color[1] << " " << color[2] << " 255 " << endl;
	}

	out << endl;

	out.close();
}
