#include "ConversionFunctions.hpp"
#include "DirectoryFunctions.hpp"


bool sort_ascending(pair<int, int_type_t> a, pair<int, int_type_t> b){
	return a < b;
}

bool sort_descending(pair<int, int_type_t> a, pair<int, int_type_t> b){
	return a > b;
}


void ReadIndividualImages(string write_directory){

	vector<double> pA(3, 0);
	vector<double> pB(3, 1000);
	ifstream in;
	vector<vector<double> > BB;
	string filename;


	string image_directory = write_directory + "/rawimage";

	vector<string> filenames;
	ReadDirectory(image_directory, filenames);

	ReconstructionStructure R_model;
	cv::Mat image_raw;


	for (int_type_t i = 0; i < filenames.size(); i++){
		filename = write_directory + "/rawimage/" + filenames[i];
		image_raw = cv::imread(filename.c_str(), CV_LOAD_IMAGE_GRAYSCALE);

		if (i == 0){
			// x = number of slices, y = rows, z = cols
			pA[0] = 0;
			pA[1] = 0;
			pA[2] = 0;
			pB[0] = filenames.size();
			pB[1] = image_raw.rows;
			pB[2] = image_raw.cols;


			BB.push_back(pA);
			BB.push_back(pB);

			// intervoxel distance is 1
			R_model.CreateGridMinimal(BB, 1);
		}

		for (int_type_t j = 0; j < R_model.number_voxels_per_dim[1]*R_model.number_voxels_per_dim[2]; j++){
			if (image_raw.data[j] == 0){
				R_model.configuration_grid[i*R_model.number_voxels_per_dim[1]*R_model.number_voxels_per_dim[2] + j] = true;
			}	else {
				R_model.configuration_grid[i*R_model.number_voxels_per_dim[1]*R_model.number_voxels_per_dim[2] + j] = false;
			}

		}


	}


	// for troubleshooting
	R_model.GenerateAndWriteSurfaceInPlyFormat(write_directory, 0, "converted_model", NULL, false);


	// write the model in our format
	std::ofstream out;

	string write_file = write_directory + "0.txt";
	out.open(write_file.c_str());

	for (int_type_t i = 0; i < R_model.number_voxels_grid; i++){
		if (R_model.configuration_grid[i] == false){
			out << i << " ";
		}
	}
	out << "-1" << endl;

	out.close();

	write_file = write_directory + "BB.txt";
	out.open(write_file.c_str());
	out << "1" << endl;
	out << "0 0 0" << endl;
	out << R_model.number_voxels_per_dim[0] << " " << R_model.number_voxels_per_dim[1] << " " << R_model.number_voxels_per_dim[2] << endl;

	out.close();

}

void WriteIndividualImages(string write_directory){

	vector<double> pA(3, 0);
	vector<double> pB(3, 1000);
	ifstream in;
	vector<vector<double> > BB;
	string filename;


	string image_directory = write_directory + "/processedimage";


	ReconstructionStructure R_model;
	string BBfile = write_directory + "BB.txt";
	string object_file = write_directory + "result.txt";
	double division;


	in.open(BBfile.c_str());

	if (!in){
		cout << "You forgot the BB files or the path is wrong" << endl << BBfile << endl;;
		exit(1);
	}

	in >> division;
	in >> pA[0] >> pA[1] >> pA[2];
	in >> pB[0] >> pB[1] >> pB[2];

	BB.push_back(pA);
	BB.push_back(pB);

	mkdir (image_directory.c_str(), S_IRWXU);



	ifstream oin;
	oin.open(object_file.c_str());


	if (oin.fail()){
		oin.close();

		cout << "file " << object_file << " returns error " << endl;
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
	}



	cv::Mat image_raw(R_model.number_voxels_per_dim[1], R_model.number_voxels_per_dim[2], CV_8UC1, cv::Scalar(0));

	for (int_type_t i = 0; i < R_model.number_voxels_per_dim[0]; i++){

		string number = ToString<int_type_t>(i);
		while (number.size() < 4){
			number = "0" + number;
		}

		filename = image_directory + "/"  + number + ".tif";


		for (int_type_t j = 0; j < R_model.number_voxels_per_dim[1]*R_model.number_voxels_per_dim[2]; j++){
			if (R_model.configuration_grid[i*R_model.number_voxels_per_dim[1]*R_model.number_voxels_per_dim[2] + j] == true){
				image_raw.data[j] = 0;
			}	else {
				image_raw.data[j] = 255;
			}
		}

		cv::imwrite(filename.c_str(), image_raw);

	}




}


