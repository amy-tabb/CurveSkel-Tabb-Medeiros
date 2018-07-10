# CurveSkel-Tabb-Medeiros

BUILDING AND RUNNING

This README file is to accompany code for curve skeletonization of elongated objects that may have noisy surfaces, produced by Amy Tabb and Henry Medeiros as a companion to their paper:
	Fast and robust curve skeletonization for real-world elongated objects

@INPROCEEDINGS{Tabb18Fast,
author={Amy Tabb and Henry Medeiros},
booktitle={2018 IEEE Winter Conference on Applications of Computer Vision (WACV)},
title={Fast and robust curve skeletonization for real-world elongated objects},
year={2018},
pages={1935-1943},
doi={10.1109/WACV.2018.00214},
month={March},}

This paper is also available from arXiv:1702.07619 [cs.CV] -- including the supplementary material.

The code may be used without restriction. If the results of the code are used as a part of a system described in a publication, we request that the authors cite a published version paper.  The citation for the code itself is: 

@electronic{tabb2018skel_code,
author = {Tabb, Amy},
year = {2018},
title = {Code from: Fast and robust curve skeletonization for real-world elongated objects},
doi = {10.15482/USDA.ADC/1399689},
owner = {Ag Data Commons},
howpublished= {\url{http://dx.doi.org/10.15482/USDA.ADC/1399689}}
} 

(Other citation styles are fine -- this is one that worked for us, but we're not totally happy with it.)

Minimal working examples are available here (third resource): http://dx.doi.org/10.15482/USDA.ADC/1399689

The model may be only used for evaluation and debugging purposes of the code, and not used in any other published, or 3D printed, work. However, no gaurantees are expressed or implied of the code or the model.

Mar 19, 2018
Comments/Bugs/Problems: amy.tabb@ars.usda.gov

This README covers instructions for building the code.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Instructions for compilation and linking:
1. This code has been written and tested on Ubuntu 14.04, using Eclipse CDT as the IDE, and is written in C/C++.


2. This code is dependent on the OpenCV-3.* libraries.  These libraries should be in the include path, or specified in your IDE.


3. Compiler flags: we use OpenMP for parallelization and the C++11 standard.  Note that Eclipse's indexer marks some items from C++11 as errors (but still compiles).  
The flags needed using the gnu compiler, openmp, and the C++11 standard are:	
		 -fopenmp  -std=gnu++11
	though depending on the compiler used, you may need different flags: https://www.dartmouth.edu/~rc/classes/intro_openmp/compile_run.html
	

4. 	libraries needed for linking are:
	gomp   [OpenMP]
	opencv_core [OpenCV]
	opencv_legacy
	opencv_highgui
	opencv_imgproc
	opencv_imgcodecs
	
 5. [July 10 update] To build, a command line example using OpenCV 4.0: 
 g++ *.cpp -o curve_skel -fopenmp -std=gnu++11 -Wall -I\usr\local\include -lgomp -lopencv_core  -lopencv_imgproc -lopencv_imgcodecs
 The executable curve_skel is created.
 In the Eclipse CDT IDE, under project properties->C/C++build->settings->Cross G++ compiler->dialect
 	select ISO C++ 11 or greater
	project properties->C/C++build->settings->Cross G++ compiler->miscellaneous
	add the openmp flag -fopenmp [platform dependent]
	project properties->C/C++build->settings->Cross G++ linker->Libraries
	add the libraries specified above or a subset of them.  The OpenCV library is undergoing some changes and the full set 		may not be necessary.
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Instructions for use:
We use a file format for representing models that is composed of a text file containing the model points "0.txt" and bounding box in file "BB.txt".  
We include samples of these two types of files in the "examples" folder.  We also have some conversion functions (see 3 in this section).

1. To run the code using the default threshold (1e-12) for rejecting spurious segments, you call with one argument, the name of the directory containing "0.txt" and "BB.txt":
./program_name directory

2. To run the code with different values of the threshold, you call with two argmuents, the directory name followed by the threshold value [0, 1]:
./program_name directory threshold

3. We have chosen the image sequence representation for conversion.  To convert models from any other format, one can use ImageJ/Fiji to load and then save as an image
sequence.  Once a sequence is created, save in a folder titled "rawimages" within the directory.  Then call with three arguments:
./program_name directory threshold 1

We have included in the "examples" folder an image sequence.  The result is saved in a folder called "processedimages", which one can then load with ImageJ/Fiji to visualize in
one's preferred environment.  The result is also saved in the format described below.



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Results format:
The results are written as .ply meshes.  One can view them with free viewer Meshlab (http://www.meshlab.net/).  A sample of results is shown in the examples, "A_Result" directory.


1. The paths_TabbMedeiros.ply file shows the paths computed by the algorithm.


2. The skel_TabbMedeiros.ply shows the 1D skeleton voxels.

3. The initial.ply file shows a mesh of the original object.  This may be important for verifying that your model files are set up correctly.




