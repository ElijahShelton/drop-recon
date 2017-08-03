# drop-recon
A computational tool for stitching together and analyzing 3D microscope images of droplets. 

Drop-Recon is an Open Source software package developed in MATLAB for stitching together and analyzing 3D microscope images of droplets. The software is being developed in the Campàs Lab at the Unversity of California Santa Barbara. The software is designed and developed by Elijah Shelton and Friedhelm Serwane with input from Otger Campàs. E.S. received support from the NSF GRFP and F.S. received support from Alexander-von-Humboldt Foundation while working on this project. The Drop-Recon source code is distributed under an Open Source license for non-commercial use (see License.md). The file DropRecon3Dv1.m contains source code owned by others and covered under other licences where indicated. Supplied c++ files are covered by the BSD license.

QUICK START:

1. Open MATLAB and navigate to directory containing DropRecon3Dv1.m
2. In the command line, type "DropRecon3Dv1.run" and press enter.
3. Read and accept the licence before continuing.
4. Select a directory containing tifs you would like to analyze. As a demo, select the directory './tifs for test/'.
5. After clicking through the subsequent dialog boxes, the tifs will be analyzed and files containing the results will be saved in subdirectories of the directory containing the tifs. Open the .mat file in matlab to access the measurements of surface coordinates and measured curvatures.

NOTES:
The current version of this software requires additional matlab toolboxes to run.
