# 2D-tti-eikonal
#Computers and Geosciences
#Calculating traveltimes in 2D general tilted transversely isotropic media using fast sweeping method
luym@sustech.edu.cn


2D-TTI-Eikonal: Calculating Traveltimes in 2D General Tilted Transversely Isotropic Media Using the Fast Sweeping Method
Introduction
This module is designed to calculate traveltimes for different wave modes (qP, qSV, and qSH) in 2D general tilted transversely isotropic (TTI) media. 
The fast sweeping method is employed to ensure computational efficiency and accuracy, making this tool suitable for geophysical applications that 
require precise traveltime computations in complex anisotropic media.

Wave Mode Calculation
The module supports three wave modes:

1.qP Wave Mode:
To calculate traveltimes for the qP wave, navigate to the 2D_tti_eikonal_qp directory:
cd /2D_tti_eikonal_qp; 'makdir data' to put the calculation results;
Run the calculation following the steps outlined below.

2.qSV Wave Mode:
For qSV wave traveltimes, navigate to the 2D_tti_eikonal_qsv directory:
cd /2D_tti_eikonal_qsv;'makdir data' to put the calculation results;
Ensure that the parameter file is copied to the corresponding directory:
cp para /2D_tti_eikonal_qsv/;
Continue with the steps for running the calculation.

3.qSH Wave Mode:
To calculate qSH wave traveltimes, navigate to the 2D_tti_eikonal_qsh directory:
cd /2D_tti_eikonal_qsh;'makdir data' to put the calculation results;
Again, copy the parameter file to the directory:
cp para /2D_tti_eikonal_qsh/;
Follow the remaining steps to compute the traveltimes.
Running the Simulation

Compile the Source Code:
First, navigate to the source directory:
cd /src;
Compile the source code using the make command:
make;
This will generate the executable file necessary to run the traveltime calculation.
Execute the Traveltime Calculation:
Run the traveltime computation by executing the following command:
./travel;

The output will be generated in the form of data files that contain the calculated traveltimes.
Visualizing the Results:

To visualize the traveltime results, you can use xcontour or MATLAB. These tools allow for the plotting of the traveltime data, enabling you to analyze the results graphically.
Modifying the Velocity Model
If you need to change the velocity model parameters:

Update the Parameter File:
Modify the para file with the desired parameters for your new velocity model.
Update the Model Parameter File:
Navigate to the src directory:
cd /src;
Modify the model_parameter.dat file to reflect the new model settings.
Recompile and Rerun:

After making the necessary changes, recompile the source code by running 'make' in the src directory, then execute the traveltime calculation as described above.
This workflow allows for flexibility in testing different wave modes and velocity models, making it a versatile tool for traveltime computation in 2D TTI media.
