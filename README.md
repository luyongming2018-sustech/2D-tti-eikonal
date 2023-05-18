# 2D-tti-eikonal
#Computers and Geosciences
#Calculating traveltimes in 2D general tilted transversely isotropic media using fast sweeping method
luym@sustech.edu.cn

#### Description
{**Traveltime calculations for three wave modes*}

1. cd /2D_tti_eikonal_qp or cd /2D_tti_eikonal_qsv or cd /2D_tti_eikonal_qsh;
you need to copy the file para to /2D_tti_eikonal_qsv and /2D_tti_eikonal_qsh.
2. cd /src, make
3. run: ./travel  
4. obtaining the traveltime results in data file, you can use xcontour or matlab to plot the traveltime.

If you change the vel model, put the parameters into para file.
Then cd src file and change the corresponding model_parameter.dat
