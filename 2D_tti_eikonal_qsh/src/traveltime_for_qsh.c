#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<time.h>
#include "subrounte_qsh.c"
#include <omp.h>
char *sdoc[] = {
"2D TravelTime Computation for TI media",
NULL};
int main()
{
	int  i,j,k,ix,iz;
	double Bvle=50000,pi=3.1415926535898;
	int nx,nz,sx,sz;
	double**vp,**vs,**gamma,**theta;
	float a;
	double Ragl,angle;

	double**T, **Tsec,**Ttr;
	double dd;
	int width;
	int total;
	int flag;
	int nx0,nz0;
	char filename[200],fnum[40],parfile[300];
	FILE *fp,fp1;
	float fa;
	double start,finish;


	//read the model parameter
	// nx,nz,shot_x,shot_z,dx,dz
	sprintf(parfile,"%s","model_parameter.dat");
	if((fp=fopen(parfile,"rb"))==NULL){
		printf("cannot open file:%s\n",parfile);
		exit(0);
	}
	fscanf(fp," model size nx and nz ");
	fscanf(fp, "%d%d", &nx0, &nz0);
	printf("nx,nz=%d,%d\n",nx0,nz0);
	fscanf(fp," shot location sx and sz ");
	fscanf(fp, "%d%d\n",&sx,&sz);
	fscanf(fp," grid interval dx dz");
	fscanf(fp, "%f", &a);
	fclose(fp);
	dd=a;//dd=10.0;
	//nx0=801;nz0=801;sx=400;sz=400;
	nx=nx0+4;
	nz=nz0+4;
	sz+=2;
	sx+=2;
	
	//sx=nx/2;sz=40+2;
	//dd=5;
	width=20;
	total=0;

	printf("2D Traveltime calculation in TI media\n");
	//printf("nx,nz=%d,%d\n",nx,nz);
	if(sx<0||sx>nx||sz<0||sz>nz){
		printf("The location of the shot is out of boundary!!!!!\n");
		exit(0);
	}

	vp=(double**)calloc(nx,sizeof(double*));
	vs=(double**)calloc(nx,sizeof(double *));
	gamma=(double**)calloc(nx,sizeof(double *));
	theta=(double **)calloc(nx,sizeof(double *));
	T=(double **)calloc(nx,sizeof(double *));
	Ttr=(double **)calloc(nx,sizeof(double *));
	Tsec=(double **)calloc(nx,sizeof(double *));

	for(i=0;i<nx;i++){
		vp[i]=(double *)calloc(nz,sizeof(double ));
		vs[i]=(double *)calloc(nz,sizeof(double ));
		gamma[i]=(double *)calloc(nz,sizeof(double ));
		theta[i]=(double *)calloc(nz,sizeof(double ));
		Ttr[i]=(double *)calloc(nz,sizeof(double ));
		Tsec[i]=(double *)calloc(nz,sizeof(double ));
	}
	for(i=0;i<nx;i++){
		T[i]=(double *)calloc(nz,sizeof(double ));
	}
 
	//read Vp0
	sprintf(filename,"../para/vpm1.bin");
	if((fp=fopen(filename,"rb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix++){
		for(iz=2;iz<nz-2;iz++){
			fread(&a,sizeof(float),1,fp);
			vp[ix][iz]=(double)a;
		}
	}
	fclose(fp); 

	//read Vs0
	sprintf(filename,"../para/vs.bin");
	if((fp=fopen(filename,"rb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix++){
		for(iz=2;iz<nz-2;iz++){
			fread(&a,sizeof(float),1,fp);
			//a=2000.0;//
			a=vp[ix][iz]/2;
			vs[ix][iz]=(double)a;
		}
	}
	fclose(fp); 

	/*sprintf(filename,"../para/vs_cx.bin");
	if((fp=fopen(filename,"wb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix=ix+2){
		for(iz=2;iz<nz-2;iz=iz+2){	
			a=(float)vs[ix][iz];
			fwrite(&a,sizeof(float),1,fp);
			//a=2000.0;//vp[ix][iz]/1.5;
		}
	}
	fclose(fp);//*/ 
	
	sprintf(filename,"../para/epsilon.bin");
	if((fp=fopen(filename,"rb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix++){
		for(iz=2;iz<nz-2;iz++){
			fread(&a,sizeof(float),1,fp);
			gamma[ix][iz]=(double)a;
		}
	}
	fclose(fp); 
	/*sprintf(filename,"../para/gamma_cx.bin");
	if((fp=fopen(filename,"wb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix=ix+2){
		for(iz=2;iz<nz-2;iz=iz+2){	
			a=(float)gamma[ix][iz];
			fwrite(&a,sizeof(float),1,fp);
			//a=2000.0;//vp[ix][iz]/1.5;
		}
	}
	fclose(fp);//*/ 

	sprintf(filename,"../para/theta.bin");
	if((fp=fopen(filename,"rb"))==NULL){
		printf("cannnot open %s\n",filename);
		exit(0);
	}
	for(ix=2;ix<nx-2;ix++){
		for(iz=2;iz<nz-2;iz++){
			fread(&a,sizeof(float),1,fp);
			//theta[ix][iz]=a*pi/180.0;
				a=30.0;
		        theta[ix][iz]=a*pi/180.0;
		}
	}
	fclose(fp); 
	printf("begin cal traverl time\n");
	for(ix=0;ix<nx;ix++){
		for(iz=0;iz<nz;iz++){
			T[ix][iz]=Bvle;
			Ttr[ix][iz]=Bvle;
			Tsec[ix][iz]=Bvle;
		}
	}	
	
	start=clock();
	firstorder_traveltime(nx,nz,dd,sx,sz,vp,vs,gamma,theta,T,width,total);

	/*for(j=2;j<nz-2;j++){
	    for(i=2;i<nx-2;i++){
	        Tsec[i][j]=T[i][j]; 
	    }
	}//*/
	printf("travetime=%f\n",T[sx][nz/2]);
	printf("travetime=%f\n",T[sx][nz-3]);

	//secondorder_traveltime(nx,nz,dd,sx,sz,vp,vs,delta,epsilon,theta,Tsec,width,total);

	finish=clock();
	printf( "Total calculation time is %f seconds\n",(finish - start) / CLOCKS_PER_SEC);
	
	fp=fopen("../data/qsh_traveltime.bin","wb");
	for(i=2;i<nx-2;i++){
	    for(j=2;j<nz-2;j++){  
		fa=(float)T[i][j];
	        fwrite(&fa,sizeof(float),1,fp);
	    }
	} 
	fclose(fp);
	/*fp=fopen("../data/second_traveltime.bin","wb");
	for(i=2;i<nx-2;i++){
	    for(j=2;j<nz-2;j++){  
	        fa=(float)Tsec[i][j];
			fwrite(&fa,sizeof(float),1,fp);
	    }
	} 
	fclose(fp);//*/

	printf("calation is over\n");

	for (i=0;i<nx;j=i++){
		free(vp[i]);
		free(vs[i]);
		free(T[i]);
		free(gamma[i]);
		free(theta[i]);
		free(Ttr[i]);
		free(Tsec[i]);
	}
	free(vp);
	free(vs);
	free(T);
	free(gamma);
	free(theta);
	free(Ttr);
	free(Tsec);

	return 0;
}
