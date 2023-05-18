#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define abs11(a) ((a)>0?(a):-(a))
#define min(a,b) ((a)<(b)?(a):(b))
#define max_2(a,b) ((a)>(b)?(a):(b))
void Test_causality(double V1, double V2, int num, int *flag)
{
	double v1,v2;
	double tmp,tmp1;

	v1=-V1;v2=-V2;
	tmp=v1*v2;
	tmp1=v2/v1;

	if(tmp>0.0){
		if(tmp1<1.0){
			if(num==6 || num==2){
				*flag=1;
			}else{
				*flag=0;
			}
		}else{
			if(num==7 || num==3){
				*flag=1;
			}else{ 
				*flag=0;
			}
		}
	}else{
		if(tmp1<-1.0){
			if(num==0|| num==4){
				*flag=1;
			}else{
				*flag=0;
			}	
		}else{
			if(num==1 || num==5){
				*flag=1;
			}else{ 
				*flag=0;
			}
		}
	}
		
		
}

//solve er ci equation using analytical method
int Solquar(double a, double b,double c, double root[2])
{

    int flag;
	double dis;
	
	dis=b*b-4*a*c;//printf("dis=%e\n",dis);//sleep(1);
	if(dis<0){
		return 0;
	}else{
		root[0]=(-b+sqrt(dis))/(2*a);//printf("root0=%e\n",root[0]);sleep(1);
		root[1]=(-b-sqrt(dis))/(2*a);
		return 1;
	}
}

void Wave_surfacetest_sh(double  px,double  pz,double  alphy,double  beta,double  gama, int *flag, double *fv)
{
	double  value;
	double  tb,terr; 

	tb=1.0/(beta*beta);
	*fv=pz*pz+(2*gama+1)*px*px-tb;
	terr=abs11(pz*pz+(2*gama+1)*px*px-tb);

	if(terr<1.0e-6){
		*flag=1;
	}else{
		*flag=0;
	}
}
void Wave_surfacetest_qsh(double  px,double  pz,double  alphy,double  beta,double  gama, int *flag)
{
	double  value;
	double  tb,terr; 

	tb=1.0/(beta*beta);
	terr=abs11(pz*pz+(2*gama+1)*px*px-tb);

	if(terr<1.0e-6){
		*flag=1;
	}else{
		*flag=0;
	}
}
void Cal_group_pxz_sh(double  Px,double  Pz,double * V1,double * V2,double  alphy,double  beta,double  gama,int num0)
{ 
    double pz,px;
	double  partialB,partialC;
	double  B,C,sqqr,slowppx,low;
	double  tp,ts;

	px=Px;
	pz=Pz;
	
	tp=alphy*alphy;
	ts=beta*beta;

	slowppx=2*(2*gama+1)*px;
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
}
void Cal_group_pxz_qsh(double  Px[2],double  Pz[2],double * V1,double * V2,double  alphy,double  beta,double  gama,int k, int num0)
{ 
    double pz,px;
	double  partialB,partialC;
	double  B,C,sqqr,slowppx,low;
	double  tp,ts;

	if(k==0){
		px=Px[0];
		pz=Pz[0];
	}else if(k==1){
		 px=Px[1];
         pz=Pz[1];
	}
	
	
	tp=alphy*alphy;
	ts=beta*beta;

	slowppx=2*(2*gama+1)*px;
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
}
//solve Tc using Ta and Tb in the local solver, and solve  quadraticequation.
int Solve_local_traveltime(int cell_num,double Tad[8],double Tc[2],double d,double Px[2],double Pz[2], double alphy, double beta,
                           double gama,double Ragl,double downt,double upt)
{
	double pz,px;
	double cp0,cp1,cp2,cp3,cp4;
	double Ta,Tb;
	double mid1,mid2;
	double root[2];
	double A,B,C,D;
	double tv,ts;
	int flag=10;
	double sa,ca,sa2,ca2,sa3,ca3,sa4,ca4;
	int i;

    switch(cell_num){ 
			case 0:{
				Ta=Tad[7];
				Tb=Tad[0];
				break;
			}
			case 1:{
				Ta=Tad[7];
				Tb=Tad[6];
				break;
			}
			case 2:{
				Ta=Tad[5];
				Tb=Tad[6];
				break;
			}
			case 3:{
				Ta=Tad[5];
				Tb=Tad[4];
				break;
			}
			case 4:{
				Ta=Tad[3];
				Tb=Tad[4];
				break;
			}
			case 5:{
				Ta=Tad[3];
				Tb=Tad[2];
				break;
			}
			case 6:{
				Ta=Tad[1];
				Tb=Tad[2];
				break;
			}
			case 7:{
				Ta=Tad[1];
				Tb=Tad[0];
				break;
			}
 	}

	ts=1.0/beta/beta;

	A=1+2*gama;
	sa=sin(Ragl);
	ca=cos(Ragl);
	sa2=sa*sa;
	ca2=ca*ca;

    if(cell_num==7||cell_num==4){
		px=-(Tb-Ta)/d;
		cp2=ca2+A*sa2;
		cp1=px*(2*A*sa*ca-2*sa*ca);
		cp0=-ts+px*px*(sa2+(2*gama+1)*ca2);//*/
		//printf("cp2,cp1,cp0=%e,%e,%e\n",cp2,cp1,cp0);
		flag=Solquar(cp2,cp1,cp0,root);
		if(flag==0){
			return flag;
		}
		//printf("halou\n");sleep(1);
		for(i=0;i<2;i++){
			Px[i]=px;
			Pz[i]=root[i];
			//printf("pz=%e\n",Pz[i]);sleep(1);
		}
		
		if(cell_num==7){
			for(i=0;i<2;i++){
				Tc[i]=-root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}return flag;
			
		}else if(cell_num==4){
			for(i=0;i<2;i++){
				Tc[i]=root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}return flag;
		}
    }
    if(cell_num==1||cell_num==6){
		pz=(Ta-Tb)/d;

		cp2=sa2+A*ca2;
		cp1=pz*(2*A*sa*ca-2*sa*ca);
		cp0=-ts+pz*pz*(ca2+A*sa2);//*/

		flag=Solquar(cp2,cp1,cp0,root);

	    if(flag==0){return flag;}
		for(i=0;i<2;i++){
				Pz[i]=pz;
				Px[i]=root[i];
		}

	    if(cell_num==1){ 
			for(i=0;i<2;i++){
				Tc[i]=root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}return flag;
		 }else if(cell_num==6){
		    for(i=0;i<2;i++){
				Tc[i]=-root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}
			return flag;
		  }
    }
    if(cell_num==0||cell_num==3){
		px=(Tb-Ta)/d;
		cp2=ca2+A*sa2;
		cp1=px*(2*A*sa*ca-2*sa*ca);
		cp0=-ts+px*px*(sa2+(2*gama+1)*ca2);//*/

		flag=Solquar(cp2,cp1,cp0,root);
	    if(flag==0){return flag;}

		for(i=0;i<2;i++){
				Px[i]=px;
				Pz[i]=root[i];
		}
		if(cell_num==0){
			for(i=0;i<2;i++){
				Tc[i]=-root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}
			return flag;
		}else if(cell_num==3){ 
			for(i=0;i<2;i++){
				Tc[i]=root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}
			return flag;
		}
	}
	if(cell_num==2||cell_num==5){
		pz=(Tb-Ta)/d;
		cp2=sa2+A*ca2;
		cp1=pz*(2*A*sa*ca-2*sa*ca);
		cp0=-ts+pz*pz*(ca2+A*sa2);//*/

		flag=Solquar(cp2,cp1,cp0,root);
		
	    if(flag==0){
	    	return flag;
	    }
		for(i=0;i<2;i++){
				Pz[i]=pz;
				Px[i]=root[i];
		}
	    if(cell_num==2){
			for(i=0;i<2;i++){
				Tc[i]=root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}
			return flag;
		}else if(cell_num==5){
			for(i=0;i<2;i++){
				Tc[i]=-root[i]*d+Tb;
				if(Tc[i]<Tb && Tc[i]<Ta){
					Tc[i]=-1;
				}
			}
			return flag;
		}
    }
}

double update_traveltime(double T,double Tad[8],double Vgr[8],double d)
{     
    double Tnow,Tcom[8];
    double tmp;
    int i;

    tmp=sqrt(2)*d;

    for(i=0;i<7;i=i+2){
        Tcom[i]=Tad[i]+d/Vgr[i];
        Tcom[i+1]=Tad[i+1]+tmp/Vgr[i+1];
    }
    Tnow=1.0e20;

    for(i=0;i<8;i++){ 
        if(Tcom[i]<Tnow){  
            Tnow=Tcom[i];
        }
    }
    return Tnow;
}
//***************************************************************************************
//function for ray tracing
void Cal_group_pxz_tracing(double Px,double Pz, double * V1, double * V2,double alphy,double beta,double delta, double eps)
{ 
    double pz,px;
	double partialB,partialC;
	double B,C,sqqr,slowppx,low;
	double tp,ts;

	px=Px;
	pz=Pz;
	
	tp=alphy*alphy;
	ts=beta*beta;

	partialB=-4*(1+delta+(eps-delta)*tp/ts)*px;
	B=1.0/tp+1.0/ts-2.0*(1.0+delta+(eps-delta)*tp/ts)*px*px;
	C=((1.0+2.0*eps)*px*px-1.0/tp)*(px*px-1.0/ts);
	partialC=2.0*(1.0+2.0*eps)*px*(px*px-1.0/ts)+((1.0+2.0*eps)*px*px-1.0/tp)*2.0*px;
	sqqr=sqrt(B*B-4*C);
	slowppx=-0.5*partialB+1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
}
//***************************************************************************************
void slowness_cal_groupvel_sh(double vp,double  vs,double gama,double  Px,double *vx,double *vz)
{
  double   a,b,c;
  double solu,pxx[2],pzz[2];
  double Vx,Vz;

   a=1;
   b=1.0/vs/vs;
   c=-(1+2*gama)*Px*Px+b;

   solu=c;
   if(solu<0){
   	 pzz[0]=0;
   }
   pxx[0]=Px;
   pzz[0]=sqrt(solu);
   
   Cal_group_pxz_qsh(pxx,pzz,&Vx,&Vz,vp,vs,gama,0,1);

   *vx=Vx;
   *vz=Vz;
}
//approximate formula for cal group vel in fixed direction
void cal_group_velocity_qsh(double Vag[8], double vp, double vs0,double gama,double Ragl)
{
	
	double a1,Cs,b1;
	double interval,Vg;
	int i,j;
	double pi=3.1415926535898;

	interval=pi/4.0;
	for(i=0;i<8;i++){
		Cs=cos(interval*i+Ragl);
		a1=vs0*sqrt(1+2*gama);
		b1=sqrt(1+2*gama*Cs*Cs);
		Vg=a1/b1;
		Vag[i]=Vg;
		//printf("vg[%d]=%f\n",i,Vg);sleep(1);
	}

}
//compute vs phase velocity in each angle[0-90]degree
void vel_s_cal(double *v,double ang, double gama,double vs0)
{
	double sina,cosa,a1,b1;
	
	sina=sin(ang);
	cosa=cos(ang);
	a1=(1.0+2.0*gama)*sina*sina;
	b1=cosa*cosa;
	(*v)=sqrt(a1+b1)*vs0;
}
void vel_sd_cal(double *v_d,double ang,double gama,double vs0)
{
	double sina,cosa,sin2a,cos2a,a1,b1;

	sina=sin(ang);
	cosa=cos(ang);

	a1=(1.0+2.0*gama)*sina*sina+cosa*cosa;
	b1=2.0*gama*sina*cosa*vs0;
	(*v_d)=b1/sqrt(a1);
}
//compute fixed direction group velocity  using iter method
void cal_group_velocity_qsh_iter(double Vag[8],double vp, double vs0,double gama,double tilt1)
{
	
	double ang_s[761];
	double vs_d,cos2a;
	double a1,b1,sina,cosa,tilt;
	double sx,cx,gx,gz;
	double vs;
	int i,j;
	double pi=3.1415926535898;
	double tmp;

	//tilt1*=-1;
	//printf("cal grou vel\n");
	//printf("tilt1=%f\n",tilt1);
	if(tilt1>0){
		for(i=0;i<8;i=i+2){
				if(i==0){
					for(j=270;j<360;j++){
						b1=360.0/180.0*pi-(pi*0.5-tilt1);

						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,gama,vs0);
						vel_sd_cal(&vs_d,a1,gama,vs0);

						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}else if(i==2){
						for(j=0;j<90;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);

							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=90.0/180.0*pi-(pi*0.5-tilt1);
							
							ang_s[j]=gz/gx-tan(b1);//printf("b1,tan=%e,%e\n",b1,tan(b1));
							//printf("a1=%f,j=%d,ans=%f\n",a1*180/pi,j,ang_s[j]);sleep(1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}

				}else if(i==4){
						for(j=90;j<180;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=180.0/180.0*pi-(pi*0.5-tilt1);
	
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
						}
				}else{
						for(j=180;j<270;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=270.0/180.0*pi-(pi*0.5-tilt1);
	
							ang_s[j]=gz/gx-tan(b1);
	
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}
				}
		}
	}else if(tilt1<0){
			for(i=0;i<8;i=i+2){
				if(i==0){
					for(j=0;j<90;j++){
						b1=pi*0.5+tilt1;
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,gama,vs0);
						vel_sd_cal(&vs_d,a1,gama,vs0);
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}else if(i==2){
						for(j=90;j<180;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=180.0/180.0*pi+tilt1;
							
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}

				}else if(i==4){
						for(j=180;j<270;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=270.0/180.0*pi+tilt1;
	
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
						}
				}else{
						for(j=270;j<360;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,gama,vs0);
							vel_sd_cal(&vs_d,a1,gama,vs0);
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=360.0/180.0*pi+tilt1;
	
							ang_s[j]=gz/gx-tan(b1);
	
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}
				}
		}
	
	}else{
		
			//vti case group vel in 0, 90, 180, 270 degree.
			Vag[0]=vs0;
			Vag[2]=vs0*sqrt(1.0+2.0*gama);
			Vag[4]=vs0;
			Vag[6]=vs0*sqrt(1.0+2.0*gama);
	
		
	}
	tilt=tilt1*180.0/pi;
	if(abs11(tilt-45.0)<1.0e-3){
		tilt=tilt1*180.0/pi+1.0;
		tilt1=tilt*pi/180.0;//printf("hello");
	}
	for(i=1;i<8;i=i+2){
		if(i==1){
			if(abs11(tilt)<46.0){
				for(j=0;j<90;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,gama,vs0);
					vel_sd_cal(&vs_d,a1,gama,vs0);
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=45.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					//printf("a1=%f,j=%d,ans=%f\n",a1*180/pi,j,ang_s[j]);sleep(1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
				for(j=0;j<90;j++){
					a1=(j-90)/180.0*pi;
					vel_s_cal(&vs,a1,gama,vs0);
					vel_sd_cal(&vs_d,a1,gama,vs0);
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=45.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			
		}else if(i==3){
			if(abs11(tilt)<46.0){
				for(j=91;j<180;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,gama,vs0);
					vel_sd_cal(&vs_d,a1,gama,vs0);
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=135.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
					for(j=0;j<90;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,gama,vs0);
						vel_sd_cal(&vs_d,a1,gama,vs0);
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=135.0/180.0*pi-tilt1;
	
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
		}else if(i==5){
			if(abs11(tilt)<46.0){
				for(j=181;j<270;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,gama,vs0);
					vel_sd_cal(&vs_d,a1,gama,vs0);
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=225.0/180.0*pi-tilt1;
				
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
				for(j=90;j<180;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,gama,vs0);
					vel_sd_cal(&vs_d,a1,gama,vs0);
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=225.0/180.0*pi-tilt1;
				
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
		}else{
				if(abs11(tilt)<46.0){
					for(j=271;j<360;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,gama,vs0);
						vel_sd_cal(&vs_d,a1,gama,vs0);
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=315.0/180.0*pi-tilt1;
				
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
				if(abs11(tilt)>=46.0){
					for(j=180;j<270;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,gama,vs0);
						vel_sd_cal(&vs_d,a1,gama,vs0);
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=315.0/180.0*pi-tilt1;
				
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz);//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
			}
	}

}
//compute soure travletime using shooting method
double Initialize_source(double alphy,double beta,double gama,int sx,int sz,int ix,int jz,double d,double angle)
{

	double pxmin,pxmax,pmid;
	double vx,vz,an1;
	double pi=3.1415926535898;
	double angsim,ang;
	double Vp;
	double dis,tmp;
	int i;
	double dsr;

	dsr=sqrt((ix-sx)*(ix-sx)+(jz-sz)*(jz-sz));
	an1=cos(angle)*(jz-sx)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz))+sin(angle)*(ix-sz)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz));
	if(an1<0){
		an1=-an1;
	}

	ang=an1; 

	pxmin=0;
	pxmax=1.0/sqrt(1+2*gama)/beta;
 
	while(abs11(pxmax-pxmin)>1.0e-10){
		pmid=(pxmin+pxmax)/2.0;
		slowness_cal_groupvel_sh(alphy,beta,gama,pmid,&vx,&vz);//printf("vx,vz=%e,%e\n",vx,vz);
		angsim=vx/sqrt(vx*vx+vz*vz);

		if(ang>angsim){
			pxmin=pmid;
		}else{
			pxmax=pmid;
		}   
	}

	Vp=sqrt(vx*vx+vz*vz);
	dis=d*dsr;

	return dis/Vp;
}

//FSM sweeping in four directions
 void firstorder_traveltime(int nx,int nz,double dd,int sx,int sz,double **alpha,double **belta,
			    double **gamma,double **theta,double **T,int width,int total)
{

	int num,flag;  
	double d,V1,V2;
	int i,j,k,l;
	double alphy, beta, deta, eps;
	double Pz[2], Px[2];
	double Tc[2],Tad[8],Tad1[8];
	double Vag[8],Vag1[8];
	double Ragl;
	double Bigvalue=50000,pi=3.14159265359;
	double mid1, mid2;
	int inputnum;
	double angle;
	double upt,downt;
	double **Tpre;
	double sum;
	double sina,cosa;
	double **group_x;
	double **group_z;
	double fa;
	double gama;
	FILE *fp;


	Tpre=(double **)calloc(nx,sizeof(double *));
	for(i=0;i<nx;i++){
	    Tpre[i]=(double *)calloc(nz,sizeof(double));
	}
	group_x=(double **)calloc(nx,sizeof(double *));
	for(i=0;i<nx;i++){
	    group_x[i]=(double *)calloc(nz,sizeof(double));
	}
	group_z=(double **)calloc(nx,sizeof(double *));
	for(i=0;i<nx;i++){
	    group_z[i]=(double *)calloc(nz,sizeof(double));
	}
	d=dd;
	for(i=0;i<nx;i++){
		for(j=0;j<nz;j++){    
		    group_x[i][j]=Bigvalue;  
			group_z[i][j]=Bigvalue;
		}
	}
    // initialize the source
	#pragma omp for private(i,j) 
 	for(i=sx-width+1;i<sx+width;i++){
		if(sz-width+1<0){
			for(j=0;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				gama=gamma[i][j];

				T[i][j]=Initialize_source(alphy,beta,gama,sx,sz,i,j,d,angle);
			} 	
		}else{
			for(j=sz-width+1;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				gama=gamma[i][j];

				T[i][j]=Initialize_source(alphy,beta,gama,sx,sz,i,j,d,angle);
			} 	
		}
    }//exit(0);

	inputnum=0;
	start: inputnum++;
	printf("%s%d\n", "interation number begin:",inputnum);
	//width-=140;
	#pragma omp for private(i,j) 
	for(i=nx-3;i>1;i--){
		for(j=nz-3;j>1;j--){    
			if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
				continue;
			}
			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			gama=gamma[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			cal_group_velocity_qsh(Vag,alphy,beta,gama,Ragl);
			//cal_group_velocity_qsh_iter(Vag,alphy,beta,gama,Ragl);
			//for(k=0;k<8;k++){
			//	printf("vag[%d]=%f,va1[%d]=%f\n",k,Vag[k],k,Vag1[k]);
			//}

			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];

			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

			upt=T[i][j];
	
      		for(num=4;num<6;num++){ 
				if(num==4&&T[i+1][j]>=T[i][j]&&T[i+1][j+1]>=T[i][j]){
					continue;
				}else if(num==5&&T[i+1][j+1]>=T[i][j]&&T[i][j+1]>=T[i][j]){
					continue;
				}
				//T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,gama,Ragl,downt,upt);
		     	flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,gama,Ragl,downt,upt);

		     	if(flag==0){continue;}
     	
         		for(k=0;k<2;k++){ 
		            if(Tc[k]<0){continue;}
		     	    mid1=Px[k];
		     	    mid2=Pz[k];
		     	    Px[k]=mid1*cosa+mid2*sina;
		     	    Pz[k]=-mid1*sina+mid2*cosa;
		     	
		     	    Wave_surfacetest_qsh(Px[k],Pz[k],alphy,beta,gama,&flag);
		     	    if(flag==0){continue;}
		     
		            Cal_group_pxz_qsh(Px,Pz,&V1,&V2,alphy,beta,gama,k,num);

					mid1=V1;
					mid2=V2;
					V1=mid1*cosa-mid2*sina;
					V2=mid1*sina+mid2*cosa;
					Test_causality(V1,V2,num,&flag);
		     	    if(flag){
		     	    	T[i][j]=min(Tc[k],T[i][j]);
							group_x[i][j]=V1;  
							group_z[i][j]=V2;
		     	    }
					
     	 		}//*/
				
     	 	}
        }   
  	}
	//exit(0);
	#pragma omp for private(i,j) //there is problem
	for(i=2;i<nx-2;i++){
		for(j=2;j<nz-2;j++){    
		      	if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
		       		continue;
	
		      	}
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				gama=gamma[i][j];
				sina=sin(Ragl);
				cosa=cos(Ragl);
  
				cal_group_velocity_qsh(Vag,alphy,beta,gama,Ragl);
				//cal_group_velocity_qsh_iter(Vag,alphy,beta,gama,Ragl);

				Tad[0]=T[i-1][j];
				Tad[1]=T[i-1][j+1];
				Tad[2]=T[i][j+1];
				Tad[3]=T[i+1][j+1];
				Tad[4]=T[i+1][j];
				Tad[5]=T[i+1][j-1];
				Tad[6]=T[i][j-1];
				Tad[7]=T[i-1][j-1]; 

				T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);
		
				for(num=0;num<2;num++){ 
				 	if(num==0&&T[i-1][j-1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
				 	    continue;
				 	}else if(num==1&&T[i-1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
				 		continue;
		     		}
					//T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,gama,Ragl,downt,upt);
					//printf("Time=%e\n",T[i][j]);//sleep(1);
                	flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,gama,Ragl,downt,upt);

					if(flag==0){continue;} 
					for(k=0;k<2;k++){ 
						if(Tc[k]<0){continue;}
						mid1=Px[k];
						mid2=Pz[k];
		
						Px[k]=mid1*cosa+mid2*sina;
						Pz[k]=-mid1*sina+mid2*cosa;

				
						Wave_surfacetest_qsh(Px[k],Pz[k],alphy,beta,gama,&flag);

						if(flag==0){continue;}

						 Cal_group_pxz_qsh(Px,Pz,&V1,&V2,alphy,beta,gama,k,num);

						mid1=V1;
						mid2=V2;
						V1=mid1*cosa-mid2*sina;
						V2=mid1*sina+mid2*cosa;

						Test_causality(V1,V2,num,&flag);

		     	    	if(flag) {
				 	    	T[i][j]=min(Tc[k],T[i][j]);
							group_x[i][j]=V1;  
							group_z[i][j]=V2;
							//printf("V1,V2=%e,%e\n",V1,V2);
							//sleep(1);
		     	    	}
		     		}//*/
		     }
  		}
  	}//exit(0);
	#pragma omp for private(i,j) //there is problem
	for(i=nx-3;i>1;i--){
	    for(j=2;j<nz-2;j++){  
      		if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
       			continue;
      		}
			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			gama=gamma[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			cal_group_velocity_qsh(Vag,alphy,beta,gama,Ragl);
			//cal_group_velocity_qsh_iter(Vag,alphy,beta,gama,Ragl);

			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];
   
			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

			for(num=2;num<4;num++){ 
		     	if(num==2&&T[i+1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
		     		continue;
		     	}else if(num==3&&T[i+1][j-1]>=T[i][j]&&T[i+1][j]>=T[i][j]){
		     		continue;
		     	}
				
				//T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,gama,Ragl,downt,upt);
				flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,gama,Ragl,downt,upt);

				if(flag==0){continue;}  //flag=0 means we cannot derive a reasonable 
     	
				for(k=0;k<2;k++){ 
					if(Tc[k]<0){
						continue;
					}
					mid1=Px[k];
					mid2=Pz[k];
					Px[k]=mid1*cosa+mid2*sina;
					Pz[k]=-mid1*sina+mid2*cosa;

					Wave_surfacetest_qsh(Px[k],Pz[k],alphy,beta,gama,&flag);

					if(flag==0){
						continue;
					}	
					 Cal_group_pxz_qsh(Px,Pz,&V1,&V2,alphy,beta,gama,k,num);

					mid1=V1;
					mid2=V2;
					V1=mid1*cosa-mid2*sina;
					V2=mid1*sina+mid2*cosa;
					Test_causality(V1,V2,num,&flag);
					if(flag){
						T[i][j]=min(Tc[k],T[i][j]);
						group_x[i][j]=V1;  
						group_z[i][j]=V2;
					}
				}//*/
				
			}
		}
  	}

	#pragma omp for private(i,j) 
 	for(i=2;i<nx-2;i++){
  		for(j=nz-3;j>1;j--){    
	      	if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
	       		continue;
	      	}

			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			gama=gamma[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			cal_group_velocity_qsh(Vag,alphy,beta,gama,Ragl);
			//cal_group_velocity_qsh_iter(Vag,alphy,beta,gama,Ragl);

			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];

			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

			upt=T[i][j];

	        for(num=6;num<8;num++){ 
     	 		if(num==6&&T[i][j+1]>=T[i][j]&&T[i-1][j+1]>=T[i][j]){
     	 			continue;
     	 		}else if(num==7&&T[i-1][j+1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
     	 			continue;
     	 		}
    			//T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,gama,Ragl,downt,upt);
				flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,gama,Ragl,downt,upt);
				if(flag==0){continue;}  
     	
         		for(k=0;k<2;k++){ 
					if(Tc[k]<0){
						continue;
					}
					mid1=Px[k];
					mid2=Pz[k];
					Px[k]=mid1*cosa+mid2*sina;
					Pz[k]=-mid1*sina+mid2*cosa;

					Wave_surfacetest_qsh(Px[k],Pz[k],alphy,beta,gama,&flag);
					if(flag==0){ 
						continue;
					}
					Cal_group_pxz_qsh(Px,Pz,&V1,&V2,alphy,beta,gama,k,num);
				
					mid1=V1;
					mid2=V2;
					V1=mid1*cosa-mid2*sina;
					V2=mid1*sina+mid2*cosa;
					Test_causality(V1,V2,num,&flag);
					if(flag){
						T[i][j]=min(Tc[k],T[i][j]);
							group_x[i][j]=V1;  
							group_z[i][j]=V2;
					}
				}//*/
				
     	 	}	
        }	
  	}

  	sum=0;

  	if(inputnum==0){
		for(j=2;j<nx-2;j++){
		  for(i=2;i<nz-2;i++){
		     Tpre[j][i]=T[j][i];
		  }
		}
		goto start;
		printf("first order iteration %d: %f  %f\n",inputnum,sum,sum/(nx-4)/(nz-4));	
  	}else{
		sum=0;
	      for(j=2;j<nx-2;j++){
			  for(i=2;i<nz-2;i++){
					sum=sum+fabs(Tpre[j][i]-T[j][i]);
			  }
	      }
	      for(j=2;j<nx-2;j++){
	      	   for(i=2;i<nz-2;i++){
		     		Tpre[j][i]=T[j][i];
		  	   }
	      }
		printf("first order iteration %d: %e  %e\n",inputnum,sum,sum/(nx-4)/(nz-4));
		//sum/=(nx-4)*(nz-4);
		//if(sum>1.0e-30){
		//	goto start;
		//}
  	}
	fp=fopen("../data/group_vel_x.bin","wb");
	for(i=2;i<nx-2;i++){
	    for(j=2;j<nz-2;j++){  
			fa=(double)group_x[i][j];
	        fwrite(&fa,sizeof(double),1,fp);
	    }
	} 
	fclose(fp);//*/
	fp=fopen("../data/group_vel_z.bin","wb");
	for(i=2;i<nx-2;i++){
	    for(j=2;j<nz-2;j++){  
			fa=(double)group_z[i][j];
	        fwrite(&fa,sizeof(double),1,fp);
	    }
	} 
	fclose(fp);//*/	
	total=1;
	if(inputnum<total){
		goto start;
	}
	  
}

