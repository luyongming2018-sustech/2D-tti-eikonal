#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define abs11(a) ((a)>0?(a):-(a))
#define min(a,b) ((a)<(b)?(a):(b))
#define max_2(a,b) ((a)>(b)?(a):(b))

//check causality
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
void Cal_group_pxz_new(double Px,double Pz, double * V1, double * V2,double alphy,double beta,double delta, double eps,int num0)
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
//qsv slowness surface
void Wave_surfacetest_vs(double  px,double  pz,double  alphy,double  beta,double  delta,double  epsilon, int *flag, double *fv)
{
	double  a,b,c;
	double  value;
	double  tb,terr; 

	a=1.0/(alphy*alphy);
	tb=1.0/(beta*beta);
	b=-2*(1+delta+(epsilon-delta)*alphy*alphy*tb)*px*px+a+tb;
	c=((1+2*epsilon)*px*px-a)*(px*px-tb);
	if(sqrt(b*b-4*c)>0){
		value=b/2+sqrt(b*b-4*c)/2;
	}else{
		value=0.0;
	}
	*fv=pz*pz-value;
	terr=abs11(pz*pz-value);

	if(terr<1.0e-6){
		*flag=1;
	}else{
		*flag=0;
	}
}
void Cal_group_pxz_new_vs(double  Px,double  Pz,double * V1,double * V2,double  alphy,double  beta,double  delta,double eps,int num0)
{ 
    double pz,px;
	double  partialB,partialC;
	double  B,C,sqqr,slowppx,low;
	double  tp,ts;

	px=Px;
	pz=Pz;
	
	tp=alphy*alphy;
	ts=beta*beta;

	partialB=-4*(1+delta+(eps-delta)*tp/ts)*px;
	B=1.0/tp+1.0/ts-2.0*(1.0+delta+(eps-delta)*tp/ts)*px*px;
	C=((1.0+2.0*eps)*px*px-1.0/tp)*(px*px-1.0/ts);
	partialC=2.0*(1.0+2.0*eps)*px*(px*px-1.0/ts)+((1.0+2.0*eps)*px*px-1.0/tp)*2.0*px;
	if(sqrt(B*B-4.0*C)>0.0){
		sqqr=sqrt(B*B-4.0*C);
	}else{
		sqqr=0;
	}
	slowppx=-0.5*partialB-1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
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
void Cal_group_pxz_vs(double  Px[2],double  Pz[2],double * V1,double * V2,double  alphy,double  beta,double  delta,double eps,int k,int num0)
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

	partialB=-4*(1+delta+(eps-delta)*tp/ts)*px;
	B=1.0/tp+1.0/ts-2.0*(1.0+delta+(eps-delta)*tp/ts)*px*px;
	C=((1.0+2.0*eps)*px*px-1.0/tp)*(px*px-1.0/ts);
	partialC=2.0*(1.0+2.0*eps)*px*(px*px-1.0/ts)+((1.0+2.0*eps)*px*px-1.0/tp)*2.0*px;
	sqqr=sqrt(B*B-4*C);
	slowppx=-0.5*partialB-1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
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
void slowness_cal_groupvel_vs(double   vp,double  vs,double  epsilon,double delta,double  Px,double *vx,double *vz)
{
   double   a,b,c;
  double solu,pxx[2],pzz[2];
  double Vx,Vz;

   a=1;
   b=-2.0*(1.0+delta+(epsilon-delta)*vp*vp/vs/vs)*Px*Px+1.0/vp/vp+1.0/vs/vs;
   c=1.0/vp/vp/vs/vs-(1.0/vp/vp+(1.0+2.0*epsilon)/vs/vs)*Px*Px+(1.0+2.0*epsilon)*Px*Px*Px*Px;

   solu=b/2.0+sqrt(b*b-4*a*c)/2.0;
   pxx[0]=Px;
   pzz[0]=sqrt(solu);
   if(solu<0){
   	 pzz[0]=0;
   }

   Cal_group_pxz_vs(pxx,pzz,&Vx,&Vz,vp,vs,delta,epsilon,0,1);

   *vx=Vx;
   *vz=Vz;
}
//source init using shooting method
double Initialize_source(double alphy,double beta,double delta,double epsilon,int sx,int sz,int ix,int jz,double d,double angle)
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
    an1=cos(angle)*(ix-sx)/dsr+sin(angle)*(jz-sz)/dsr;
	//an1=cos(angle)*(jz-sx)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz))+sin(angle)*(ix-sz)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz));
	if(an1<0){
		an1=-an1;
	}

	ang=an1; 

	pxmin=0;
	//pxmax=1.0/sqrt(1+2*epsilon)/beta;
	pxmax=1.0/beta;
	while(abs11(pxmax-pxmin)>1.0e-8){
		pmid=(pxmin+pxmax)/2.0;
		slowness_cal_groupvel_vs(alphy,beta,epsilon,delta,pmid,&vx,&vz);//printf("vx,vz=%e,%e\n",vx,vz);
		angsim=vx/sqrt(vx*vx+vz*vz);

		if(ang>angsim){
			pxmin=pmid;
		}else{
			pxmax=pmid;
		}   
	}

	Vp=sqrt(vx*vx+vz*vz);
	dis=d*dsr;//printf("vp,dis.vp=%e,%e\n",Vp,dis/Vp);sleep(1);

	return dis/Vp;
}
//compute vs phase velocity in each angle[0-90]degree
void vel_s_cal(double *v,double ang,double delta,double epsi,double ratio)
{
	double f,sina,cosa,a1,b1,cos2a;

	f=1-1.0/ratio/ratio;
	
	sina=sin(ang);
	cosa=cos(ang);
	cos2a=cosa*cosa-sina*sina;

	a1=(1+epsi*sina*sina-f*0.5);
	b1=0.5*f*sqrt(1+4*sina*sina/f*(2*delta*cosa*cosa-epsi*cos2a)+4*epsi*epsi*sina*sina*sina*sina/f/f);
	(*v)=a1-b1;
}
void vel_sd_cal(double *v_d,double ang,double delta,double epsi,double ratio)
{
	double f,sina,cosa,sin2a,cos2a,a1,b1;

	f=1-1.0/ratio/ratio;
	
	sina=sin(ang);
	cosa=cos(ang);
	cos2a=cosa*cosa-sina*sina;
	sin2a=2*sina*cosa;

	a1=(1+2*epsi*sina*sina/f);
	a1=a1*a1-2*(epsi-delta)*sin2a*sin2a/f;

	b1=2*(1+2*epsi*sina*sina/f)*2*epsi/f*2*sina*cosa-2*(epsi-delta)/f*2*sin2a*cos2a*2;
	(*v_d)=2*epsi*sina*cosa-f/4*b1/sqrt(a1);
}
//compute fixed direction group velocity using approximation formula
void cal_group_velocity_vs(double  Vag[8],double  vp,double  vs,double  epsilon,double delta,double  Ragl)
{ 
	int i;
	double  pi=3.1415926535898,interval;
	double  a0,a1,a2,Vg;
	double  len,Sn,s2,s4;

	interval=pi/4.0;
	for(i=0;i<8;i++){
		Sn=sin(interval*i+Ragl);
		s2=Sn*Sn;
		s4=s2*s2;
		a0=1.0/vs/vs;	
		a1=2.0*a0*vp*vp/vs/vs*(delta-epsilon);
		a2=-2.0*a0*vp*vp/vs/vs*(delta-epsilon);
		
		len=a0+a1*s2+a2*s4;//printf("len=%e\n",len);
		Vg=1.0/sqrt(len);
		Vag[i]=Vg;
	}
}
void cal_group_velocity_qsv(double Vag[8], double vp, double vs0,double epsi,double delta,double tilt1)
{
	
	double ang_s[761];
	double f,ratio,vs_d,cos2a;
	double a1,b1,sina,cosa,tilt;
	double sx,cx,gx,gz;
	double vs;
	int i,j;
	double pi=3.1415926535898;
	double tmp;

	ratio=vp/vs0;
	//tilt1*=-1;
	//printf("cal grou vel\n");
	//printf("tilt1=%f\n",tilt1);
	if(tilt1>0){
		for(i=0;i<8;i=i+2){
				if(i==0){
					for(j=270;j<360;j++){
						b1=360.0/180.0*pi-(pi*0.5-tilt1);

						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,delta,epsi,ratio);
						vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
						vs=sqrt(vs);

						vs_d=vs_d/2/vs;
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}else if(i==2){
						for(j=0;j<90;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=90.0/180.0*pi-(pi*0.5-tilt1);
							
							ang_s[j]=gz/gx-tan(b1);//printf("b1,tan=%e,%e\n",b1,tan(b1));
							//printf("a1=%f,j=%d,ans=%f\n",a1*180/pi,j,ang_s[j]);sleep(1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}

				}else if(i==4){
						for(j=90;j<180;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=180.0/180.0*pi-(pi*0.5-tilt1);
	
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
								break;
							}
						}
				}else{
						for(j=180;j<270;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=270.0/180.0*pi-(pi*0.5-tilt1);
	
							ang_s[j]=gz/gx-tan(b1);
	
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
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
						vel_s_cal(&vs,a1,delta,epsi,ratio);
						vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
						vs=sqrt(vs);

						vs_d=vs_d/2/vs;
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}else if(i==2){
						for(j=90;j<180;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=180.0/180.0*pi+tilt1;
							
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}

				}else if(i==4){
						for(j=180;j<270;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=270.0/180.0*pi+tilt1;
	
							ang_s[j]=gz/gx-tan(b1);
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
								break;
							}
						}
				}else{
						for(j=270;j<360;j++){
							a1=(j)/180.0*pi;
							vel_s_cal(&vs,a1,delta,epsi,ratio);
							vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
							vs=sqrt(vs);

							vs_d=vs_d/2/vs;
							gx=sin(a1)*vs+cos(a1)*vs_d;
							gz=cos(a1)*vs-sin(a1)*vs_d;
							b1=360.0/180.0*pi+tilt1;
	
							ang_s[j]=gz/gx-tan(b1);
	
							if(ang_s[j-1]*ang_s[j]<0){
								Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
								break;
							}
					}
				}
		}
	
	}else{
		
			//vti case group vel in 0, 90, 180, 270 degree.
			Vag[0]=vs0;
			Vag[4]=vs0;
			Vag[2]=vs0;
			Vag[6]=vs0;
	
		
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
					vel_s_cal(&vs,a1,delta,epsi,ratio);
					vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
					vs=sqrt(vs);

					vs_d=vs_d/2/vs;
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=45.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					//printf("a1=%f,j=%d,ans=%f\n",a1*180/pi,j,ang_s[j]);sleep(1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
				for(j=0;j<90;j++){
					a1=(j-90)/180.0*pi;
					vel_s_cal(&vs,a1,delta,epsi,ratio);
					vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
					vs=sqrt(vs);

					vs_d=vs_d/2/vs;
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=45.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			
		}else if(i==3){
			if(abs11(tilt)<46.0){
				for(j=91;j<180;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,delta,epsi,ratio);
					vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
					vs=sqrt(vs);

					vs_d=vs_d/2/vs;
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=135.0/180.0*pi-tilt1;
	
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
					for(j=0;j<90;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,delta,epsi,ratio);
						vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
						vs=sqrt(vs);

						vs_d=vs_d/2/vs;
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=135.0/180.0*pi-tilt1;
	
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
		}else if(i==5){
			if(abs11(tilt)<46.0){
				for(j=181;j<270;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,delta,epsi,ratio);
					vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
					vs=sqrt(vs);

					vs_d=vs_d/2/vs;
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=225.0/180.0*pi-tilt1;
				
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
			if(abs11(tilt)>=46.0){
				for(j=90;j<180;j++){
					a1=(j)/180.0*pi;
					vel_s_cal(&vs,a1,delta,epsi,ratio);
					vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
					vs=sqrt(vs);

					vs_d=vs_d/2/vs;
					gx=sin(a1)*vs+cos(a1)*vs_d;
					gz=cos(a1)*vs-sin(a1)*vs_d;
					b1=225.0/180.0*pi-tilt1;
				
					ang_s[j]=gz/gx-tan(b1);
					if(ang_s[j-1]*ang_s[j]<0){
						Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
						break;
					}
				}
			}
		}else{
				if(abs11(tilt)<46.0){
					for(j=271;j<360;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,delta,epsi,ratio);
						vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
						vs=sqrt(vs);

						vs_d=vs_d/2/vs;
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=315.0/180.0*pi-tilt1;
				
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
				if(abs11(tilt)>=46.0){
					for(j=180;j<270;j++){
						a1=(j)/180.0*pi;
						vel_s_cal(&vs,a1,delta,epsi,ratio);
						vel_sd_cal(&vs_d,a1,delta,epsi,ratio);
						vs=sqrt(vs);

						vs_d=vs_d/2/vs;
						gx=sin(a1)*vs+cos(a1)*vs_d;
						gz=cos(a1)*vs-sin(a1)*vs_d;
						b1=315.0/180.0*pi-tilt1;
				
						ang_s[j]=gz/gx-tan(b1);
						if(ang_s[j-1]*ang_s[j]<0){
							Vag[i]=sqrt(gx*gx+gz*gz)*vp;//printf("vag_1=%f\n",Vag[i]);
							break;
						}
					}
				}
			}
	}

}




//-------------------------------------------------------------------------------------------------
double newton(double  px,double pbeg, double alphy,double  beta,double  delta,double  eps,double tmp)
{
	double x,x0;
	double f,f1;
	double  a,b,c;
	double  value;
	double  tb,terr; 
	double  partialB,partialC;
	double  tp,ts;
	double y1,y;
	double thet=60.0;
	double lamda=1.0;
	double pxx,pzz;
	double dpx,dpz,bpz,cpz;
	double pz;

	x=pbeg;
	//pxx=px;
	int nn;
	nn=0;//printf("netown numer=%d\n",nn);	
	do{
		x0=x;
		pz=x0;
		pxx=px*cos(tmp)+pz*sin(tmp);
		pzz=-px*sin(tmp)+pz*cos(tmp);
			
		a=1.0/(alphy*alphy);
		tb=1.0/(beta*beta);
		b=a+tb-2*(1+delta+(eps-delta)*alphy*alphy*tb)*pxx*pxx;
		c=((1+2*eps)*pxx*pxx-a)*(pxx*pxx-tb);

		
		dpx=2.0*px*sin(tmp)*cos(tmp)+2.0*pz*sin(tmp)*sin(tmp);
		dpz=-2.0*px*sin(tmp)*cos(tmp)+2.0*pz*cos(tmp)*cos(tmp);

		bpz=-2.0*(1.0+delta+(eps-delta)*alphy*alphy/beta/beta)*dpx;
		cpz=(1.0+2.0*eps)*dpx*(pxx*pxx-1.0/beta/beta)+((1.0+2.0*eps)*pxx*pxx-1.0/alphy/alphy)*dpx;
		//pd=2.0*dpz-bpz+(b*bpz-2.0*cpz)/sqrt(b*b-4.0*c);
		value=b+sqrt(b*b-4.0*c);

		f=2.0*pzz*pzz-value;//printf("f,px,pz,x0,x=%e,%e,%e,%e,%e\n",f,px,pz,x0,x);sleep(1);
		//printf("f,x0,x=%e,%e,%e\n",f,x0,x);sleep(1);
		if(abs11(f)<1.0e-20){
			break;
		}else{	
			f1=2.0*dpz-bpz-(b*bpz-2.0*cpz)/sqrt(b*b-4.0*c);
			if(abs11(f1)<1.0e-20){
				break;
			}	
			x=x0-f/f1;
		}
		nn++;
	}while(abs11(x-x0)>1e-11);//printf("f1=%e\n",f1);
	//if(nn>4){
		
	//}
	return x;	
}
double newton_px(double  pz,double pbeg, double alphy,double  beta,double  delta,double  eps,double tmp)
{
	double x,x0;
	double f,f1;
	double  a,b,c;
	double  value;
	double  tb,terr; 
	double  partialB,partialC;
	double  tp,ts;
	double y1,y;
	double thet=60.0;
	double lamda=1.0;
	double pxx,pzz;
	double dpx,dpz,bpx,cpx;
	double px;

	x=pbeg;
	//pxx=px;
	int nn;
	nn=0;//printf("netown numer=%d\n",nn);	
	do{
		x0=x;
		px=x0;
		pxx=px*cos(tmp)+pz*sin(tmp);
		pzz=-px*sin(tmp)+pz*cos(tmp);
			
		a=1.0/(alphy*alphy);
		tb=1.0/(beta*beta);
		b=a+tb-2*(1+delta+(eps-delta)*alphy*alphy*tb)*pxx*pxx;
		c=((1+2*eps)*pxx*pxx-a)*(pxx*pxx-tb);

		
		dpz=-2.0*pz*sin(tmp)*cos(tmp)+2.0*px*sin(tmp)*sin(tmp);
		dpx=2.0*pz*sin(tmp)*cos(tmp)+2.0*px*cos(tmp)*cos(tmp);

		bpx=-2.0*(1.0+delta+(eps-delta)*alphy*alphy/beta/beta)*dpx;
		cpx=(1.0+2.0*eps)*dpx*(pxx*pxx-1.0/beta/beta)+((1.0+2.0*eps)*pxx*pxx-1.0/alphy/alphy)*dpx;
		//pd=2.0*dpz-bpz+(b*bpz-2.0*cpz)/sqrt(b*b-4.0*c);
		value=b+sqrt(b*b-4.0*c);

		f=2.0*pzz*pzz-value;//printf("f,px,pz,x0,x=%e,%e,%e,%e,%e\n",f,px,pz,x0,x);sleep(1);
		//printf("f,x0,x=%e,%e,%e\n",f,x0,x);sleep(1);
		if(abs11(f)<1.0e-20){
			break;
		}else{	
			f1=2.0*dpz-bpx-(b*bpx-2.0*cpx)/sqrt(b*b-4.0*c);
			if(abs11(f1)<1.0e-20){
				break;
			}	
			x=x0-f/f1;
		}
		nn++;
	}while(abs11(x-x0)>1e-11);//printf("f1=%e\n",f1);
	return x;	
}
//sovle qsv tti slowness equation using the newton method in the local solver
double  Solve_local_traveltime_news2(double   Tim, int cell_num,double   Tad[8],double   d,double   alphy, double   beta,
                           double   deta,double   epsilon,double   Ragl,double  downt,double   upt)
{
	double  pz,px;
	double   cp0,cp1,cp2,cp3,cp4;
	double   Ta,Tb;
	double   mid1,mid2;
	double   root[2];
	double   A,B,C,D;
	double   tv,ts;
	int flag=0;
	double   sa,ca,sa2,ca2,sa3,ca3,sa4,ca4;
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

	tv=1.0/alphy/alphy;
	ts=1.0/beta/beta;

	A=1+2*epsilon;
	B=-(tv+ts);
	C=-(tv+(1+2*epsilon)*ts);
	D=2*(1+deta+(epsilon-deta)*ts/tv);
	sa=sin(Ragl);
	ca=cos(Ragl);

    double   T_beg;
	double   md1,md2,V1,V2,Px,Pz;
	double   Tc,T_c;
	double   Tcc;
	double T_end,tup,tdown,fv;
	double td;
	Tcc=Tim;
	double val=1.0e-6;
	int tag,nn_s;
	double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3;
	double da;
	double da1,da2,da3;
	double a,b,c;
	T_beg=min(Ta,Tb);
	T_end=upt;
	tup=T_end;tdown=T_beg;
    if(cell_num==7||cell_num==4){
		if(cell_num==7){
			

			b1=tup;
			b2=(tup+tdown)/2.0;
			b3=tdown;
		
			T_c=b1;
			px=-(Tb-Ta)/d;
			pz=(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;

			T_c=b2;
			px=-(Tb-Ta)/d;
			pz=(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d2=fv;

			T_c=b3;
			px=-(Tb-Ta)/d;
			pz=(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-----------------------------------------------------------------
			//newton method
			px=-(Tb-Ta)/d;
	
			T_c=0.5*(T_end+T_beg);
			//if(d1*d3>0){
			//	if(d1*d2<0 || d2*d3<0)
			//	printf("d1*d2,d2*d3=%e,%e\n",d1*d2,d2*d3);
			//}
			if(d1*d3<0){
				pz=newton(px,(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);//printf("netown\n");	
			}
			md1=px;
			md2=pz;
			Px=md1*ca+md2*sa;
			Pz=-md1*sa+md2*ca;
			Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
			md1=V1;
			md2=V2;
			V1=md1*ca-md2*sa;
			V2=md1*sa+md2*ca;
			Test_causality(V1,V2,cell_num,&flag);
			if(flag){
				Tc=-pz*d+Tb;
				Tcc=min(Tc,Tim);	
			}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}else if(cell_num==4){
			T_beg=min(Ta,Tb);
			T_end=upt;
			tup=T_end;tdown=T_beg;

			b1=tup;
			b2=(tup+tdown)/2.0;
			b3=tdown;
		
			T_c=b1;
			px=-(Tb-Ta)/d;
			pz=-(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;

			T_c=b3;
			px=-(Tb-Ta)/d;
			pz=-(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-------------
			//newton method
			px=-(Tb-Ta)/d;
			T_c=0.5*(T_end+T_beg);
			//if(d1*d3<0){
			if(d1*d3<0)
				pz=newton(px,-(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
			md1=px;
			md2=pz;
			Px=md1*ca+md2*sa;
			Pz=-md1*sa+md2*ca;
			Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
			md1=V1;
			md2=V2;
			V1=md1*ca-md2*sa;
			V2=md1*sa+md2*ca;
			Test_causality(V1,V2,cell_num,&flag);
			if(flag){
				Tc=pz*d+Tb;
				Tcc=min(Tc,Tim);	
			}	
				if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
				return Tcc;
		}
		
    }
    if(cell_num==1||cell_num==6){
	    if(cell_num==1){ 
			//--------------------------------------
			T_beg=min(Ta,Tb);
			T_end=upt;
			tup=T_end;tdown=T_beg;
			b1=tup;
			b2=(tup+tdown)/2.0;
			b3=tdown;

	
			T_c=b1;
			px=(T_c-Tb)/d;
			pz=(Ta-Tb)/d;
			mid1=px;
			mid2=pz;
			px=mid1*ca+mid2*sa;
			pz=-mid1*sa+mid2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;

			T_c=b3; 
			px=(T_c-Tb)/d;
			pz=(Ta-Tb)/d;
			mid1=px;
			mid2=pz;
			px=mid1*ca+mid2*sa;
			pz=-mid1*sa+mid2*ca;
			Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-------------
			//newton method
				pz=(Ta-Tb)/d;
			
				T_c=0.5*(T_end+T_beg);
				if(d1*d3<0)
				px=newton_px(pz,-(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=px*d+Tb;
					Tcc=min(Tc,Tim);	
				}
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			return Tcc;//*/
		 }else if(cell_num==6){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
		
				T_c=b1;
				px=-(T_c-Tb)/d;
				pz=(Ta-Tb)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3; 
				px=-(T_c-Tb)/d;
				pz=(Ta-Tb)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;
				//newton method
				pz=(Ta-Tb)/d;
			
				T_c=0.5*(T_end+T_beg);
				if(d1*d3<0)
				px=newton_px(pz,(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=-px*d+Tb;
					Tcc=min(Tc,Tim);	
				}
				if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
				return Tcc;
		 }
    }
    if(cell_num==0||cell_num==3){
		if(cell_num==0){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
	
				T_c=b1;
				px=(Tb-Ta)/d;
				if(cell_num==0){
					pz=(Tb-T_c)/d;
				}else{
					pz=-(Tb-T_c)/d;
				}
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3;
				px=(Tb-Ta)/d;
				if(cell_num==0){
					pz=(Tb-T_c)/d;
				}else{
					pz=-(Tb-T_c)/d;
				}
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;
					//-----------------------------------------------------------------
				//newton method
				px=(Tb-Ta)/d;
			
				T_c=0.5*(T_end+T_beg);
				//if(d1*d3<0){
				if(d1*d3<0)
					pz=newton(px,(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=-pz*d+Tb;
					Tcc=min(Tc,Tim);	
				}
			 	if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			 return Tcc;
		}else if(cell_num==3){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;

				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
	
				T_c=b1;
				px=(Tb-Ta)/d;
				if(cell_num==0){
					pz=(Tb-T_c)/d;
				}else{
					pz=-(Tb-T_c)/d;
				}
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3;
				px=(Tb-Ta)/d;
				if(cell_num==0){
					pz=(Tb-T_c)/d;
				}else{
					pz=-(Tb-T_c)/d;
				}
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;
				//-----------------------------------------------------------------
				//newton method
				px=(Tb-Ta)/d;
				if(d3>d1){
					T_c=T_end;
				}else{
					T_c=T_beg;
				}
				T_c=0.5*(T_end+T_beg);
				//if(d1*d3<0){
				if(d1*d3<0)
					pz=newton(px,-(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=pz*d+Tb;
					Tcc=min(Tc,Tim);	
				}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}
	}
	if(cell_num==2||cell_num==5){
	    if(cell_num==2){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
	
				T_c=b1;
	
				px=(T_c-Tb)/d;
	
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3;
				px=(T_c-Tb)/d;
				
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;

				//newton method
				pz=-(Ta-Tb)/d;
				
				T_c=0.5*(T_end+T_beg);
				if(d1*d3<0)
				px=newton_px(pz,-(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=px*d+Tb;
					Tcc=min(Tc,Tim);	
				}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}else if(cell_num==5){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
				T_c=b1;
			
				px=-(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;


				T_c=b3;
				px=-(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;

				//newton method
				pz=-(Ta-Tb)/d;
				T_c=0.5*(T_end+T_beg);
				if(d1*d3<0)
				px=newton_px(pz,(Tb-T_c)/d,alphy,beta,deta,epsilon,Ragl);
				md1=px;
				md2=pz;
				Px=md1*ca+md2*sa;
				Pz=-md1*sa+md2*ca;
				Cal_group_pxz_new(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
				md1=V1;
				md2=V2;
				V1=md1*ca-md2*sa;
				V2=md1*sa+md2*ca;
				Test_causality(V1,V2,cell_num,&flag);
				if(flag){
					Tc=-px*d+Tb;
					Tcc=min(Tc,Tim);	
				}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}
    }
}
 void firstorder_traveltime(int nx,int nz,double dd,int sx,int sz,double **alpha,double **belta,
			    double **delta,double **epsilon,double **theta,double **T,int width,int total)
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
	double (*fix_vg)[8];
	FILE *fp;


	fix_vg=(double(*)[8])calloc(nx*nz,sizeof(*fix_vg));
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
				angle=Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}else{
			for(j=sz-width+1;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}
    }//exit(0);
	/*for(i=2;i<nx-2;i++){
		for(j=2;j<nz-2;j++){  
			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			deta=delta[i][j];
			eps=epsilon[i][j];

			cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			for(k=0;k<8;k++){
				fix_vg[i*nz+j][k]=Vag[k];
			}
		}
	}*/
	inputnum=0;
	start: inputnum++;
	printf("%s%d\n", "interation number begin:",inputnum);

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
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			cal_group_velocity_vs(Vag,alphy,beta,eps,deta,Ragl);
			/*for(k=0;k<8;k++){
				Vag[k]=fix_vg[i*nz+j][k];
			}//*/
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
			downt=-10;
		    upt=10;

      		for(num=4;num<6;num++){ 
				if(num==4&&T[i+1][j]>=T[i][j]&&T[i+1][j+1]>=T[i][j]){
					continue;
				}else if(num==5&&T[i+1][j+1]>=T[i][j]&&T[i][j+1]>=T[i][j]){
					continue;
				}
				T[i][j]=Solve_local_traveltime_news2(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
		     
				
     	 	}
        }   
  	}
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
				deta=delta[i][j];
				eps=epsilon[i][j];
				sina=sin(Ragl);
				cosa=cos(Ragl);
  
				//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
				cal_group_velocity_vs(Vag,alphy,beta,eps,deta,Ragl);
				/*for(k=0;k<8;k++){
					Vag[k]=fix_vg[i*nz+j][k];
				}//*/
				Tad[0]=T[i-1][j];
				Tad[1]=T[i-1][j+1];
				Tad[2]=T[i][j+1];
				Tad[3]=T[i+1][j+1];
				Tad[4]=T[i+1][j];
				Tad[5]=T[i+1][j-1];
				Tad[6]=T[i][j-1];
				Tad[7]=T[i-1][j-1]; 

				T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

				downt=-10;
				upt=T[i][j];
				upt=10;
		
				for(num=0;num<2;num++){ 
				 	if(num==0&&T[i-1][j-1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
				 	    continue;
				 	}else if(num==1&&T[i-1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
				 		continue;
		     		}
					T[i][j]=Solve_local_traveltime_news2(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
					
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
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			cal_group_velocity_vs(Vag,alphy,beta,eps,deta,Ragl);
			/*for(k=0;k<8;k++){
				Vag[k]=fix_vg[i*nz+j][k];
			}//*/
			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];
   
			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

			downt=-10;
			upt=T[i][j];
			upt  =10;

			for(num=2;num<4;num++){ 
		     	if(num==2&&T[i+1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
		     		continue;
		     	}else if(num==3&&T[i+1][j-1]>=T[i][j]&&T[i+1][j]>=T[i][j]){
		     		continue;
		     	}
				
				T[i][j]=Solve_local_traveltime_news2(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				
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
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			cal_group_velocity_vs(Vag,alphy,beta,eps,deta,Ragl);
			/*for(k=0;k<8;k++){
				Vag[k]=fix_vg[i*nz+j][k];
			}//*/
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

			downt=-10;
			upt=10;

	        for(num=6;num<8;num++){ 
     	 		if(num==6&&T[i][j+1]>=T[i][j]&&T[i-1][j+1]>=T[i][j]){
     	 			continue;
     	 		}else if(num==7&&T[i-1][j+1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
     	 			continue;
     	 		}
				T[i][j]=Solve_local_traveltime_news2(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				
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
		//printf("first order iteration %d: %f  %f\n",inputnum,sum,sum/(nx-4)/(nz-4));	
  	}else{
		sum=0;
	      for(j=2;j<nx-2;j++){
			  for(i=2;i<nz-2;i++){
					sum=sum+abs11(Tpre[j][i]-T[j][i]);
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
	total=2;
	if(inputnum<total){
		goto start;
	}
	  
}

