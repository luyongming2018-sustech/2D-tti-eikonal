#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double abs11(double a)
{ 
	if(a<0){
		return -a;
	}else{
		return a;
	}
}
double max11(double a,double b,double c,double d)
{
  	double mid;

	if(a>b){
		mid=a;
	}else{ 
    	mid=b;
	}

	if(mid<c){
		mid=c;
	}
	if(mid<d){
		mid=d;
	}

	return mid;
}
double min4(double a,double b,double c,double d)
{
  	double mid;

	if(a<b){
		mid=a;
	}else{ 
    	mid=b;
	}

	if(mid>c){
		mid=c;
	}
	if(mid>d){
		mid=d;
	}

	return mid;
}
double min(double a,double b)
{
      if(a<b){
            return a;
      }else{ 
            return b;
      }
}
double max_2(double a,double b)
{
      if(a>b){
            return a;
      }else{ 
            return b;
      }
}
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
void sort(double a[3])
{
    double mid;

    if(a[1]<a[0]){
        mid=a[0];
        a[0]=a[1];
        a[1]=mid;
    }
    if(a[2]<a[0]){
        mid=a[0];
        a[0]=a[2];
        a[2]=mid;
    }
    if(a[2]<a[1]){
        mid=a[1];
        a[1]=a[2];
        a[2]=mid;
    }
}

//---------------------------------------------------------------------------------------------------------
void Cal_fun(double a,double b,double c,double d,double e,double x, double*fx)
{
    *fx=a*x*x*x*x+b*x*x*x+c*x*x+d*x+e;
}
//er fen fa
double Bisection(double db,double ub,double a,double b,double c,double d,double e)
{
    double up; 
    double down; 
    double mid;
    double fv;
    int num,i;

    num=0;
    up=ub;
    down=db;

    Cal_fun(a,b,c,d,e,down,&fv);
  
    if(fv<0.0){
        up=ub;
        down=db;
    }else{
        up=db;
        down=ub;
    }
    while(abs11(down-up)>1.0e-8){
	//for(i=0;i<limit;i++){
        mid=(up+down)/2.0;
        Cal_fun(a,b,c,d,e,mid,&fv);
        if(fv<0.0){
            down=mid;
        }else{
            up=mid;
        }
    }
    return mid;
}
//solve san ci fang cheng
int Solcub(double a, double b, double c,double sol[3])
{
    double Q,R;
    double mid;
    double theta;
    double pi=3.14159265359;
    double mid1;

    Q=(a*a-3*b)/9;
    R=(2*a*a*a-9*a*b+27*c)/54;

    mid=Q*Q*Q-R*R;
  
    if(mid<0){
        sol[0]=0;
        sol[1]=0;
        sol[2]=0;
        return 1;
    }else{
        mid1=R/sqrt(Q*Q*Q);
        theta=acos(mid1);

        sol[0]=-2*sqrt(Q)*cos(theta/3)-a/3;
        sol[1]=-2*sqrt(Q)*cos((theta+2*pi)/3)-a/3;
        sol[2]=-2*sqrt(Q)*cos((theta-2*pi)/3)-a/3;

        sort(sol);
        return 0;
    }

}
//solve si ci equation using numerical method
int Solquar(double downt, double upt,double a, double b, double c,double d,double e,double root[2])
{

    double pa,pb,pc,pd;
    int flag;
    double sol[3];
    double mid1,mid2,mid3;
    double mida,midb;
    double aa,bb;
    double x1,x2,x3;

    pa=a*4;
    pb=b*3;
    pc=c*2;
    pd=d;

    root[0]=0;
    root[1]=0;


    pb=pb/pa;
    pc=pc/pa;
    pd=pd/pa;
    aa=downt;
    bb=upt;

    flag=Solcub(pb,pc,pd,sol);

//  root number of the cubic equation is less than 3   
    if(flag !=0){return 0;}

    Cal_fun(a,b,c,d,e,sol[0],&mid1);
    Cal_fun(a,b,c,d,e,sol[1],&mid2);
    Cal_fun(a,b,c,d,e,sol[2],&mid3);
    Cal_fun(a,b,c,d,e,aa,&mida);
    Cal_fun(a,b,c,d,e,bb,&midb);
    x1= sol[0];
    x2= sol[1];
    x3= sol[2];

//  root number is equal to 3

    if(abs11(mid1)<=1.0e-20 && mid3<0){
       if(midb<0 && mida>0){
            root[0] = Bisection(x2,x3,a,b,c,d,e);
            root[1] = 0;
            return 1;
       }
       if(bb>x3 && midb>0 && aa<x3 && mida>=0){
            root[0] = Bisection(x2,x3,a,b,c,d,e);
            root[1] = 0;
            return 1;
       }
    }

    if(abs11(mid2)<=1.0e-20){
       if(aa<sol[1] && bb>sol[1]){
            root[0] = sol[1];
            root[1] = 0;
            return 1;
       }
    }
    if(abs11(mid3)<=1.0e-20 && mid1<0){
       if( mida<0 && midb>0){
            root[0] = Bisection(x1,x2,a,b,c,d,e);
            root[1] = 0;
            return 1;
       }
       if(aa<x1 && mida>0 && bb>x1 && midb>=0){
            root[0] = Bisection(x1,x2,a,b,c,d,e);
            root[1] = 0;
            return 1;
       }
    }

//  root number is equal to 4
    if(mid1<0 && mid3<0){
        if(aa>x1 && aa<x3 && mida>0){
            if(midb<0 || bb>x3){
                root[0] = Bisection(x2,x3,a,b,c,d,e);
                root[1] = 0;
                return 1;
            }
        }

        if(bb>x1 && bb<x3 && midb>0){
            if(mida<0 || aa<x1){
                root[0] = Bisection(x1,x2,a,b,c,d,e);
                root[1] = 0;
                return 1;
            }
        }
        if(mida<0 && aa<x2){
            if(midb<0 && bb>x2){
                root[0] = Bisection(x1,x2,a,b,c,d,e);
                root[1] = Bisection(x2,x3,a,b,c,d,e);
                return 2;
            }
            if(bb>x3){
                root[0] = Bisection(x1,x2,a,b,c,d,e);
                root[1] = Bisection(x2,x3,a,b,c,d,e);
                return 2;
            }
        }

        if(mida>0 && aa<x1){
            if(midb<0 && bb>x2){
                root[0] = Bisection(x1,x2,a,b,c,d,e);
                root[1] = Bisection(x2,x3,a,b,c,d,e);
                return 2;
            }
            if(bb>x3){
                root[0] = Bisection(x1,x2,a,b,c,d,e);
                root[1] = Bisection(x2,x3,a,b,c,d,e);
                return 2;
            }
        }
    }
     return 0;
}
void Wave_surfacetest(double px,double pz,double alphy,double beta,double delta,double epsilon, int *flag)
{
	double a,b,c;
	double value;
	double tb,terr; 

	a=1.0/(alphy*alphy);
	tb=1.0/(beta*beta);
	b=-2*(1+delta+(epsilon-delta)*alphy*alphy*tb)*px*px+a+tb;
	c=((1+2*epsilon)*px*px-a)*(px*px-tb);
	value=b/2-sqrt(b*b-4*c)/2;

	terr=abs11(pz*pz-value);

	if(terr<1.0e-9){
		*flag=1;
	}else{
		*flag=0;
	}
}
void Wave_surfacetest_vp(double  px,double  pz,double  alphy,double  beta,double  delta,double  epsilon, int *flag,double *fv)
{
	double  a,b,c;
	double  value;
	double  tb,terr; 

	a=1.0/(alphy*alphy);
	tb=1.0/(beta*beta);
	b=-2*(1+delta+(epsilon-delta)*alphy*alphy*tb)*px*px+a+tb;
	c=((1+2*epsilon)*px*px-a)*(px*px-tb);
	value=b/2.0-sqrt(b*b-4.0*c)/2.0;

	terr=abs11(pz*pz-value);
	*fv=pz*pz-value;
	if(terr<1.0e-6){
		*flag=1;
	}else{
		*flag=0;
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
void Wave_surfacetest_vs(double  px,double  pz,double  alphy,double  beta,double  delta,double  epsilon, int *flag, double *fv)
{
	double  a,b,c;
	double  value;
	double  tb,terr; 

	a=1.0/(alphy*alphy);
	tb=1.0/(beta*beta);
	b=-2*(1+delta+(epsilon-delta)*alphy*alphy*tb)*px*px+a+tb;
	c=((1+2*epsilon)*px*px-a)*(px*px-tb);
	value=b/2+sqrt(b*b-4*c)/2;
	*fv=pz*pz-value;
	terr=abs11(pz*pz-value);

	if(terr<1.0e-9){
		*flag=1;
	}else{
		*flag=0;
	}
}
void Wave_surfacetest_sh(double  px,double  pz,double  alphy,double  beta,double  delta,double  epsilon, int *flag)
{
	double  gama=0.4;
	double  value;
	double  tb,terr; 

	tb=1.0/(beta*beta);
	terr=abs11(pz*pz+(2*gama+1)*px*px-tb);

	if(terr<1.0e-9){
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
	sqqr=sqrt(B*B-4*C);
	slowppx=-0.5*partialB-1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
}
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
		value=b-sqrt(b*b-4.0*c);
		if(b*b-4.0*c<0)printf("b*b-4.0*c=%e\n");
		f=2.0*pzz*pzz-value;
		if(abs11(f)<1.0e-20){
			break;
		}else{	
			f1=2.0*dpz-bpz+(b*bpz-2.0*cpz)/sqrt(b*b-4.0*c);
			if(abs11(f1)<1.0e-20){
				break;
			}	
			x=x0-f/f1;
		}
		nn++;
	}while(abs11(x-x0)>1e-6);//printf("f1=%e\n",f1);
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
		if(b*b-4.0*c<0)printf("b*b-4.0*c=%e\n");
		value=b-sqrt(b*b-4.0*c);

		f=2.0*pzz*pzz-value;//printf("f,px,pz,x0,x=%e,%e,%e,%e,%e\n",f,px,pz,x0,x);sleep(1);
		//printf("f,x0,x=%e,%e,%e\n",f,x0,x);sleep(1);
		if(abs11(f)<1.0e-20){
			break;
		}else{	
			f1=2.0*dpz-bpx+(b*bpx-2.0*cpx)/sqrt(b*b-4.0*c);
			if(abs11(f1)<1.0e-20){
				break;
			}	
			x=x0-f/f1;
		}
		nn++;
	}while(abs11(x-x0)>1e-6);//printf("f1=%e\n",f1);
	//printf("iter =%d\n",nn-1);
	return x;	
}
double  Solve_local_traveltime_news(double   Tim, int cell_num,double   Tad[8],double   d,double   alphy, double   beta,
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
    if(cell_num==7||cell_num==4){
		if(cell_num==7){
			T_beg=max_2(Ta,Tb);
			T_end=upt;
			tup=T_end;tdown=T_beg;

			b1=tup;
			b2=(tup+tdown)/2.0;
			b3=tdown;
		
			T_c=b1;
			px=-(Tb-Ta)/d;
			pz=(Tb-T_c)/d;//printf("pz1=%e\n",pz);
			a=pz;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;


			T_c=b3;
			px=-(Tb-Ta)/d;
			pz=(Tb-T_c)/d;
			b=pz;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-----------------------------------------------------------------
			//newton method
			px=-(Tb-Ta)/d;
			T_c=0.5*(T_end+T_beg);
			//if(d1*d3>0)printf("pz1,pz2=%e,%e\n",a,b);sleep(1);
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
				Tcc=min(Tc,Tcc);	
			}
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			return Tcc;
		}else if(cell_num==4){
			T_beg=max_2(Ta,Tb);
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
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;

			T_c=b3;
			px=-(Tb-Ta)/d;
			pz=-(Tb-T_c)/d;
			md1=px;
			md2=pz;
			px=md1*ca+md2*sa;
			pz=-md1*sa+md2*ca;
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-------------
			//newton method
			px=-(Tb-Ta)/d;
			if(d3>d1){
				T_c=T_end;
			}else{
				T_c=T_beg;
			}
			T_c=0.5*(T_end+T_beg);
			//if(d1*d3<0){
			if(d1*d3<0)
				pz=newton(px,(T_c-Tb)/d,alphy,beta,deta,epsilon,Ragl);

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
				if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
				return Tcc;
		}
		
    }
    if(cell_num==1||cell_num==6){
	    if(cell_num==1){ 
			//--------------------------------------
			T_beg=max_2(Ta,Tb);
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
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d1=fv;

			T_c=b3; 
			px=(T_c-Tb)/d;
			pz=(Ta-Tb)/d;
			mid1=px;
			mid2=pz;
			px=mid1*ca+mid2*sa;
			pz=-mid1*sa+mid2*ca;
			Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
			d3=fv;
			//-------------
			//newton method
				pz=(Ta-Tb)/d;
				if(d3>d1){
					T_c=T_end;
				}else{
					T_c=T_beg;
				}
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
				}//*/
				/*while(abs11(tup-tdown)>val){
				T_c=(tup+tdown)*0.5;
				px=(T_c-Tb)/d;pz=(Ta-Tb)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Cal_group_pxz_new(px,pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Px=(T_c-Tb)/d;
					if(flag){
						Tc=Px*d+Tb;
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}//*/
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			//if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;//*/
		 }else if(cell_num==6){
				T_beg=max_2(Ta,Tb);
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
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3; 
				px=-(T_c-Tb)/d;
				pz=(Ta-Tb)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;
				//newton method
				pz=(Ta-Tb)/d;
				if(d3>d1){
					T_c=T_end;
				}else{
					T_c=T_beg;
				}
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
				if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
				return Tcc;
		 }
    }
    if(cell_num==0||cell_num==3){
		if(cell_num==0){
				T_beg=max_2(Ta,Tb);
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
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
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
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
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
				}//*/
				/*while(abs11(tup-tdown)>val){
					T_c=(tup+tdown)*0.5;
					px=(Tb-Ta)/d;
					pz=-(T_c-Tb)/d;
					md1=px;
					md2=pz;

					px=md1*ca+md2*sa;
					pz=-md1*sa+md2*ca;
	
					Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
					if(fv<0.0){
						tdown=T_c;
					}else{
						tup=T_c;
					}
					if(tag==1){
						Cal_group_pxz_new(px,pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
						md1=V1;
						md2=V2;
						V1=md1*ca-md2*sa;
						V2=md1*sa+md2*ca;
						Test_causality(V1,V2,cell_num,&flag);
						Pz=-(T_c-Tb)/d;
						if(flag){
							Tc=-Pz*d+Tb;
							Tcc=min(Tc,Tim);
							//break;
						}
				
					}
				}//*/
			 	//if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
				if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
				 return Tcc;
		}else if(cell_num==3){
				T_beg=max_2(Ta,Tb);
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
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
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
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
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
				}//*/
			
				/*while(abs11(tup-tdown)>val){
				T_c=(tup+tdown)*0.5;
				px=(Tb-Ta)/d;
				pz=(T_c-Tb)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Cal_group_pxz_new(px,pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Pz=(T_c-Tb)/d;
					if(flag){
						Tc=Pz*d+Tb;//printf("Ta,Tb,Tc,Tim=%e,%e,%e,%e\n",Ta,Tb,Tc,Tim);sleep(1);
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}//*/
			//if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			return Tcc;
		}
	}
	if(cell_num==2||cell_num==5){
	    if(cell_num==2){
				T_beg=max_2(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
	
				T_c=b1;
				if(cell_num==2){
					px=(T_c-Tb)/d;
				}else{
					px=-(T_c-Tb)/d;
				}
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;

				T_c=b3;
				if(cell_num==2){
					px=(T_c-Tb)/d;
				}else{
					px=-(T_c-Tb)/d;
				}
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;

				//newton method
				pz=-(Ta-Tb)/d;
				if(d3>d1){
					T_c=T_end;
				}else{
					T_c=T_beg;
				}
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
				}//*/
				/*nn_s=0;T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				while(abs11(tup-tdown)>val){
				T_c=(tup+tdown)*0.5;
				px=(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Px=px;Pz=pz;
					Cal_group_pxz_new(px,pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Px=(T_c-Tb)/d;
					if(flag){
						Tc=Px*d+Tb;
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}//*/
			//if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			return Tcc;
		}else if(cell_num==5){
				T_beg=max_2(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				b1=tup;
				b2=(tup+tdown)/2.0;
				b3=tdown;
				T_c=b1;
				if(cell_num==2){
					px=(T_c-Tb)/d;
				}else{
					px=-(T_c-Tb)/d;
				}
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d1=fv;


				T_c=b3;
				if(cell_num==2){
					px=(T_c-Tb)/d;
				}else{
					px=-(T_c-Tb)/d;
				}
				pz=(Tb-Ta)/d;
				mid1=px;
				mid2=pz;
				px=mid1*ca+mid2*sa;
				pz=-mid1*sa+mid2*ca;
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);
				d3=fv;

				//newton method
				pz=-(Ta-Tb)/d;
				if(d3>d1){
					T_c=T_end;
				}else{
					T_c=T_beg;
				}
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
				}//*/
				/*nn_s=0;T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
				while(abs11(tup-tdown)>val){
				T_c=(tup+tdown)*0.5;
				px=-(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vp(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Cal_group_pxz_new(px,pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Px=-(T_c-Tb)/d;
					if(flag){
						Tc=-Px*d+Tb;
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}//*/
			//if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			if(Tcc<Ta&&Tcc<Tb){Tcc=Tim;}
			return Tcc;
		}
    }
}
double  Solve_local_traveltime_modify_tti(double   Tim, int cell_num,double   Tad[8],double   d,double   alphy, double   beta,
                           double   deta,double   epsilon,double   Ragl,double  downt,double   upt)
{
	double  pz,px;
	double   cp0,cp1,cp2,cp3,cp4;
	double   Ta,Tb;
	double   mid1,mid2;
	double   root[2];
	double   A,B,C,D;
	long double   tv,ts;
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


	sa=sin(Ragl);
	ca=cos(Ragl);

    double   T_beg;
	double   md1,md2,V1,V2,Px,Pz;
	double   Tc,T_c;
	double   Tcc;
	double T_end,tup,tdown,fv;
	Tcc=Tim;
	int tag=100;
	long double b,c,value;
	double px1,pz1;
    if(cell_num==7||cell_num==4){
				px=-(Tb-Ta)/d;
				//pz=(Tb-T_c)/d;
				b=-2*(1+deta+(epsilon-deta)*alphy*alphy*ts)*px*px+tv+ts;
				c=((1+2*epsilon)*px*px-tv)*(px*px-ts);
				if((b*b-4.0*c)>0.0){
					value=b*0.5+sqrt(b*b-4.0*c)*0.5;
				}else{
					value=0;
				}
				if(value>0){
					root[0]=-sqrt(value);
					root[1]=sqrt(value);
				}else{
					root[0]=0;
					root[1]=0;
				}
				
				for(i=0;i<2;i++){
					pz=root[i];
					md1=px;
					md2=pz;

					px1=md1*ca+md2*sa;
					pz1=-md1*sa+md2*ca;
					//px1=md1*ca-md2*sa;
					//pz1=md1*sa+md2*ca;
					Wave_surfacetest_vs(px1,pz1,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
					if(tag==1){
						Px=px1;Pz=pz1;
						Cal_group_pxz_new_vs(px1,pz1,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
						md1=V1;
						md2=V2;
						V1=md1*ca-md2*sa;
						V2=md1*sa+md2*ca;
						Test_causality(V1,V2,cell_num,&flag);
						if(flag){
							if(cell_num==7){
								Tc=-pz*d+Tb;
							}else if(cell_num==4){
								Tc=pz*d+Tb;
							}
							Tcc=min(Tc,Tim);//printf("Tc=%e\n",Tc);printf("i,px,pz=%d,%e,%e\n",i,Px,Pz); //sleep(1);
						}
	
					}
				}
				if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
				return Tcc;
	}
	 if(cell_num==1||cell_num==6){
				pz=(Ta-Tb)/d;
				b=tv+(1.0+2.0*epsilon)*ts-2.0*(1.0+deta+(epsilon-deta)*alphy*alphy*ts)*pz*pz;
				b=b/(1.0+2.0*epsilon);
				c=pz*pz*pz*pz+tv*ts-(tv+ts)*pz*pz;
				c=c/(1.0+2.0*epsilon);
				if((b*b-4.0*c)>0.0){
					value=b*0.5+sqrt(b*b-4.0*c)*0.5;
				}else{
					value=0;
				}
				if(value>0.0){
					root[0]=-sqrt(value);//printf("px=%f\n",pz);sleep(1);
					root[1]=sqrt(value);
				}else{
					root[0]=0;
					root[1]=0;
				}
				
				for(i=0;i<2;i++){
					px=root[i];
					md1=px;
					md2=pz;
		
					px1=md1*ca+md2*sa;
					pz1=-md1*sa+md2*ca;
					//px1=md1*ca-md2*sa;
					//pz1=md1*sa+md2*ca;
					Wave_surfacetest_vs(px1,pz1,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
					if(tag==1){
						Px=px1;Pz=pz1;
						Cal_group_pxz_new_vs(px1,pz1,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
						md1=V1;
						md2=V2;
						V1=md1*ca-md2*sa;
						V2=md1*sa+md2*ca;
						Test_causality(V1,V2,cell_num,&flag);
						if(flag){
							if(cell_num==1){
								Tc=px*d+Tb;
							}else{
								Tc=-px*d+Tb;
							}
							Tcc=min(Tc,Tim);
						}
				
					}
				}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;//*/
    }
    if(cell_num==0||cell_num==3){
			px=(Tb-Ta)/d;
			//pz=-(T_c-Tb)/d;
			b=-2*(1+deta+(epsilon-deta)*alphy*alphy*ts)*px*px+tv+ts;
			c=((1+2*epsilon)*px*px-tv)*(px*px-ts);
			if((b*b-4.0*c)>0.0){
				value=b*0.5+sqrt(b*b-4.0*c)*0.5;
			}else{
				value=0;
			}
	
			if(value>0.0){
				root[0]=-sqrt(value);//printf("px=%f\n",pz);sleep(1);
				root[1]=sqrt(value);
			}else{
				root[0]=0;
				root[1]=0;
			}
			for(i=0;i<2;i++){
				pz=root[i];
				md1=px;
				md2=pz;

				px1=md1*ca+md2*sa;
				pz1=-md1*sa+md2*ca;
				//px1=md1*ca-md2*sa;
				//pz1=md1*sa+md2*ca;

				Wave_surfacetest_vs(px1,pz1,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(tag==1){
					Px=px1;Pz=pz1;
					Cal_group_pxz_new_vs(px1,pz1,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					if(flag){
						if(cell_num==0){
							Tc=-pz*d+Tb;
						}else{
							Tc=pz*d+Tb;
						}
						Tcc=min(Tc,Tim);
					}
				}
			}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
	
	}
	/*if(cell_num==2||cell_num==5){
				//px=(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				b=tv+(1.0+2.0*epsilon)*ts-2.0*(1.0+deta+(epsilon-deta)*alphy*alphy*ts)*pz*pz;
				b=b/(1.0+2.0*epsilon);
				c=pz*pz*pz*pz+tv*ts-(tv+ts)*pz*pz;
				c=c/(1.0+2.0*epsilon);
				if((b*b-4.0*c)>0.0){
					value=b*0.5+sqrt(b*b-4.0*c)*0.5;
				}else{
					value=0;
				}
				//printf("pz,b,c,val=%f,%e,%e,%e\n",pz,b,c,value);
				if(value>0.0){
					root[0]=-sqrt(value);
					root[1]=sqrt(value);
					//printf("r1,r2=%e,%e\n",root[0],root[1]);//sleep(1);
				}else{
					root[0]=0.0;
					root[1]=0.0;
				}
				for(i=0;i<2;i++){
					px=root[i];
					md1=px;
					md2=pz;

					px1=md1*ca+md2*sa;
					pz1=-md1*sa+md2*ca;
					//px1=md1*ca-md2*sa;
					//pz1=md1*sa+md2*ca;
					Wave_surfacetest_vs(px1,pz1,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
					if(tag==1){
						Px=px1;Pz=pz1;
						Cal_group_pxz_new_vs(px1,pz1,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
						md1=V1;
						md2=V2;
						V1=md1*ca-md2*sa;
						V2=md1*sa+md2*ca;
						Test_causality(V1,V2,cell_num,&flag);
						if(flag){
							if(cell_num==2){
								Tc=px*d+Tb;Tcc=min(Tc,Tcc);
								//printf("tc=%f\n",Tc);//sleep(1);
								//printf("r1,r2=%f,%f\n",root[0],root[1]);//sleep(1);
							}else if(cell_num==5){
								Tc=-px*d+Tb;Tcc=min(Tc,Tcc);
								//printf("tc=%f\n",Tc);//sleep(1);
								//printf("r1,r2=%f,%f\n",root[0],root[1]);//sleep(1);
							}
							
						}
				
					}
				}
			if(Tcc<Ta&& Tcc<Tb){Tcc=Tim;}
			return Tcc;
    }//*/
	if(cell_num==2||cell_num==5){
	    if(cell_num==2){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
			//for(i=0;i<n_m;i++){
			while(abs11(tup-tdown)>1.0e-6){
				T_c=(tup+tdown)*0.5;
				px=(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Px=px;Pz=pz;
					Cal_group_pxz_new_vs(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Px=(T_c-Tb)/d;
					if(flag){
						Tc=Px*d+Tb;
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}else if(cell_num==5){
				T_beg=min(Ta,Tb);
				T_end=upt;
				tup=T_end;tdown=T_beg;
			while(abs11(tup-tdown)>1.0e-10){
				T_c=(tup+tdown)*0.5;
				px=-(T_c-Tb)/d;
				pz=(Tb-Ta)/d;
				md1=px;
				md2=pz;

				px=md1*ca+md2*sa;
				pz=-md1*sa+md2*ca;
	
				Wave_surfacetest_vs(px,pz,alphy,beta,deta,epsilon,&tag,&fv);//printf("tag,fv=%d,%e\n",tag,fv);//sleep(1);
				if(fv<0.0){
					tdown=T_c;
				}else{
					tup=T_c;
				}
				if(tag==1){
					Px=px;Pz=pz;
					Cal_group_pxz_new_vs(Px,Pz,&V1,&V2,alphy,beta,deta,epsilon,cell_num);
					md1=V1;
					md2=V2;
					V1=md1*ca-md2*sa;
					V2=md1*sa+md2*ca;
					Test_causality(V1,V2,cell_num,&flag);
					Px=-(T_c-Tb)/d;
					if(flag){
						Tc=-Px*d+Tb;
						Tcc=min(Tc,Tim);
						//break;
					}
				
				}
			}
			if(Tcc<max_2(Ta,Tb)){Tcc=Tim;}
			return Tcc;
		}
    }//*/
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
void Cal_group_pxz_vp(double Px[2],double Pz[2], double * V1, double * V2,double alphy,double beta,double delta, double eps,int k,int num0)
{ 
    double pz,px;
	double partialB,partialC;
	double B,C,sqqr,slowppx,low;
	double tp,ts;

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
	slowppx=-0.5*partialB+1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
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
void Solve_local_time(int cell_num,double Tad[8],double Tc,double d,double *px,double *pz)
{
	double cp0,cp1,cp2,cp3,cp4;
	double Ta,Tb;
	double mid1,mid2;
	double root[2];
	double A,B,C,D;
	double tv,ts;
	int flag=0;
	double sa,ca,sa2,ca2,sa3,ca3,sa4,ca4;
	int i;

    switch(cell_num){ 
			case 0:{
				Ta=Tad[7];
				Tb=Tad[6];
				break;
			}
			case 1:{
				Ta=Tad[7];
				Tb=Tad[0];
				break;
			}
			case 2:{
				//Ta=Tad[5];
				//Tb=Tad[6];
				Ta=Tad[1];
				Tb=Tad[0];
				break;
			}
			case 3:{
				//Ta=Tad[5];
				//Tb=Tad[4];
				Ta=Tad[1];
				Tb=Tad[2];
				break;
			}
			case 4:{
				//Ta=Tad[3];
				//Tb=Tad[4];
				Ta=Tad[3];
				Tb=Tad[2];
				break;
			}
			case 5:{
				//Ta=Tad[3];
				//Tb=Tad[2];
				Ta=Tad[3];
				Tb=Tad[4];
				break;
			}
			case 6:{
				//Ta=Tad[1];
				//Tb=Tad[2];
				Ta=Tad[5];
				Tb=Tad[4];
				break;
			}
			case 7:{
				//Ta=Tad[1];
				//Tb=Tad[0];
				Ta=Tad[5];
				Tb=Tad[6];
				break;
			}
 	}


    if(cell_num==7||cell_num==4){
		*px=-(Tb-Ta)/d;
		if(cell_num==7){
			*pz=(Tc-Tb)/d;

		}else if(cell_num==4){
			*pz=(Tb-Tc)/d;
		}
    }
    if(cell_num==1||cell_num==6){
		*pz=(Tb-Ta)/d;
		if(cell_num==1){ 
			*px=(Tc-Tb)/d;
		 }else if(cell_num==6){
			*px=(Tb-Tc)/d;
		 }
	
    }
    if(cell_num==0||cell_num==3){
		*px=(Tb-Ta)/d;
		if(cell_num==0){
			*pz=(Tc-Tb)/d;
		}else if(cell_num==3){ 
			*pz=(Tb-Tc)/d;
		}
	}
	if(cell_num==2||cell_num==5){
		*pz=(Ta-Tb)/d;
		if(cell_num==2){ 
			*px=(Tc-Tb)/d;
		 }else if(cell_num==5){
			*px=(Tb-Tc)/d;
		 }
    }
}
//***************************************************************************************
void cal_group_velocity_vp(double Vag[8],double vp,double vs,double epsilon, double delta,double Ragl)
{ 
	int i;
	double pi=3.1415926,interval;
	double yita,Q,A,C,E,D,Vg;
	double len,Sn,Cs;

	interval=pi/4.0;
	for(i=0;i<8;i++){
		Sn=sin(interval*i+Ragl);
		Cs=cos(interval*i+Ragl);
		yita=(epsilon-delta)/(1.0+2.0*delta);
		C=1.0/vp/vp;
		Q=1.0+2.0*yita;
		A=1.0/vp/vp/(1.0+2.0*epsilon);
		E=A*Sn*Sn+C*Cs*Cs;
		len=E*E+4.0*(Q*Q-1.0)*A*C*Sn*Sn*Cs*Cs;
		D=sqrt(len);
		Vg=(1.0+2.0*Q)/2.0/(1.0+Q)*E+1.0/2.0/(1.0+Q)*D;
		Vg=1.0/Vg;
		Vg=sqrt(Vg);
		Vag[i]=Vg;
	}
}

void slowness_cal_groupvel_vp(double vp,double vs,double epsilon, double delta,double Px,double*vx,double*vz)
{
   double a,b,c;
   double solu,pxx[2],pzz[2];
   double Vx,Vz;

   a=1;
   b=-2.0*(1.0+delta+(epsilon-delta)*vp*vp/vs/vs)*Px*Px+1.0/vp/vp+1.0/vs/vs;
   c=1.0/vp/vp/vs/vs-(1.0/vp/vp+(1.0+2.0*epsilon)/vs/vs)*Px*Px+(1.0+2.0*epsilon)*Px*Px*Px*Px;

   solu=b/2.0-sqrt(b*b-4*a*c)/2.0;
   pxx[0]=Px;
   pzz[0]=sqrt(solu);
   if(solu<0){
   	 pzz[0]=0;
   }

   Cal_group_pxz_vp(pxx,pzz,&Vx,&Vz,vp,vs,delta,epsilon,0,1);

   *vx=Vx;
   *vz=Vz;
}
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
double Initialize_source_vs(double alphy,double beta,double delta,double epsilon,int sx,int sz,int ix,int jz,double d,double angle)
{

	double pxmin,pxmax,pmid;
	double vx,vz,an1;
	double pi=3.1415926535898;
	double angsim,ang;
	double Vp;
	double dis,tmp;
	int i;
	double dsr;

	//dsr=sqrt((ix-sx)*(ix-sx)+(jz-sz)*(jz-sz));
    //an1=cos(angle)*(ix-sx)/dsr+sin(angle)*(jz-sz)/dsr;
	an1=cos(angle)*(jz-sx)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz))+sin(angle)*(ix-sz)/sqrt((jz-sx)*(jz-sx)+(ix-sz)*(ix-sz));
	if(an1<0){
		an1=-an1;
	}

	ang=an1; 

	pxmin=0;
	//pxmax=1.0/sqrt(1+2*epsilon)/beta;
	pxmax=1.0/beta;
    for(i=0;i<1500;i++){
	//while(abs11(pxmax-pxmin)>1.0e-14){
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
double Initialize_source_vp(double alphy,double beta,double delta,double epsilon,int sx,int sz,int ix,int jz,double d,double angle)
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
	pxmax=1.0/sqrt(1+2*epsilon)/alphy;
    for(i=0;i<500;i++){
	//while(abs11(pxmax-pxmin)>1.0e-14){
		pmid=(pxmin+pxmax)/2.0;
		slowness_cal_groupvel_vp(alphy,beta,epsilon,delta,pmid,&vx,&vz);//printf("vx,vz=%e,%e\n",vx,vz);
		angsim=vx/sqrt(vx*vx+vz*vz);

		if(ang>angsim){
			pxmin=pmid;
		}else{
			pxmax=pmid;
		}   
	}

	Vp=sqrt(vx*vx+vz*vz);
	dis=d*dsr;//printf("vp,dis.vp=%e,%e\n",Vp,dis/Vp);

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
//compute fixed direction group velocity 
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
				angle=Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source_vp(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}else{
			for(j=sz-width+1;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source_vp(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}
    }//exit(0);

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
			cal_group_velocity_vp(Vag,alphy,beta,eps,deta,Ragl);
			//for(k=0;k<8;k++){
					//printf("vag[%d]=%f\n",k,Vag[k]);
					//printf("vag1[%d]=%f\n",k,Vag1[k]);
					//sleep(1);
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
			downt=-10;
		   // upt=10;

      		for(num=4;num<6;num++){ 
				if(num==4&&T[i+1][j]>=T[i][j]&&T[i+1][j+1]>=T[i][j]){
					continue;
				}else if(num==5&&T[i+1][j+1]>=T[i][j]&&T[i][j+1]>=T[i][j]){
					continue;
				}
				//T[i][j]=Solve_local_traveltime_modify_tti(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
		     	/*flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);

		     	if(flag==0){continue;}
     	
         		for(k=0;k<2;k++){ 
		            if(Tc[k]<0){continue;}
		     	    mid1=Px[k];
		     	    mid2=Pz[k];
		     	    Px[k]=mid1*cosa+mid2*sina;
		     	    Pz[k]=-mid1*sina+mid2*cosa;
		     	
		     	    Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);
		     	    if(flag==0){continue;}
		     
		            Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num);

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
				cal_group_velocity_vp(Vag,alphy,beta,eps,deta,Ragl);
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
				//upt=10;
		
				for(num=0;num<2;num++){ 
				 	if(num==0&&T[i-1][j-1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
				 	    continue;
				 	}else if(num==1&&T[i-1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
				 		continue;
		     		}
					//T[i][j]=Solve_local_traveltime_modify_tti(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
					T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
					//printf("Time=%e\n",T[i][j]);//sleep(1);
                	/*flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);

					if(flag==0){continue;} 
					for(k=0;k<2;k++){ 
						if(Tc[k]<0){continue;}
						mid1=Px[k];
						mid2=Pz[k];
		
						Px[k]=mid1*cosa+mid2*sina;
						Pz[k]=-mid1*sina+mid2*cosa;

				
						Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);

						if(flag==0){continue;}

						Cal_group_pxz_vp(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num);

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
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			cal_group_velocity_vp(Vag,alphy,beta,eps,deta,Ragl);
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
			//upt  =10;

			for(num=2;num<4;num++){ 
		     	if(num==2&&T[i+1][j-1]>=T[i][j]&&T[i][j-1]>=T[i][j]){
		     		continue;
		     	}else if(num==3&&T[i+1][j-1]>=T[i][j]&&T[i+1][j]>=T[i][j]){
		     		continue;
		     	}
				
				//T[i][j]=Solve_local_traveltime_modify_tti(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				/*flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);

				if(flag==0){continue;}  //flag=0 means we cannot derive a reasonable 
     	
				for(k=0;k<2;k++){ 
					if(Tc[k]<0){
						continue;
					}
					mid1=Px[k];
					mid2=Pz[k];
					Px[k]=mid1*cosa+mid2*sina;
					Pz[k]=-mid1*sina+mid2*cosa;

					Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);

					if(flag==0){
						continue;
					}	
					Cal_group_pxz_vp(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num);

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
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);

			//cal_group_velocity_qsv(Vag,alphy,beta,eps,deta,Ragl);
			cal_group_velocity_vp(Vag,alphy,beta,eps,deta,Ragl);
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
			//upt=10;

	        for(num=6;num<8;num++){ 
     	 		if(num==6&&T[i][j+1]>=T[i][j]&&T[i-1][j+1]>=T[i][j]){
     	 			continue;
     	 		}else if(num==7&&T[i-1][j+1]>=T[i][j]&&T[i-1][j]>=T[i][j]){
     	 			continue;
     	 		}
    			//T[i][j]=Solve_local_traveltime_modify_tti(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				T[i][j]=Solve_local_traveltime_news(T[i][j],num,Tad,d,alphy,beta,deta,eps,Ragl,downt,upt);
				/*flag=Solve_local_traveltime(num,Tad,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
				if(flag==0){continue;}  
     	
         		for(k=0;k<2;k++){ 
					if(Tc[k]<0){
						continue;
					}
					mid1=Px[k];
					mid2=Pz[k];
					Px[k]=mid1*cosa+mid2*sina;
					Pz[k]=-mid1*sina+mid2*cosa;

					Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);
					if(flag==0){ 
						continue;
					}
					Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num);
				
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
		     Tpre[i][j]=T[i][i];
		  }
		}
		goto start;
		printf("first order iteration %d: %f  %f\n",inputnum,sum,sum/(nx-4)/(nz-4));	
  	}else{
		sum=0;
	      for(j=2;j<nx-2;j++){
			  for(i=2;i<nz-2;i++){
					sum=sum+fabs(Tpre[i][j]-T[i][i]);
			  }
	      }
	      for(j=2;j<nx-2;j++){
	      	   for(i=2;i<nz-2;i++){
		     		Tpre[i][j]=T[i][i];
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
//calculate traveltime in the local sover by 2 oder 


// function: if px and pz are known, we can compute vx and vz 
void Cal_group_pxz(double Px[2],double Pz[2], double * V1, double * V2,double alphy,double beta,double delta, double eps,int k,int num0)
{ 
    double pz,px;
	double partialB,partialC;
	double B,C,sqqr,slowppx,low;
	double tp,ts;

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
	slowppx=-0.5*partialB+1.0/4.0/sqqr*(2.0*B*partialB-4.0*partialC);
	low=px*slowppx+2.0*pz*pz;
	*V1=slowppx/low;
	*V2=2*pz/low;
}
//function: limit the range of Tc using the adjacent points
int HighTest(int num,double Tad[8],double Tad1[8],double T,double Bigvalue)
{
	double Ta,Tb,Ta1,Tb1;
	switch(num){  
		case 0:{  
			Ta=Tad[7];
			Tb=Tad[0];
			Ta1=Tad1[7];
			Tb1=Tad1[0];
			break;
     	}
     	case 1:{
			Ta=Tad[7];
			Tb=Tad[6];
			Ta1=Tad1[7];
			Tb1=Tad1[6];
			break;
     	}
     	case 2:{
			Ta=Tad[5];
			Tb=Tad[6];
			Ta1=Tad1[5];
			Tb1=Tad1[6];
			break;
   		}
    	case 3:{
			Ta=Tad[5];
			Tb=Tad[4];
			Ta1=Tad1[5];
			Tb1=Tad1[4];
			break;
    	}
    	case 4:{   
			Ta=Tad[3];
			Tb=Tad[4];
			Ta1=Tad1[3];
			Tb1=Tad1[4];
			break;
    	}
    	case 5:{
			Ta=Tad[3];
			Tb=Tad[2];
			Ta1=Tad1[3];
			Tb1=Tad1[2];
			break;
    	}
    	case 6:{
			Ta=Tad[1];
			Tb=Tad[2];
			Ta1=Tad1[1];
			Tb1=Tad1[2];
			break;
    	}
    	case 7:{
			Ta=Tad[1];
			Tb=Tad[0];
			Ta1=Tad1[1];
			Tb1=Tad1[0];
			break;
    	}
    }

    if(Ta1>=Ta || Ta1>=T || Tb1>=Tb || Tb1>=T){
    	return 0;
    }
	if(Ta>=T || Tb>=T){
    	return 0;
    }
    if(Ta==Bigvalue || Tb==Bigvalue || Tb1==Bigvalue || Ta1==Bigvalue){
    	return 0;
    }

    return 1;
}
//function: the construction of the local solver for second oder 
int slownessHigh(int num,double Tad[8],double Tad1[8],double Tc[2],double d,double Px[2],double Pz[2],double alphy,
                double beta,double deta,double eps,double Rag,double downt,double upt)
{
	double A,B,C,D,E,F,G,H,I,J,K,L,M;
	double Ta,Tb,Ta1,Tb1;
	double a,b,c,f,e;
	int flag;
	int i;
	double tmp;
	double root[2];
	double sina,cosa;
	
	tmp=2.0*d;

    switch(num){  
	    case 0:{  
			Ta=Tad[7];
			Tb=Tad[0];
			Ta1=Tad1[7];
			Tb1=Tad1[0];
			break;
     	}
     	case 1:{
			Ta=Tad[7];
			Tb=Tad[6];
			Ta1=Tad1[7];
			Tb1=Tad1[6];
			break;
     	}
     	case 2:{
			Ta=Tad[5];
			Tb=Tad[6];
			Ta1=Tad1[5];
			Tb1=Tad1[6];
			break;
   		}
    	case 3:{
			Ta=Tad[5];
			Tb=Tad[4];
			Ta1=Tad1[5];
			Tb1=Tad1[4];
			break;
    	}
    	case 4:{   
			Ta=Tad[3];
			Tb=Tad[4];
			Ta1=Tad1[3];
			Tb1=Tad1[4];
			break;
    	}
    	case 5:{
			Ta=Tad[3];
			Tb=Tad[2];
			Ta1=Tad1[3];
			Tb1=Tad1[2];
			break;
    	}
    	case 6:{
			Ta=Tad[1];
			Tb=Tad[2];
			Ta1=Tad1[1];
			Tb1=Tad1[2];
			break;
    	}
    	case 7:{
			Ta=Tad[1];
			Tb=Tad[0];
			Ta1=Tad1[1];
			Tb1=Tad1[0];
			break;
    	}
 	}

	switch(num){  
	   	case 0:{
			F=0;
			G=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
			H=-3.0/tmp;
			I=-(Tb1-4*Tb)/tmp;
			break;
	   	}
	   	case 1:{
		    F=3.0/tmp;
		    G=(Tb1-4*Tb)/tmp;
		    H=0;
		    I=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 2:{
		    F=3.0/tmp;
		    G=(Tb1-4*Tb)/tmp;
		    H=0;
		    I=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 3:{
		    F=0;
		    G=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
		    H=3.0/tmp;
		    I=(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 4:{
			F=0;
			G=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			H=3.0/tmp;
			I=(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 5:{
			F=-3.0/tmp;
			G=-(Tb1-4*Tb)/tmp;
			H=0;
			I=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 6:{
			F=-3.0/tmp;
			G=-(Tb1-4*Tb)/tmp;
			H=0;
			I=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 7:{
			F=0;
			G=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			H=-3.0/tmp;
			I=-(Tb1-4*Tb)/tmp;
			break;
	   }
	}

	sina=sin(Rag);
	cosa=cos(Rag);
	J=F*cosa+H*sina;
	K=G*cosa+I*sina;
	L=H*cosa-F*sina;
	M=I*cosa-G*sina;

	A=1+2*eps;
	B=-1.0/alphy/alphy-1.0/beta/beta;
	C=-1.0/alphy/alphy-A/beta/beta;
	D=2.0*(1.0+deta+(eps-deta)*alphy*alphy/beta/beta);
	E=1.0/alphy/alphy/beta/beta;

	a=L*L*L*L+A*J*J*J*J+D*J*J*L*L;
	b=4*L*L*L*M+4*A*J*J*J*K+2*D*J*K*L*L+2*D*L*M*J*J;
	c=6*L*L*M*M+6*A*J*J*K*K+B*L*L+C*J*J+D*J*J*M*M+D*L*L*K*K+4*D*J*K*L*M;
	f=4*L*M*M*M+4*A*J*K*K*K+2*B*L*M+2*J*K*C+2*D*J*K*M*M+2*D*L*M*K*K;
	e=M*M*M*M+A*K*K*K*K+B*M*M+C*K*K+D*K*K*M*M+E;

	for(i=0;i<2;i++){
		root[i]=0;
	}
	flag=Solquar(downt,upt,a,b,c,f,e,root);

	if(flag==1){
		Tc[0]=root[0];
		Tc[1]=-1;
		Px[0]=J*root[0]+K;
		Pz[0]=L*root[0]+M;
		if(root[0]<max11(Ta,Ta1,Tb,Tb1)){
			Tc[0]=-1;
			return 0;
		}
		return 1;
	}
	if(flag==2){
		for(i=0;i<2;i++){
			Tc[i]=root[i];
			Px[i]=J*root[i]+K;
			Pz[i]=L*root[i]+M;
			if(root[i]<max11(Ta,Ta1,Tb,Tb1)){
				Tc[i]=-1;
			}
		}
	  	return 2;
	}

    return 0;
}
void slownessHigh_new(int num,double Tad[8],double Tad1[8],double Tc[2],double d,double Px[2],double Pz[2],double alphy,
                double beta,double deta,double eps,double Rag,double downt,double upt)
{
	double A,B,C,D,E,F,G,H,I,J,K,L,M;
	double Ta,Tb,Ta1,Tb1;
	int i;
	double tmp;
	double root[2];
	double sina,cosa;
	
	tmp=2.0*d;

    switch(num){  
	    case 0:{  
			Ta=Tad[7];
			Tb=Tad[0];
			Ta1=Tad1[7];
			Tb1=Tad1[0];
			break;
     	}
     	case 1:{
			Ta=Tad[7];
			Tb=Tad[6];
			Ta1=Tad1[7];
			Tb1=Tad1[6];
			break;
     	}
     	case 2:{
			Ta=Tad[5];
			Tb=Tad[6];
			Ta1=Tad1[5];
			Tb1=Tad1[6];
			break;
   		}
    	case 3:{
			Ta=Tad[5];
			Tb=Tad[4];
			Ta1=Tad1[5];
			Tb1=Tad1[4];
			break;
    	}
    	case 4:{   
			Ta=Tad[3];
			Tb=Tad[4];
			Ta1=Tad1[3];
			Tb1=Tad1[4];
			break;
    	}
    	case 5:{
			Ta=Tad[3];
			Tb=Tad[2];
			Ta1=Tad1[3];
			Tb1=Tad1[2];
			break;
    	}
    	case 6:{
			Ta=Tad[1];
			Tb=Tad[2];
			Ta1=Tad1[1];
			Tb1=Tad1[2];
			break;
    	}
    	case 7:{
			Ta=Tad[1];
			Tb=Tad[0];
			Ta1=Tad1[1];
			Tb1=Tad1[0];
			break;
    	}
 	}

	switch(num){  
	   	case 0:{
			F=0;
			G=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
			H=-3.0/tmp;
			I=-(Tb1-4*Tb)/tmp;
			break;
	   	}
	   	case 1:{
		    F=3.0/tmp;
		    G=(Tb1-4*Tb)/tmp;
		    H=0;
		    I=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 2:{
		    F=3.0/tmp;
		    G=(Tb1-4*Tb)/tmp;
		    H=0;
		    I=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 3:{
		    F=0;
		    G=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
		    H=3.0/tmp;
		    I=(Tb1-4*Tb)/tmp;
		    break;
	    }
	    case 4:{
			F=0;
			G=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			H=3.0/tmp;
			I=(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 5:{
			F=-3.0/tmp;
			G=-(Tb1-4*Tb)/tmp;
			H=0;
			I=(Ta1-4*Ta)/tmp-(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 6:{
			F=-3.0/tmp;
			G=-(Tb1-4*Tb)/tmp;
			H=0;
			I=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			break;
	    }
	    case 7:{
			F=0;
			G=-(Ta1-4*Ta)/tmp+(Tb1-4*Tb)/tmp;
			H=-3.0/tmp;
			I=-(Tb1-4*Tb)/tmp;
			break;
	   }
	}

	sina=sin(Rag);
	cosa=cos(Rag);
	//J=F*cosa-H*sina;
	//K=G*cosa-I*sina;

	//L=H*cosa-F*sina;
	//M=I*cosa-G*sina;

	double T_beg,T_end,tp;
	double tup,tdown;
	double val=1.0e-4;
	double px,pz,T_c;
	int tag=100,flag=-100;
	double fv,Tim;
	double V1,V2,md1,md2;
	double b,c,value;
	double tv,ts;

	//Tim=Tc;
	//analytical solution
	tv=1.0/alphy/alphy;
	ts=1.0/beta/beta;
	px=F*T_c+G;
	pz=H*T_c+I;
	if(abs11(H)>0){
		b=-2*(1+deta+(eps-deta)*alphy*alphy*ts)*px*px+tv+ts;
		c=((1+2*eps)*px*px-tv)*(px*px-ts);
		if((b*b-4.0*c)>0.0){
			value=b*0.5-sqrt(b*b-4.0*c)*0.5;
		}else{
			value=0;
		}
		if(value>0){
			Pz[0]=-sqrt(value);
			Pz[1]=sqrt(value);
			Px[0]=px;
			Px[1]=px;
		}else{
			Pz[0]=0;
			Pz[1]=0;
		}
		for(i=0;i<2;i++){

			Tc[i]=(Pz[i]-I)/H;
			if(Tc[i]<max11(Ta,Ta1,Tb,Tb1)){
				Tc[i]=-1;
			}
		}
		
	}
	if(abs11(F)>0){
		b=tv+(1.0+2.0*eps)*ts-2.0*(1.0+deta+(eps-deta)*alphy*alphy*ts)*pz*pz;
		b=b/(1.0+2.0*eps);
		c=pz*pz*pz*pz+tv*ts-(tv+ts)*pz*pz;
		c=c/(1.0+2.0*eps);
		if((b*b-4.0*c)>0.0){
			value=b*0.5-sqrt(b*b-4.0*c)*0.5;
		}else{
			value=0;
		}
		if(value>0){
			Px[0]=-sqrt(value);
			Px[1]=sqrt(value);
			Pz[0]=pz;Pz[1]=pz;
		}else{
			Px[0]=0;
			Px[1]=0;
		}
		for(i=0;i<2;i++){
			Tc[i]=(Px[i]-G)/F;
			if(Tc[i]<max11(Ta,Ta1,Tb,Tb1)){
				Tc[i]=-1;
			}
		}
	}//*/
	
	//iter method
	//T_beg=min(Ta,Tb);//
	/*T_beg=0;//min4(Ta,Ta1,Tb,Tb1);
	T_end=10;//upt;
	tup=T_end;
	tdown=T_beg;
	while(abs11(tup-tdown)>val){
		T_c=(tup+tdown)*0.5;
		px=F*T_c+G;
		pz=H*T_c+I;
		
		md1=px;
		md2=pz;
		px=md1*cosa+md2*sina;
		pz=-md1*sina+md2*cosa;
		Wave_surfacetest_vp(px,pz,alphy,beta,deta,eps,&tag,&fv);
		//printf("px,pz=%e,%e\n",px,pz); 
		//if(tag==1){
		//	printf("tup=%f,tdown=%f,tag,fv=%d,%e\n",tup,tdown,tag,fv);sleep(1);
		//}
		if(fv<0.0){
			tdown=T_c;
		}else{
			tup=T_c;
		}
		if(tag==1){
			Tc[0]=T_c;
			Px[0]=px;
			Pz[0]=pz;
			//break;
		}
	}
	//printf("Tc=%f\n",T_c);
	if(Tc[0]<max11(Ta,Ta1,Tb,Tb1)){
		Tc[0]=-1.0;
	}
	Tc[1]=-1.0;
	Px[1]=0;
	Pz[1]=0;
	//printf("Tc=%f\n",T_c);*/
}
// cal group velicty in fixed direction by a appro formula
void cal_group_velocity(double Vag[8],double vp,double vs,double epsilon, double delta,double Ragl)
{ 
	int i;
	double pi= 3.1415926535898,interval;
	double yita,Q,A,C,E,D,Vg;
	double len,Sn,Cs;

	interval=pi/4.0;
	for(i=0;i<8;i++){
		Sn=sin(interval*i+Ragl);
		Cs=cos(interval*i+Ragl);
		yita=(epsilon-delta)/(1.0+2.0*delta);
		C=1.0/vp/vp;
		Q=1.0+2.0*yita;
		A=1.0/vp/vp/(1.0+2.0*epsilon);
		E=A*Sn*Sn+C*Cs*Cs;
		len=E*E+4.0*(Q*Q-1.0)*A*C*Sn*Sn*Cs*Cs;
		D=sqrt(len);
		Vg=(1.0+2.0*Q)/2.0/(1.0+Q)*E+1.0/2.0/(1.0+Q)*D;
		Vg=1.0/Vg;
		Vg=sqrt(Vg);
		Vag[i]=Vg;
	}
}

//------------------------------------------------------------------------------------------------------/
void secondorder_traveltime(int nx,int nz,double dd,int sx,int sz,double**alpha,double**belta,double**delta,
			    double**epsilon,double**theta,double**T,int width,int total)
{
	int num0,flag;  
	double d,V1,V2;
	int i,j,k,l;
	double alphy, beta, deta, eps;
	double Pz[2], Px[2];
	double Tc[2],Tc2,Tad[8],Tad1[8];
	double Vag[8];
	double Ragl;
	double Bigvalue=50000,pi=3.1415926535898;
	double mid1, mid2;
	int inputnum;
	int number=0;
	double angle;
	double upt,downt;
	double Vp,VX,VZ;
	double sina,cosa;

	double **Tpre;
	double sum;


	Tpre=(double **)calloc(nx,sizeof(double *));
	for(i=0;i<nx;i++){
		Tpre[i]=(double *)calloc(nz,sizeof(double ));
	}

	d=dd;
	for(i=sx-width+1;i<sx+width;i++){
		if(sz-width+1<0){
			for(j=0;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source_vp(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}else{
			for(j=sz-width+1;j<sz+width;j++){
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];

				T[i][j]=Initialize_source_vp(alphy,beta,deta,eps,sx,sz,i,j,d,angle);
			} 	
		}
    }
	inputnum=-1;
	start: inputnum++;
  
    for(i=2;i<nx-2;i++){
  		for(j=2;j<nz-2;j++){    

	     		if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
	       			continue;
	     		}

	     		for(l=0;l<8;l++){
	     			Vag[l]=0;
	     		}
	       
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];
				sina=sin(Ragl);
				cosa=cos(Ragl);

				cal_group_velocity(Vag,alphy,beta,eps,deta,Ragl);

				Tad[0]=T[i-1][j];
				Tad[1]=T[i-1][j+1];
				Tad[2]=T[i][j+1];
				Tad[3]=T[i+1][j+1];
				Tad[4]=T[i+1][j];
				Tad[5]=T[i+1][j-1];
				Tad[6]=T[i][j-1];
				Tad[7]=T[i-1][j-1];

				Tad1[0]=T[i-2][j];
				Tad1[1]=T[i-2][j+2];
				Tad1[2]=T[i][j+2];
				Tad1[3]=T[i+2][j+2];
				Tad1[4]=T[i+2][j];
				Tad1[5]=T[i+2][j-2];
				Tad1[6]=T[i][j-2];
				Tad1[7]=T[i-2][j-2];

				for(k=0;k<7;k++){
					if(k==0){
						downt=min(Tad[k],Tad[k+1]);
					}else{
						downt=min(downt,Tad[k+1]);
					}
				}

				T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);
				//downt=-1000;
				upt=T[i][j];
				//upt=10;
		        for(num0=0;num0<2;num0++){ 
				 	flag=HighTest(num0,Tad,Tad1,T[i][j],Bigvalue);
				 	if(flag==0){
				 		continue;
				 	}
				 	/*flag=slownessHigh(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
				 	if(flag==0){
				 		continue;
				 	}*/
					slownessHigh_new(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
					for(k=0;k<2;k++){
							if(Tc[k]<0){
								continue;
							}
							Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);
							if(flag==0){
								continue;
							} 
								
							Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num0);

							mid1=V1;
							mid2=V2;
							V1=mid1*cosa-mid2*sina;
							V2=mid1*sina+mid2*cosa;
							Test_causality(V1,V2,num0,&flag);
						
							if(flag){ 
						  		T[i][j]=min(Tc[k],T[i][j]);
							}
				   }//*/
			}
		}   
	}
    for(i=nx-3;i>1;i--){
  		for(j=2;j<nz-2;j++){    
		     	if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
		     		continue;
		     	}	
		     	for(l=0;l<8;l++){
		     		Vag[l]=0;
			}
			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);
			cal_group_velocity(Vag,alphy,beta,eps,deta,Ragl);

			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];

			Tad1[0]=T[i-2][j];
			Tad1[1]=T[i-2][j+2];
			Tad1[2]=T[i][j+2];
			Tad1[3]=T[i+2][j+2];
			Tad1[4]=T[i+2][j];
			Tad1[5]=T[i+2][j-2];
			Tad1[6]=T[i][j-2];
			Tad1[7]=T[i-2][j-2];

			for(k=0;k<7;k++){
				if(k==0){
					downt=min(Tad[k],Tad[k+1]);
				}else{
					downt=min(downt,Tad[k+1]);
				}
			}
			//downt=-1000;
			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);
			upt=T[i][j];
			//upt=10;
			for(num0=2;num0<4;num0++){ 
				flag=HighTest(num0,Tad,Tad1,T[i][j],Bigvalue);
				if(flag==0){
					continue;
				}
				// solve quartic equation
		       	flag=slownessHigh(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
				if(flag==0){
					continue;
				}//*/
				
				// solve quadratic equation
				//slownessHigh_new(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
				for(k=0;k<2;k++){
				    if(Tc[k]<0){
				    	continue;
				    }

			        Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);
				    if(flag==0){
				    	continue;
				    }
				    Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num0);

				    mid1=V1;
				    mid2=V2;
				    V1=mid1*cosa-mid2*sina;
				    V2=mid1*sina+mid2*cosa;
				    Test_causality(V1,V2,num0,&flag);
				    
				    if(flag){ 
						T[i][j]=min(Tc[k],T[i][j]);
			    	}
			    }//*/
			}
		}
	}
  	for(i=nx-3;i>1;i--){
  		for(j=nz-3;j>1;j--){    
	     		if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
	       			continue;
	     		}
	     		for(l=0;l<8;l++){
	     			Vag[l]=0;
	     		}
				Ragl=theta[i][j];
				angle=-Ragl;
				alphy=alpha[i][j];
				beta=belta[i][j];
				deta=delta[i][j];
				eps=epsilon[i][j];
				sina=sin(Ragl);
				cosa=cos(Ragl);
				cal_group_velocity(Vag,alphy,beta,eps,deta,Ragl);

				Tad[0]=T[i-1][j];
				Tad[1]=T[i-1][j+1];
				Tad[2]=T[i][j+1];
				Tad[3]=T[i+1][j+1];
				Tad[4]=T[i+1][j];
				Tad[5]=T[i+1][j-1];
				Tad[6]=T[i][j-1];
				Tad[7]=T[i-1][j-1];

				Tad1[0]=T[i-2][j];
				Tad1[1]=T[i-2][j+2];
				Tad1[2]=T[i][j+2];
				Tad1[3]=T[i+2][j+2];
				Tad1[4]=T[i+2][j];
				Tad1[5]=T[i+2][j-2];
				Tad1[6]=T[i][j-2];
				Tad1[7]=T[i-2][j-2];

				for(k=0;k<7;k++){
					if(k==0){
						downt=min(Tad[k],Tad[k+1]);
					}else{
						downt=min(downt,Tad[k+1]);
					}
				}
				
				T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);

				upt=T[i][j];
				//downt=-1000;
				///upt=10;
     		    for(num0=4;num0<6;num0++){ 
					flag=HighTest(num0,Tad,Tad1,T[i][j],Bigvalue);
					if(flag==0){
						continue;
					}
					flag=slownessHigh(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
					if(flag==0){
						continue;
					}//*/
					//slownessHigh_new(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
		 			for(k=0;k<2;k++){
						if(Tc[k]<0){
							continue;
						}

						Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);
						if(flag==0){
							continue;
						}

						Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num0);

						mid1=V1;
						mid2=V2;
						V1=mid1*cosa-mid2*sina;
						V2=mid1*sina+mid2*cosa;
						Test_causality(V1,V2,num0,&flag);

						if(flag){
							T[i][j]=min(Tc[k],T[i][j]);
						}
		 		}//*/
		    }
		}
	}
 	for(i=2;i<nx-2;i++){
	 	for(j=nz-3;j>1;j--){    
			if(i>=sx-width+1 && i<sx+width && j>=sz-width+1 && j<sz+width){
	       			continue;
			}
			for(l=0;l<8;l++){
				Vag[l]=0;
			}
			Ragl=theta[i][j];
			angle=-Ragl;
			alphy=alpha[i][j];
			beta=belta[i][j];
			deta=delta[i][j];
			eps=epsilon[i][j];
			sina=sin(Ragl);
			cosa=cos(Ragl);
			cal_group_velocity(Vag,alphy,beta,eps,deta,Ragl);

			Tad[0]=T[i-1][j];
			Tad[1]=T[i-1][j+1];
			Tad[2]=T[i][j+1];
			Tad[3]=T[i+1][j+1];
			Tad[4]=T[i+1][j];
			Tad[5]=T[i+1][j-1];
			Tad[6]=T[i][j-1];
			Tad[7]=T[i-1][j-1];

			Tad1[0]=T[i-2][j];
			Tad1[1]=T[i-2][j+2];
			Tad1[2]=T[i][j+2];
			Tad1[3]=T[i+2][j+2];
			Tad1[4]=T[i+2][j];
			Tad1[5]=T[i+2][j-2];
			Tad1[6]=T[i][j-2];
			Tad1[7]=T[i-2][j-2];


			for(k=0;k<7;k++){
				if(k==0){
					downt=min(Tad[k],Tad[k+1]);
				}else{
					downt=min(downt,Tad[k+1]);
				}
			}

			T[i][j]=update_traveltime(T[i][j],Tad,Vag,d);
			upt=T[i][j];
			//upt=10;
			//downt=-1000;
			for(num0=6;num0<8;num0++){ 
				flag=HighTest(num0,Tad,Tad1,T[i][j],Bigvalue);
				if(flag==0){
					continue;
				}
				flag=slownessHigh(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
				if(flag==0){
					continue;
				}//*/
				slownessHigh_new(num0,Tad,Tad1,Tc,d,Px,Pz,alphy,beta,deta,eps,Ragl,downt,upt);
	 		    for(k=0;k<2;k++){
		    		if(Tc[k]<0){
		    			continue;
		    		}
					Wave_surfacetest(Px[k],Pz[k],alphy,beta,deta,eps,&flag);

					if(flag==0){continue;}

					Cal_group_pxz(Px,Pz,&V1,&V2,alphy,beta,deta,eps,k,num0);

					mid1=V1;
					mid2=V2;
					V1=mid1*cosa-mid2*sina;
					V2=mid1*sina+mid2*cosa;
					Test_causality(V1,V2,num0,&flag);

	  				if(flag){ 
		  				T[i][j]=min(Tc[k],T[i][j]);
	    		    }
		 		}//*/
	        }
	    }
	}
	
  	sum=0;
  	if(inputnum==0){
	      	for(j=2;j<nx-2;j++){
			  	for(i=2;i<nz-2;i++){
				 		Tpre[i][j]=T[i][i];
			  	}
	      	}
  	}else{
		sum=0;
      		for(j=2;j<nx-2;j++){
 				for(i=2;i<nz-2;i++){
          			sum=sum+fabs(Tpre[i][j]-T[i][i]);
          		}  
      		}
        	for(j=2;j<nx-2;j++){
        		for(i=2;i<nz-2;i++){
                		Tpre[i][j]=T[i][i];
               }
        	}  
		printf("Second order iteration %d: %e %f\n",inputnum,sum,sum/(nx-4.0)/(nz-4.0));
    }
	if(inputnum<total){
  		goto start;
	}

}
