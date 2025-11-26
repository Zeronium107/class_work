#include <iostream>
#include <math.h>
#include <fstream>

const long double b = 0.2;
const long double w0 = 1;
//const long double w0 = 0.5

using namespace std;

struct Dyna{
    long double RKx;
    long double RKw;
};

Dyna RungeKutta(long double x, long double w, long double t, long double dt);

long double a(long double x, long double t);

int main() {
    const long double pi = 3.14159265358979323846; 
     const int period = 300*pi;//for a simple pendulum where L = 1 and g=1
     long double dt=0.02;
     int N;
     N=period/dt;
     long double* t=new long double[N];
     long double* x=new long double [N];
     long double* w=new long double [N];
     Dyna y1;
     x[0] = 0.314;
     w[0]= 0;
     t[0]= 0;
    
    ofstream Example;
    Example.open("driven_pendulum_1.txt");
    Example << "time theta w" << endl;
    Example << t[0] << " " << x[0] << " " << w[0] << endl;
    
     for (int i = 0; i< N-1; i++)
     {
         y1 = RungeKutta(x[i], w[i], t[i], dt);
         x[i+1]=y1.RKx;
         w[i+1]=y1.RKw;
         t[i+1]=dt*(i+1);
         Example << t[i+1] << " "<< x[i+1] << " " << w[i+1] << endl;
     }
     delete [] t;
     delete [] x;
     delete [] w;
    Example.close();
}

Dyna RungeKutta(long double x1, long double v1, long double t1,
long double dt1)
{
long double g1[2]={};
long double g2[2]={};
long double g3[2]={};
long double g4[2]={};
     Dyna y2;
     long double tempx;
     long double tempv;
     long double tempt;

     g1[0]=v1;
     g1[1]=a(x1, t1);

     tempx=x1+dt1*g1[0]/2;
     tempv=v1+dt1*g1[1]/2;
     tempt = t1+dt1/2;
     g2[0]=tempv;
     g2[1]=a(tempx, tempt);//need to update the time for each RK function

     tempx=x1+dt1*g2[0]/2;
     tempv=v1+dt1*g2[1]/2; 
     tempt = t1+dt1/2;
     g3[0]=tempv;
     g3[1]=a(tempx, tempt);
    
     tempx=x1+dt1*g3[0];
     tempv=v1+dt1*g3[1];
     tempt = t1+dt1;
     g4[0]=tempv;
     g4[1]=a(tempx, tempt);
    
     y2.RKx=x1+dt1*(g1[0]+2*g2[0]+2*g3[0]+g4[0])/6;
     y2.RKw=v1+dt1*(g1[1]+2*g2[1]+2*g3[1]+g4[1])/6;
    
return y2;
}

long double a(long double x2, long double t1) //acceleration
{
     long double a1;
     a1=-sin(x2)+b*cos(w0*t1);
return a1;
}

