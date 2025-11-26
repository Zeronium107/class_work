#include <iostream>
#include <math.h>

//basketball thrown upwards using relaxation method

using namespace std;

const double t=4;
const int N=30;
const double h=t/N;
const double p=1.5; 
const double g=9.8;



struct bench{
    double height[N+1];
    double diff;
};

bench relax(double u[N+1], double f[N+1]);

int main() {
    int nmax=10000;
    double del=1e-10;
    double* y=new double[N+1];
    double* f=new double[N+1];
    
    //calculate force
    for (int i=0; i<N+1; i++)
    {
       f[i]=-g; //downward force due to gravity
    }
    
    //initial guess with a sin function
    for (int i=0; i<N+1; i++)
    {
        y[i]=0.03*sin(3.141593*h*i/t);
    }

    int n=0;
    bench B;
    do
      {
          B=relax(y, f);
          for (int i=0; i<N+1; i++)
          {
          y[i]=B.height[i];
          double time = i*h;
          cout << time << "\t" <<y[i] << endl ;
          }
           //cout << "error is " << B.diff << endl;
          n++;
       }
    while ( n< nmax && B.diff > del);
    return 0;
}

bench relax(double y[N+1], double f[N+1])
{
    bench b1;
    double diff=0;
    b1.height[0]=0;
    b1.height[N]=0;
    for (int k=1; k< N; k++)
      {
          double uk_old=y[k];
          double uk_new=0.5*y[k-1]+0.5*y[k+1]-0.5*h*h*f[k];
          y[k]=(1-p)*uk_old+p*uk_new;
          b1.height[k]=y[k];
          diff+=(y[k]-uk_old)*(y[k]-uk_old);
       };
    b1.diff=diff;
    return b1;
}