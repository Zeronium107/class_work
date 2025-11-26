#include <iostream>
#include <math.h>

//a person sitting on the bench with a dog and we find the bending of the bench

using namespace std;

const double L=3;
const int N=30;
const double h=L/N;
const double Y=1e9;
const double I=0.03*0.03*0.03*0.2/3;
const double p=1.5; 
const double rho=3.0;
const double g=9.8;
const double x0=0.25;
const double x1=0.15; //for dog


struct bench{
    double curve[N+1];
    double diff;
};

bench relax(double u[N+1], double f[N+1]);

int main() {
    int nmax=10000;
    double del=1e-10;
    double* u=new double[N+1];
    double* f=new double[N+1];
    
    //calculate force
    for (int i=0; i<N+1; i++)
    {
       double x=i*h;
    if (abs(x-L/2)<=x0)
        f[i]=-200*(exp(-(x-L/2)*(x-L/2)/x0/x0)-exp(-1))-rho*g;
    else if(abs(x-L/4)<=x1)
        f[i]=-100*(exp(-(x-L/4)*(x-L/4)/x1/x1)-exp(-1))-rho*g; //the donward strain due to the dog 
    else
         f[i]=-rho*g;
    }
    
    //initial guess with a sin function
    for (int i=0; i<N+1; i++)
    {
        u[i]=0.03*sin(3.141593*h*i/L);
    }

    int n=0;
    bench B;
    do
      {
          B=relax(u, f);
          for (int i=0; i<N+1; i++)
          {
          u[i]=B.curve[i];
          cout << u[i]  << endl ;
          }
           //cout << "error is " << B.diff << endl;
          n++;
       }
    while ( n< nmax && B.diff > del);
    return 0;
}

bench relax(double u[N+1], double f[N+1])
{
    bench b1;
    double diff=0;
    b1.curve[0]=0;
    b1.curve[N]=0;
    for (int k=1; k< N; k++)
      {
          double uk_old=u[k];
          double uk_new=0.5*u[k-1]+0.5*u[k+1]-0.5*h*h/Y/I*f[k];
          u[k]=(1-p)*uk_old+p*uk_new;
          b1.curve[k]=u[k];
          diff+=(u[k]-uk_old)*(u[k]-uk_old);
       };
    b1.diff=diff;
    return b1;
}