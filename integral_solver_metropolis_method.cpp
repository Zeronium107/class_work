#include <iostream>
#include <ctime>
#include <math.h>

struct Sample{
long double sxpoint;
long double sypoint;
int xcount;
int ycount;
};

Sample Metropolis(double, double, double, double, double);

double weight(double, double);

double g(double, double);

const double z=1.0/4.0;
using namespace std;

int main() {
    int N= 200000;
    int n1=7;
    int nskip=15;
    double* x=new double[N];
    double* y=new double[N];
    double* G=new double[N/nskip];
    double* dx=new double[N-1];
    double* dy=new double[N-1];
    double* wi=new double[N-1];
    Sample sp;
    double Sum=0;
    double yaccept=0;
    double xaccept=0;
    srand(time(0));
    for (int i=0; i< N-1; i++)
    {
        dx[i]=rand()/double(RAND_MAX);//[0,1]], actual dx range [-h,h]
        dy[i]=rand()/double(RAND_MAX);//[0,1]], actual dx range [-h,h]
        wi[i]=rand()/double(RAND_MAX);
    }
    //create a distribution of sampling points according to W
    x[0]=rand()/double(RAND_MAX);
    y[0]=rand()/double(RAND_MAX);

    double xold=x[0];
    double yold = y[0];
    for (int i=0; i<N-1; i++)
    {
        sp=Metropolis(xold, yold, dx[i], dy[i], wi[i]);
        x[i+1]=sp.sxpoint;
        y[i+1] = sp.sypoint;
        xold=x[i+1];
        yold=y[i+1];
        xaccept=xaccept+sp.xcount;
        yaccept=yaccept+sp.ycount;
    }
    //integrating g using the sampling points
    for (int j=0; j<N/nskip; j++)
    {
        double x1=x[n1+j*nskip];//skip n1 points and use one of every nskip (15) points
        double y1=y[n1+j*nskip];
        G[j]=g(x1, y1);
        Sum=Sum+G[j];
    }
    double Ave=Sum/(N/nskip);
    double Raccept=xaccept/(N-1);
    cout << "The result is " << Ave << endl;
    cout << "The acceptance rate is " << Raccept << endl;
    return 0;
}


Sample Metropolis(double xold, double yold, double dx1, double dy1, double w1)
{
    const double h=0.35;
    double p;
    Sample spnew;
    spnew.xcount=0;
    spnew.sxpoint=xold+h*2*(dx1-0.5);
    spnew.ycount=0;
    spnew.sypoint=yold+h*2*(dy1-0.5);
    if (spnew.sxpoint<0 || spnew.sypoint <0 || spnew.sxpoint>1 || spnew.sypoint>1)
    {
        spnew.sxpoint=xold;
        spnew.sypoint=yold;
    }
    else
    {
        p=weight(spnew.sxpoint, spnew.sypoint)/weight(xold, yold);
        if (p >= w1)
        {    
            spnew.xcount=1;
            spnew.ycount=1;
        }
        else
        {
            spnew.sxpoint=xold;
            spnew.sypoint=yold;
        }
    }

    return spnew;
}

double weight(double x1, double y1)
{
    return x1*y1;
}

double g(double x2, double y2)
{
    double g_func = exp(sin(x2*y2));
    double weight_func = x2*y2;
    //cout << "this is the weight " <<weight_func << endl; //for debugging
    return z*(g_func/weight_func);
}
