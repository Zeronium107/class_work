#include <iostream>
#include <ctime>
#include <cmath>
struct Sample {
    long double rpoint;
    int xcount;
};

double weight(double, double);
Sample Metropolis(long double [], int, double);

const double d = 1.122462; //optimal spacing
const long double kb = 1.0; //boltzmann constant
float epsilon = 1.0;
const float sigma = 1.0;

using namespace std;

int main() {
    long double Tvals[4]= {0.1*epsilon/kb, 0.2*epsilon/kb, 0.3*epsilon/kb, 0.4*epsilon/kb }; //for looping through the b's
    long double bondlegth[100]={}; //bonds legth array
    int N= 500000;
    int m = 100;
    Sample sp;
    const int meas = 50000;
    int nskip=1000;

    srand(time(0));

    long* length=new long[N];
    double* e=new double[meas/nskip];
    double xaccept=0;

    for (int kti = 0; kti < 4; ++kti) 
    {
        double T = Tvals[kti];
        double xaccept=0;
        for (int i=1; i<100; i++) {
            bondlegth[i]=d; //set the length of the bonds, so the distance from bond zero bondslenght[0] = 0  to bondslenght[1] is d?
        }
        bondlegth[0] = 0;

        for (int k=0; k<N; k++)//to see if lenght change is accepted
        {
            int x1 = rand()%(m-2) + 1;
            sp = Metropolis(bondlegth, x1, T); //1st arguement is for the separation between bonds, second arguement is to choose a bond radomly
            bondlegth[x1] = sp.rpoint;                                    //third is to change that legth fourth is just the temp, fith argument is to decide wether to accept change
            xaccept=xaccept+sp.xcount;

            length[k]=0;//to total length of the system
            double sumE = 0.0;
            for (int i = 1; i < 100; ++i)
                sumE += bondlegth[i];
            length[k] = sumE;

        }
        int start = N - meas;
        int nmeas = meas/nskip;

        for (int u = 0; u < nmeas; ++u) 
        {
            int idx = start + (u+1)*nskip - 1;
            e[u] = length[idx];
            //cout << m[u] << endl;
        }
        double avg = 0;
        for (int u = 0; u < nmeas; ++u) 
            avg += e[u];
        avg /= nmeas;
        cout << T << " " << avg <<endl;

        double Raccept=xaccept/(N);
        //cout << "The acceptance rate is " << Raccept << endl;

    }
    delete[] length;
    delete[] e;

    return 0;
}
double weight(double dH, double T)
{
    return (exp(-dH/T)); //just T since the temps are in units of kb^-1
}

Sample Metropolis(long double bondlegth[], int x1, double T) 
{
    const double h=0.05; //could not get acceptance rate to 0.5
    double dx = (rand()/double(RAND_MAX))*2*h-h;// uniform in [–h,h]
    double w  =  rand()/double(RAND_MAX);
    long double old = bondlegth[x1];
    Sample spold;
    Sample spnew;
    spold.xcount = 0;
    spnew.xcount = 0;
    spold.rpoint = bondlegth[x1];
    double A = sigma/bondlegth[x1];
    
    //find U for the bond we are dealing with
    double U_old = 4*epsilon*(pow(A, 12)-pow(A, 6));

    //change bondlength
    bondlegth[x1] += dx;
    
    //find new u
    double A_new = sigma/bondlegth[x1];
    double U_new= 4*epsilon*(pow(A_new, 12)-pow(A_new, 6));
    double dE= U_new - U_old;
    spnew.rpoint = bondlegth[x1];

    if (dE <= 0) 
    {
        spnew.xcount = 1;
        return spnew;
    }
    if (weight(dE, T) >= w)
    {
        spnew.xcount=1;
        return spnew;
    }
    else
    {
        bondlegth[x1] = old;
        return spold; 
    }
    }
