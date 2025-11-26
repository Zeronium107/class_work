#include <iostream>
#include <cmath>
#include <fstream>

int target_nodes = 0;

using namespace std;

    const double e = 1.6022e-19;
    const double Ang= 1e-10;
    const double hbar=1.0546e-34;
    const double m = 9.1094e-31;
    const long double a_nbr = 1e-11;


struct Dyna{
    double RKu;
    double RKv;
};

double secant(int, double, double, double);

double f(double);

Dyna RungeKutta(double, double, double, double, double);

double a(double, double, double);

int main() {
    double del = 1e-4;
    int n_iter = 30;
    
    //ground state
    target_nodes = 0;
    double E0_ground = 200;
    double d_E_ground = 1;
    double root_ground = secant(n_iter, del*e, E0_ground*e, d_E_ground*e);
    cout << "Ground state energy: " << root_ground/e << " eV" << endl;
    //first excited state
    target_nodes = 1;
    double E0_first = 400;
    double d_E_first = 1;
    double root_first = secant(n_iter, del*e, E0_first*e, d_E_first*e);
    cout << "First excited state energy: " << root_first/e << " eV" << endl;
    //second excited state
    target_nodes = 2;
    double E0_second = 900;
    double d_E_second = 1;
    double root_second = secant(n_iter, del*e, E0_second*e, d_E_second*e);
    cout << "Second excited state energy: " << root_second/e << " eV" << endl;
    
    return 0;
}


double secant(int n, double del, double A, double dA){
int k=0;
double A1=A+dA;
while (abs(dA)> del && (k<n))
{
    double d=f(A1)-f(A);
    double A2=A1-f(A1)*(A1-A)/d;
    A = A1;
    A1 = A2;
    dA = A1-A;
    k++;
    cout << "dE = " << dA/e << " eV" << endl;
}
return A;
}

double f(double E){
    double L = 20 * a_nbr;
    int N = 2000;
    double dx = L / N;
    double* u = new double[N+1];
    double* v = new double[N+1];
    double x0 = -10 * a_nbr;
    u[0] = 0.0;
    v[0] = 1E-5;
    string filename = "HW_wavefunction_" + to_string(target_nodes) + ".txt";
    ofstream wf(filename);    
    wf << x0 << "\t" << u[0] << endl;
    Dyna y2;
    for (int i = 0; i < N; i++) {
        double x = x0 + i * dx;
        y2 = RungeKutta(u[i], v[i], x, dx, E);
        u[i+1] = y2.RKu;
        v[i+1] = y2.RKv;
        double tempx = x0+(i+1)*dx;
        wf << tempx << "\t" << u[i+1] << endl;
    }
    wf.close();
    double result = u[N];
    delete[] u;
    delete[] v;
    return result;
}


Dyna RungeKutta(double u, double v, double x, double dx, double E)
{
    double c1[2]={};
    double c2[2]={};
    double c3[2]={};
    double c4[2]={};
    Dyna y3;
    double tempu;
    double tempv;

    c1[0]=v;
    c1[1]=a(u, E, x);

    tempu=u+dx*c1[0]/2;
    tempv=v+dx*c1[1]/2;
    c2[0]=tempv;
    c2[1]=a(tempu, E, x+dx/2);

    tempu=u+dx*c2[0]/2;
    tempv=v+dx*c2[1]/2;
    c3[0]=tempv;
    c3[1]=a(tempu, E, x+dx/2);

    tempu=u+dx*c3[0];
    tempv=v+dx*c3[1];

    c4[0]=tempv;
    c4[1]=a(tempu, E, x+dx);

    y3.RKu=u+dx*(c1[0]+2*c2[0]+2*c3[0]+c4[0])/6;
    y3.RKv=v+dx*(c1[1]+2*c2[1]+2*c3[1]+c4[1])/6;
    return y3;
}

double a(double psi1, double E, double x1)
{
    double V0 = 50*e; //ev
    double x2 = x1 * x1;
    double a2 = a_nbr * a_nbr;
    double V = V0 * (x2 / a2) * (x2 / a2);  // V = V0 * x^4 / a^4
    return 2*m/hbar/hbar*(V-E)*psi1;
}