#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
const long double G = 6.6738e-11;
const long double M = 1.9891e30;//sun
const long double Me = 5.972e24;//Earth
const long double steps = 30000;
const long double Mp = 1.89819e27;//jupiter

int main() {
    //const long double pi = 3.14159265358979323846;
    long double  dt=14400;
    //long double  dt=pi/m;
    int N=steps+1;
    long double* t=new long double[N];
    //earth
    long double* ye=new long double [N];
    long double* vye=new long double [N];
    long double* xe=new long double [N];
    long double* vxe=new long double [N];
    long double* Re=new long double[N];
    //jupiter
    long double* yp=new long double [N];
    long double* vyp=new long double [N];
    long double* xp=new long double [N];
    long double* vxp=new long double [N];
    long double* Rp=new long double[N];
    //earth
    t[0]= 0;
    ye[0]= 0;
    vye[0]= 3.0287e4;
    xe[0] = 1.4710e11;
    vxe[0] = 0;
    Re[0] =sqrt((xe[0]*xe[0])+(ye[0]*ye[0]));
    //jupiter
    yp[0]= 0;
    vyp[0]= 1.372e4;
    xp[0] = 7.40522e11;
    vxp[0] = 0;
    Rp[0] =sqrt((xp[0]*xp[0])+(yp[0]*yp[0]));
    ofstream orbit;
    orbit.open("orbit_earth_jupiter.txt");
    orbit << "     Earth            Jupiter" << endl;
    orbit << "time y x vx vy Radius y x vx vy Radius "<< endl;
    orbit << t[0] << " " << ye[0] << "  " << xe[0] << "  " <<  vxe[0] <<
    "  " << vye[0] << " " << Re[0] << " " << yp[0] << "  " << xp[0] << "  " <<  vxp[0] << "  " << vyp[0] << " " << Rp[0]<<endl;
    for(int i = 1; i<=steps; i++)
    {
        t[i]=dt*i;
        long double Rej = sqrt((xe[i-1] - xp[i-1]) * (xe[i-1] - xp[i-1]) + 
        (ye[i-1] - yp[i-1]) * (ye[i-1] - yp[i-1])); //earth jupiter distance
        //earth jupiter interactions
        long double Fej = G * Mp * Me / (Rej * Rej);
        long double Fejx = Fej * ((xp[i-1] - xe[i-1]) / Rej);
        long double Fejy = Fej * ((yp[i-1] - ye[i-1]) / Rej);
        
        //Earth
        //Euler method
        Re[i-1] = sqrt(xe[i-1] * xe[i-1] + ye[i-1] * ye[i-1]);
        ye[i] = ye[i-1] + vye[i-1]*dt;
        vye[i] = vye[i-1]-(G*M*ye[i-1])/(Re[i-1]*Re[i-1]*Re[i-1])*dt-Fejy/Me*dt;
        xe[i] = xe[i-1] + vxe[i-1]*dt;
        vxe[i] = vxe[i-1] - (G*M*xe[i-1])/(Re[i-1]*Re[i-1]*Re[i-1])*dt-Fejx/Me*dt;
        //Jupiter
        //Euler method
        Rp[i-1] = sqrt(xp[i-1] * xp[i-1] + yp[i-1] * yp[i-1]);
        yp[i] = yp[i-1] + vyp[i-1]*dt;
        vyp[i] = vyp[i-1]-(G*M*yp[i-1])/(Rp[i-1]*Rp[i-1]*Rp[i-1])*dt-Fejy/Mp*dt;
        xp[i] = xp[i-1] + vxp[i-1]*dt;
        vxp[i] = vxp[i-1] - (G*M*xp[i-1])/(Rp[i-1]*Rp[i-1]*Rp[i-1])*dt-Fejx/Mp*dt;
    

        //corrector method
        //update earth jupiter interaction
        long double Rej_up = sqrt((xe[i]-xp[i])*(xe[i]-xp[i])+(ye[i]-yp[i])*(ye[i]-yp[i]));
        long double Fej_up = G*Mp*Me/(Rej_up*Rej_up);
        long double Fejx_up = Fej_up*((xp[i]-xe[i])/Rej_up);
        long double Fejy_up = Fej_up*((yp[i]-ye[i])/Rej_up);
        
        //earth
        ye[i] = ye[i-1] + (vye[i-1]+vye[i])*dt/2;
        vye[i] = vye[i-1] - ((G*M*(ye[i-1]+ye[i])) / (Re[i-1]*Re[i-1]*Re[i-1])+(Fejy+Fejy_up)/Me)*dt/2;
        xe[i] = xe[i-1] + (vxe[i-1]+vxe[i])*dt/2;
        vxe[i] = vxe[i-1] - ((G*M*(xe[i-1]+xe[i])) / (Re[i-1]*Re[i-1]*Re[i-1])+ (Fejx+Fejx_up) / Me) * dt / 2;
        Re[i] =sqrt((xe[i]*xe[i])+(ye[i]*ye[i]));

        //jupiter
        yp[i] = yp[i-1] + (vyp[i-1]+vyp[i])*dt/2;
        vyp[i] = vyp[i-1] - ((G*M*(yp[i-1]+yp[i])) / (Rp[i-1]*Rp[i-1]*Rp[i-1])+(Fejy+Fejy_up)/Mp)*dt/2;
        xp[i] = xp[i-1] + (vxp[i-1]+vxp[i])*dt/2;
        vxp[i] = vxp[i-1] - ((G*M*(xp[i-1]+xp[i])) / (Rp[i-1]*Rp[i-1]*Rp[i-1])+ (Fejx+Fejx_up) / Mp) * dt / 2;
        Rp[i] =sqrt((xp[i]*xp[i])+(yp[i]*yp[i]));
        
        orbit << t[i] << " " << ye[i] << "  " << xe[i] << "  " <<  vxe[i] << "  " << vye[i] << " " << Re[i] << " " << yp[i] << "  " << xp[i] << "  " <<  vxp[i] << "  " << vyp[i] << " " << Rp[i]<<endl;
    }
    orbit.close();
    delete [] t;
    delete [] xe;
    delete [] ye;
    delete [] vxe;
    delete [] vye;
    delete [] xp;
    delete [] yp;
    delete [] vxp;
    delete [] vyp;
    return 0;
}