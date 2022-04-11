#include "fonctions.h"

double f1(double x,double y,double t){
    return 2*(y - y*y + x - x*x);
}

double g1(double x,double y,double t){
    return 0.0;
}

double h1(double x,double y,double t){
    return 0.0;
}

double f2(double x,double y,double t){
    return (sin(x)+ cos(y));
}

double g2(double x,double y,double t){
    return (sin(x)+ cos(y));
}

double h2(double x,double y,double t){
    return (sin(x)+ cos(y));
}

double f3(double x,double y,double t){
    return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(PI*t/2);
}

double g3(double x,double y,double t){
    return 0.0;
}

double h3(double x,double y,double t){
    return 1.0;
}


