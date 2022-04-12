#include <stdio.h>
#include<stdlib.h>
#include<grad.h>

void matvec(double*A,double*U,double*P, int Nx, int Ny, double a, double b, double c);
void remplissage_second_membre(double* F,double* U0, double a,double b,int Nx,int Ny,double t, double dt);
double f(double x,double y,double t);
double g(double x,double y,double t);
double h(double x,double y,double t);

int main(int argc, char* argv[])
{
     
  return 0;
}

 void matvec(double*A,double*U,double*P, int Nx, int Ny, double a, double b, double c)
      {

        for (int i = 0; i < Nx+1; i++)
        {
            for (int j = 0; j < Ny+1; i++)
            {
                if (i==1 && (j!=1 || j!=Ny))
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
              }
                if (i==Nx && (j!=1 || j!=Ny))
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)])+c*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
              }
                if (j==Ny && (i!=1 || i!=Nx))
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]);
              }
                if (j==1 && (i!=1 || i!=Nx))
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+c*(U[(j)*Nx+(i-1)]);
              }
                if (i==1 && j==1)
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j)*Nx+(i-1)]);
              }
                if (i==1 && j==Ny)
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]);
              }
                if (i==Nx && j==Ny)
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)])+c*(U[(j-2)*Nx+(i-1)]);
              }
                if (i==Nx && j==1)
              {
              P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)])+c*(U[(j)*Nx+(i-1)]);
              }
              else
              {
                P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
              }
              
              
              
            }
            
        }

      }

      void remplissage_second_membre(double* F,double* U0, double a,double b,int Nx,int Ny,double t, double dt)
     {
          double dx=1/(Nx+1);
          double dy=1/(Ny+1);
          for (int i = 0; i < Nx+1; i++)
          {
              for (int j = 0; j < Ny+1; j++)
              {
                  F[i-1+Nx*(j-1)]=U0[i-1+Nx*(j-1)]+dt*f((i+1)*dx,(j+1)*dy,t);

                    if (i==1)
                    {
                    F[i-1+Nx*(j-1)]-=a*h((i-1)*dx,j*dy,t);
                    }
                    if (i==Nx)
                    {
                    F[i-1+Nx*(j-1)]-=a*h((i+1)*dx,j*dy,t);
                    }
                    if (j==Ny)
                    {
                    F[i-1+Nx*(j-1)]-=b*g(i*dx,(j+1)*dy,t);
                    }
                    if (j==1)
                    {
                    F[i-1+Nx*(j-1)]-=b*g(i*dx,(j-1)*dy,t);
                    }
              }
              
          }
          
        
      }
double f(double x,double y,double t){
    return 2*(y - y*y + x - x*x);
}

double g(double x,double y,double t){
    return 0.0;
}

double h(double x,double y,double t){
    return 0.0;
}     
