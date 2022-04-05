#include <stdio.h>
#include<stdlib.h>
#include"matrice.h"
#include"fonctions.h"

    
    void remplissage_matrice(double*A,double a, double,b, double c, int Nx,int Ny)
     {
        
      }
     void remplissage_second_membre(double* F,double* U0, double a,double b,int Nx,int Ny,double dt)
     {
          for (int i = 0; i < Nx; i++)
          {
              for (int j = 0; j < Ny; j++)
              {
                  F[i-1+Nx*(j-1)]=U0[i-1+Nx*(j-1)]+dt*
              }
              
          }
          
        
      }
       
