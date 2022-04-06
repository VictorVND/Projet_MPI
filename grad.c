#include<stdlib.h>
#include<stdio.h>
#include"grad.h"



void matvec(double*A,double*U,double*P, int Nx, int Ny)
{

  for (int i = 0; i < Nx; i++)
  {
      for (int j = 0; j < Ny; i++)
      {
          if (i=1 && (j¡=1 || j¡=Ny))
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
        }
          if (i=Nx && (j¡=1 || j¡=Ny))
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)])+c*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
        }
          if (j=Ny && (i¡=1 || i¡=Nx))
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]);
        }
          if (j=1 && (i¡=1 || i¡=Nx))
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+c*(U[(j)*Nx+(i-1)]);
        }
          if (i=1 && j=1)
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j)*Nx+(i-1)]);
        }
           if (i=1 && j=Ny)
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i)])+c*(U[(j-2)*Nx+(i-1)]);
        }
           if (i=Nx && j=Ny)
        {
         P[(j-1)*Nx+(i-1)]=a*U[(j-1)*Nx+(i-1)]+b*(U[(j-1)*Nx+(i-2)])+c*(U[(j-2)*Nx+(i-1)]);
        }
           if (i=Nx && j=1)
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
