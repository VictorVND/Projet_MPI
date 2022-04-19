#include <stdio.h>
#include<stdlib.h>
#include <string.h>
#include<math.h>

#define PI 3.14

void matvec(double*U,double*P, int Nx, int Ny, double a, double b, double c);
void remplissage_second_membre(double* F,double* U0, double a,double b,int Nx,int Ny,double t, double dt);
double sol_exacte(double x,double y,double t);
void grad_c(double a, double b, double c, double *F,double *U,int Nx,int Ny);
double norm(double*U, int Nx, int Ny);
double f(double x,double y,double t);
double g(double x,double y,double t);
double h(double x,double y,double t);

int main(int argc, char* argv[])
{
     //Lecture des paramètres [TODO] Faire une fonction qui lit les paramètres dans un fichier
    int Nx,Ny,m,n, nbr_itr;
    double Lx,Ly,D,dt,dx,dy,t, a,b,c;
    double *F, *U0, *P, *U, *U_exact,*D1;
    FILE *fichier;
    
    n=0;
    nbr_itr=100;

    Nx = 100;
    Ny = 100;
    D = 1.0;
    
    m=Nx*Ny;
    //Initialisation des structures de données
    
     dx=1./(Nx+1);
     dy=1./(Ny+1);
     dt=0.1;

    F=(double *)malloc(m*sizeof(double));
    D1=(double *)malloc(m*sizeof(double));
    U0=(double *)malloc(m*sizeof(double));
    U_exact=(double *)malloc(m*sizeof(double));
    U=(double *)malloc(m*sizeof(double));
    P=(double *)malloc(m*sizeof(double));
      a=-dt/(dx*dx);
      b=-dt/(dy*dy);
      c=1-2*(a+b);


     for (int i = 0; i < m; i++)
     {
         U0[i]=1.0;
         P[i]=0.0;
         F[i]=0.0;
         U[i]=1.0;
        // printf("U[%d]=%f\n",i, U[i]);

     }
     t=1.0;
     printf("a= %f et b=%f\n",a,b);

     //remplissage_second_membre(F,U0, a, b, Nx, Ny,n*dt, dt);

    //Résolution
    //------------------------//
         for (int i = 0; i < m; i++)
       {
         //printf("F[%d]=%f\n",i, F[i]);
       }
      // matvec(U,P,Nx,Ny,a,b,c);
       for (int i = 0; i < m; i++)
       {
        // printf("P[%d]=%f\n",i, P[i]);
       }
       //printf("a+b+c=%f\n", a+b+c);
      // printf("2a+b+c=%f\n", 2*a+b+c);
       //printf("a+2b+c=%f\n", a+2*b+c);
       
  
     while (n<nbr_itr)
     {  
       char* base="sol";
       FILE *fichier;
       char title[10];
        sprintf(title, "%s%d.dat",base, n);

         double err=0; 
         remplissage_second_membre(F,U0, a, b, Nx, Ny,n*dt, dt);
            for (int i = 0; i < m; i++)
       {
         //printf("F[%d]=%f\n",i, F[i]);
       }
         grad_c(a,b, c,F,U,Nx,Ny);
                  for (int i = 0; i < m; i++)
       {
        // printf("U[%d]=%f\n",i, U[i]);
       }
        // for (int k = 0; k < m; k++)
        //  {
        //      int i = k%Nx + 1;
        //      int j = k/Nx + 1;
        //      U_exact[k]=sol_exacte(dx*i,dy*j,t);
        //      D1[k]=U[k]-U_exact[k];
        //  }
           for (int i = 0; i < m; i++)
       {
         //printf("U_exact[%d]=%f et U[%d]=%f et D1[%d]=%f\n",i, U_exact[i], i, U[i], i, D1[i]);
         printf(" U[%d]=%f\n", i, U[i]);
       }
         err=norm(D1, Nx,Ny);
        // printf("U au pas de temps %d est %f et %f et uexact est %f et %f \n",n, U[0], U[11], U_exact[0], U_exact[11]);
         printf("l'erreur est de %f\n", err);
         U0=U;
           fichier=fopen(title, "w+");

     for (int k = 0; k < m; k++)
     {
        int i = k%Nx + 1 ;
        int j = k/Nx + 1;
       fprintf(fichier, "%f %f %f \n " , i*dx, j*dy, U[k]);

     }

     fclose(fichier);

         n+=1;
     }
     

    //Gradient Conjugué
   
  return 0;
}

 void matvec(double*U,double*P, int Nx, int Ny, double a, double b, double c)
      {

        for (int i = 1; i < Nx+1; i++)
        {
            for (int j = 1; j < Ny+1; j++)
            {

                if ((i==1 && j!=1) && (i==1 && j!=Ny))
              {
                P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+i])+b*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);      
              }
                else if ((i==Nx && j!=1) && (i==Nx && j!=Ny))
              {
              P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i-2)])+b*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
              }
                else if ((j==Ny && i!=1) && (j==Ny && i!=Nx))
              {
              P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+b*(U[(j-2)*Nx+(i-1)]);
              }
                else if ((j==1 && i!=1) && (j==1 && i!=Nx))
              {
              P[(j-1)*Nx+(i-1)]=c*U[(j-1)*Nx+(i-1)]+a*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+b*(U[(j)*Nx+(i-1)]);
              }
                else if (i==1 && j==1)
              {
                 
                 P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i)])+b*(U[(j)*Nx+(i-1)]);
              }
                else if (i==1 && j==Ny)
              {
                 P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i)])+b*(U[(j-2)*Nx+(i-1)]);
              }
                else if (i==Nx && j==Ny)
              {

                P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i-2)])+b*(U[(j-2)*Nx+(i-1)]);
              
              }
                else if (i==Nx && j==1)
              {
                P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i-2)])+b*(U[(j)*Nx+(i-1)]);
              }
              else
              {
                P[(j-1)*Nx+(i-1)]=c*(U[(j-1)*Nx+(i-1)])+a*(U[(j-1)*Nx+(i-2)]+U[(j-1)*Nx+(i)])+b*(U[(j-2)*Nx+(i-1)]+U[(j)*Nx+(i-1)]);
                
              }
              
              
              
            }
            
        }

      }

      void remplissage_second_membre(double* F,double* U0, double a,double b,int Nx,int Ny,double t, double dt)
     {
          double dx=1./(Nx+1);
          double dy=1./(Ny+1);
          for (int i = 1; i < Nx+1; i++)
          {
              for (int j = 1; j < Ny+1; j++)
              {
                  F[i-1+Nx*(j-1)]=U0[i-1+Nx*(j-1)]+dt*f(i*dx,j*dy,t);

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


void grad_c(double a,double b,double c,double*F,double*x,int Nx,int Ny){
    int kmax = 20000;
    int k=0;
    double *r_k,*r_kpun,*p,*Ax;
    double alpha,beta,denom,gamma;
    double norm_r,res_min;

    res_min = 1e-6;

    r_k = (double *)malloc((Nx*Ny)*sizeof(double));
    r_kpun = (double *)malloc((Nx*Ny)*sizeof(double));
    p = (double *)malloc((Nx*Ny)*sizeof(double));
    Ax = (double *)malloc((Nx*Ny)*sizeof(double));

        //Calcul du résidu initial:
    matvec(x,Ax,Nx,Ny,a,b,c);
    for(int i=0; i<Nx*Ny;i++){
        r_k[i] = F[i] - Ax[i];
        p[i] = F[i] - Ax[i];
    }
    //Calcul de la norme de r
    norm_r = 0.0;
    for(int i=0;i<Nx*Ny;i++){
        norm_r += r_k[i]*r_k[i];
    }
    norm_r = sqrt(norm_r);

    while(norm_r>res_min && k<kmax+1){
        //Calcul de alpha
        alpha = 0.0;
        denom = 0.0;
        gamma=0.0;
        matvec(p,Ax,Nx,Ny,a,b,c);
        for(int i=0;i<Nx*Ny;i++){
            alpha += r_k[i]*r_k[i];
            denom += Ax[i]*p[i];
        }
        alpha = alpha/denom;

        //Calcul des nouveaux x et r
        for(int i=0;i<Nx*Ny;i++){
            x[i] = x[i] + alpha*p[i];
            r_kpun[i] = r_k[i] - alpha*Ax[i];
        }

        //Calcul de beta
        beta = 0;
        denom = 0;
        for(int i=0;i<Nx*Ny;i++){
            beta += r_kpun[i]*r_kpun[i];
            denom += r_k[i]*r_k[i];
        }
        beta = beta/denom;

        //Clacul du nouveau p
        for(int i=0;i<Nx*Ny;i++){
            p[i] = r_kpun[i] + beta*p[i];
        }

        //r_k prend la valeur de r_kpun
        for(int i=0;i<Nx*Ny;i++){
            r_k[i] = r_kpun[i];
        }

        k = k + 1;

        norm_r = 0.0;
        for(int i=0;i<Nx*Ny;i++){
            norm_r += r_k[i]*r_k[i];
        }
        norm_r = sqrt(norm_r);

        if(k>kmax){
            printf("ERREUR : Residu non atteint. Norme du residu : %f\n",norm_r);
        }

    }

    free(r_k),free(r_kpun),free(p),free(Ax);

    return;    
    
}
double norm(double *r_k, int Nx, int Ny){
    double norm_r = 0.0;
    for(int i=0;i<Nx*Ny;i++){
        norm_r += r_k[i]*r_k[i];
    }
    norm_r = sqrt(norm_r);
    return norm_r;                                                                                                                                                                                                        
}   

double f(double x,double y,double t){
   //return 2*(y - y*y + x - x*x);
    //return sin(x) + cos(y);
    return exp(-((x-0.5)*(x-0.5)))*exp(-((y-0.5)*(y-0.5)))*cos(PI*t/2);
}

double g(double x,double y,double t){
    return 0.0;
    //return sin(x) + cos(y);
}

double h(double x,double y,double t){
   // return 0.0;
    return 1.0;
     //return sin(x) + cos(y);
}     
double sol_exacte(double x,double y,double t){
    return x*(1.-x)*y*(1.-y);
  //return sin(x) + cos(y);
}  
