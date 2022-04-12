#include <stdio.h>
#include<stdlib.h>
#include<math.h>

void matvec(double *A,double *x,double *Ax,int Nx,int Ny);
void grad_c(double a,double b,double c,double*F,double*x,int Nx,int Ny);
int main(int argc, char* argv[])
{



    return 0;
}



void matvec(double *A,double *x,double *Ax,int Nx,int Ny){
    int n = Nx*Ny;
    for(int i=0;i<n;i++){
        Ax[i] = 0.0;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            Ax[i] += A[j*n+i]*x[j];
        }
    }
    return;
}

void grad_c(double a,double b,double c,double*F,double*x,int Nx,int Ny){
    int kmax = 1000,k=0;
    double *r_k,*r_kpun,*p,*Ax;
    double alpha,beta,denom;
    double norm_r,res_min = 1e-3;

    r_k = (double *)malloc((Nx*Ny)*sizeof(double));
    r_kpun = (double *)malloc((Nx*Ny)*sizeof(double));
    p = (double *)malloc((Nx*Ny)*sizeof(double));
    Ax = (double *)malloc((Nx*Ny)*sizeof(double));

        //Calcul du résidu initial:
    matvec(A,x,Ax,Nx,Ny);
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

    while(norm_r>res_min && k<kmax){
        //Calcule de alpha
        alpha = 0.0;
        denom = 0.0;
        matvec(A,p,Ax,Nx,Ny);
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

        if(k==kmax-1){
            printf("ERREUR : Residu non atteint. Norme du residu : %f\n",norm_r);
        }

    }

    free(r_k),free(r_kpun),free(p),free(Ax);

    return;    
    
}