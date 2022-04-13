#include <stdio.h>
#include<stdlib.h>
#include<math.h>

void matvec(double *A,double *x,double *Ax,int Nx,int Ny);
void grad_c(double *A,double*F,double*x,int Nx,int Ny);
int main(int argc, char* argv[])
{
    int Nx = 3;
    int Ny = 1;
    double *A,*x,*b;
    A = (double *)malloc((Nx*Nx)*sizeof(double));
    x = (double *)malloc((Nx)*sizeof(double));
    b = (double *)malloc((Nx)*sizeof(double));

    for(int i=1;i<Nx+1;i++){
        x[i-1] = 1.0;
        // for(int j=1;j<Nx+1;j++){
        //     if(j==i){
        //         A[(j-1)*Nx+(i-1)] = -2.0;
        //     }
        //     else if(j==i+1 || j==i-1){
        //         A[(j-1)*Nx+(i-1)] = 1.0;
        //     }
        //     else{
        //         A[(j-1)*Nx+(i-1)] = 0.0;
        //     }
        // }
    }

    A[0] = 25.0, A[1] = 15.0, A[2] = -5.0;
    A[3] = 15.0, A[4] = 18.0, A[5] = 0;
    A[6] = -5.0, A[7] = 0, A[8] = 11.0;

    b[0] = -60.0;
    b[1] = -9;
    b[2] = 48.0;

    // x[0] = 20.0;
    // x[1] = 1.0;
    // x[2] = 6.0;
    // for(int i = 0;i<Nx;i++){
    //     printf("%f  %f  %f\n", A[i*Nx], A[i*Nx+1], A[i*Nx+2]);
    // }

    //matvec(A,x,b,Nx,Ny);

    grad_c(A,b,x,Nx,Ny);
    printf("%f  %f  %f\n", x[0], x[1], x[2]);
    free(A),free(x);
    


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

void grad_c(double *A,double*F,double*x,int Nx,int Ny){
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