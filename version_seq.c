#include <stdio.h>
#include<stdlib.h>
#include<math.h>

void matvec_opti(double *Ax, double *x, double a, double b, double c, int Nx, int Ny);
void grad_c(double*x, double*F, double a, double b, double c, int Nx, int Ny);
void save_sol(double *x, int Nx, int Ny, double Lx, double Ly);

int main(int argc, char* argv[])
{
    double *x,*F;
    double a,b,c;
    int Nx=3,Ny=3;

    a = 1.0;
    b = 8.0;
    c = -2.0;

    x = (double *)malloc((Nx*Ny)*sizeof(double));
    F = (double *)malloc((Nx*Ny)*sizeof(double));

    for(int i=0;i<Nx*Ny;i++){
        x[i] = 0.0;
    }
    F[0] = 32.0;
    F[1] = 40.0;
    F[2] = 44.0;
    F[3] = 61.0;
    F[4] = 80.0;
    F[5] = 89.0;
    F[6] = 26.0;
    F[7] = 40.0;
    F[8] = 38.0;

    
    // for(int i=0;i<Nx*Ny;i++){
    //     printf("F[%d] = %f\n",i,F[i]);
    // }

    // for(int i=0;i<Nx*Ny;i++){
    //     printf("x[%d] = %f\n",i,x[i]);
    // }

    grad_c(x,F,a,b,c,Nx,Ny);

    for(int i=0;i<Nx*Ny;i++){
        printf("x[%d] = %f\n",i,x[i]);
    }
    


    return 0;
}


void matvec_opti(double *Ax, double *x, double a, double b, double c, int Nx, int Ny){
    for(int i=1;i<Nx+1;i++){
        for(int j=1;j<Ny+1;j++){
            //Cas generaux 

            //Cas general du bloc general
            if(i!=1 && i!=Nx && j!=1 && j!=Ny){
                Ax[(j-1)*Nx+(i-1)] = b*(x[(j-2)*Nx+(i-1)] + x[j*Nx+(i-1)]) + a*(x[(j-1)*Nx+(i-2)] + x[(j-1)*Nx+i]) + c*x[(j-1)*Nx+(i-1)];
            }
            //Cas general bloc haut
            else if(i!=1 && i!=Nx && j==1){
                Ax[(j-1)*Nx+(i-1)] = b*x[j*Nx+(i-1)] + a*(x[(j-1)*Nx+(i-2)] + x[(j-1)*Nx+i]) + c*x[(j-1)*Nx+(i-1)];
            }
            //Cas general bloc bas
            else if(i!=1 && i!=Nx && j==Ny){
                Ax[(j-1)*Nx+(i-1)] = b*x[(j-2)*Nx+(i-1)] + a*(x[(j-1)*Nx+(i-2)] + x[(j-1)*Nx+i]) + c*x[(j-1)*Nx+(i-1)];
            }

            //Cas particuliers du bloc general

            //Cas haut du bloc general
            else if(i==1 && j!=1 && j!=Ny){
                Ax[(j-1)*Nx+(i-1)] = b*(x[(j-2)*Nx+(i-1)] + x[j*Nx+(i-1)]) + a*x[(j-1)*Nx+i] + c*x[(j-1)*Nx+(i-1)];
            }
            //Cas bas du bloc general
            else if(i==Nx && j!=1 && j!=Ny){
                Ax[(j-1)*Nx+(i-1)] = b*(x[(j-2)*Nx+(i-1)] + x[j*Nx+(i-1)]) + a*x[(j-1)*Nx+(i-2)] + c*x[(j-1)*Nx+(i-1)];
            }

            //Cas particuliers du bloc haut

            //Cas haut du bloc haut
            else if(i==1 && j==1){
                Ax[(j-1)*Nx+(i-1)] = b*x[j*Nx+(i-1)] + a*x[(j-1)*Nx+i] + c*x[(j-1)*Nx+(i-1)];
            }
            //Cas bas du bloc haut
            else if(i==Nx && j==1){
                Ax[(j-1)*Nx+(i-1)] = b*x[j*Nx+(i-1)] + a*x[(j-1)*Nx+(i-2)] + c*x[(j-1)*Nx+(i-1)];
            }

            //Cas particuliers du bloc bas

            //Cas haut du bloc bas
            else if(i==1 && j==Ny){
                Ax[(j-1)*Nx+(i-1)] = b*x[(j-2)*Nx+(i-1)] + a*x[(j-1)*Nx+i] + c*x[(j-1)*Nx+(i-1)];
            }
            //Cas bas du bloc bas
            else if(i==Nx && j==Ny){
                Ax[(j-1)*Nx+(i-1)] = b*x[(j-2)*Nx+(i-1)] + a*x[(j-1)*Nx+(i-2)] + c*x[(j-1)*Nx+(i-1)];
            }
            else{
                printf("ERREUR : produit matrice vecteur pour i = %d  et  j = %d \n",i,j);
            }
        }
    }
}

void grad_c(double*x, double*F, double a, double b, double c, int Nx, int Ny){
    int kmax = 1000,k=0;
    double *r_k,*r_kpun,*p,*Ax;
    double alpha,beta,denom;
    double norm_r,res_min = 1e-3;

    r_k = (double *)malloc((Nx*Ny)*sizeof(double));
    r_kpun = (double *)malloc((Nx*Ny)*sizeof(double));
    p = (double *)malloc((Nx*Ny)*sizeof(double));
    Ax = (double *)malloc((Nx*Ny)*sizeof(double));

    //Calcul du résidu initial:
    matvec_opti(Ax,x,a,b,c,Nx,Ny); //Calcul de A*x = Ax
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
        matvec_opti(Ax,p,a,b,c,Nx,Ny);
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
            return;
        }

    }

    printf("Residu : %f\n", norm_r);

    free(r_k),free(r_kpun),free(p),free(Ax);

    return;    
}

void save_sol(double *x, int Nx, int Ny,double Lx, double Ly){
    double sol_exacte,dx,dy;
    dx = Lx/(Nx+1);
    dy = Ly/(Ny+1);
    FILE *out_file = fopen("output.dat", "w");
    
    if (out_file == NULL) 
            {   
              printf("Error! Could not open file\n"); 
              exit(-1);
            } 

    for(int i=0;i<Nx*Ny;i++){
        fprintf(out_file, "%f  %f  %f %f\n", x[i],sol_exacte,i*dx,i*dy);
    }

    return;
}