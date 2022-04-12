#include <stdio.h>
#include<stdlib.h>



int main(int argc, char* argv[])
{

    //Lecture des paramètres [TODO] Faire une fonction qui lit les paramètres dans un fichier
    int Nx,Ny;
    double Lx,Ly,D,dt;
    

    Nx = 4;
    Ny = 3;
    D = 1.0;

    //Initialisation des structures de données


    //Résolution
    //------------------------//


    //Gradient Conjugué
    int kmax,k;
    double* r,p,x,Ax;
    double alpha,beta;
    //Déclarer tous les vecteurs

    //Calcul du résidu initial:
    for(int i=0; i<Nx*Ny;i++){
        //faire matvec qui sort Ax
        r[i] = b[i] - Ax[i];
        p[i] = b[i] - Ax[i];
    }

    k = 0;
    while(){
        alpha = 0
        for(int i=0; i<Nx*Ny;i++){
            alpha += r[i]*r[i];
        }
    }








        return 0;
}