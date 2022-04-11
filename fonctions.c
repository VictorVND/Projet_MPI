

double fct_f(double x,double y,double t){
    return 2*(y - y*y + x - x*x);
}

double fct_g(double x,double y,double t){
    return 0.0;
}

double fct_h(double x,double y,double t){
    return 0.0;
}

void read_data(int*Nx, int*Ny,double*Lx,double*Ly,double*D,double*dt){
    
    FILE *fptr;

    fptr = fopen("donnees.dat", "r");

    Nx = getw(fptr);
    Ny = getw(fptr);
    Lx = getw(fptr);
    Ly = getw(fptr);
    D = getw(fptr);
    dt = getw(fptr);

    fclose(fptr);

    return;
}


