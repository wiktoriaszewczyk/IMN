#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;


const int nx = 400;
const int ny = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 0.1;   // 10 * delta; - rozmycie
const double xA = 0.45; // polozenie srodka gestosci
const double yA = 0.45;
const int IT_MAX = 11000;

// t = tn = dt * n  n = 0,1,..
// x = delta * i    i = 0,1,...,nx
// y = delta * j    j = 0,1,...,ny

void getPsiFromFile(double psi[nx+1][ny+1]){
    ifstream filePsi("data/psi.dat");
    
    int psiI, psiJ;
    double psiP;

    if(filePsi.is_open()){
        while(filePsi>>psiI>>psiJ>>psiP){
            // cout<<psiI<<"\t"<<psiJ<<"\t"<<psiP<<endl;
            psi[psiI][psiJ] = psiP;
        }
    }
    else{
        cout<<"Nie udalo sie otworzyc pliku! (psi.dat)";
    }
}

void setV(double vx[nx+1][ny+1], double vy[nx+1][ny+1], double psi[nx+1][ny+1]){
    for(int i=1; i<=nx-1; i++){
        for(int j=1; j<=ny-1; j++){
            vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2 * delta);
            vy[i][j] = - (psi[i+1][j] - psi[i-1][j]) / (2 * delta);
        }
    }

    // zastawka
    for(int i=i_1; i<=i_2; i++){
        for(int j=0; j<=j_1; j++){
            vx[i][j] = 0;
            vy[i][j] = 0;
        }
    }

    // dolny i gorny brzeg
    for(int i=1; i<=nx-1; i++){
        vx[i][0]  = 0; 
        vy[i][ny] = 0;
    }

    //lewy i prawy brzeg
    for(int j=0; j<=ny; j++){
        vx[0][j]  = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    }
    
}

void setU0(double u[nx+1][ny+1]){
    double x, y;
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            x = i * delta;
            y = j * delta;
            u[i][j] = (1 / (2 * M_PI * pow(sigma,2))) * exp(-( (pow((x - xA),2) + pow((y - yA),2)) / (2 * pow(sigma,2)) ));
        }
    }
}

void AD(double vx[nx+1][ny+1], double vy[nx+1][ny+1], double psi[nx+1][ny+1], double dt, double D){
    string filename = "data/u_D="+to_string(D)+"_";
    ofstream fileCXsr, fileu1, fileu2, fileu3, fileu4, fileu5;
    fileCXsr.open("data/t_C_Xsr_D="+to_string(D)+".dat");
    fileu1.open(filename+"1.dat");
    fileu2.open(filename+"2.dat");
    fileu3.open(filename+"3.dat");
    fileu4.open(filename+"4.dat");
    fileu5.open(filename+"5.dat");

    double u0[nx+1][ny+1];
    double u1[nx+1][ny+1];

    setU0(u0);

    int nrMap = 1;

    for(int it=1; it<=IT_MAX; it++){
        // u1 = u0
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                u1[i][j] = u0[i][j];
            }
        }

        // 
        for(int k=1; k<=20; k++){
            for(int i=0; i<=nx; i++){
                for(int j=1; j<=ny-1; j++){
                    if(j<=j_1 && i<=i_2 && i>=i_1){
                        continue;
                    }
                    else if(i == 0){
                        // periodyczne WB -lewy sasiad (0,j) to (nx,j)
                        // zamiana i-1 na nx we wzorze 9
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[i + 1][j] - u0[nx][j]) / (2.0 * delta)) + (u1[i + 1][j] - u1[nx][j]) / (2.0 * delta)) 
                                - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) 
                                + dt / 2.0 * D * ((u0[i + 1][j] + u0[nx][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) 
                                                    + (u1[i + 1][j] + u1[nx][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    }
                    else if(i == nx){
                        // periodyczne WB - prawy sasiad (nx,j) to (0,j)
                        // zamiana i+1 na 0 we wzorze 9
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[0][j] - u0[i - 1][j]) / (2.0 * delta)) + (u1[0][j] - u1[i - 1][j]) / (2.0 * delta)) 
                                - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) 
                                + dt / 2.0 * D * ((u0[0][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) 
                                                    + (u1[0][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    }
                    else{
                        // wzor 9
                        u1[i][j] = (1.0 / (1 + (2 * D * dt / pow(delta, 2)))) * (u0[i][j] - (dt / 2.0) * vx[i][j] * (((u0[i + 1][j] - u0[i - 1][j]) / (2.0 * delta)) + (u1[i + 1][j] - u1[i - 1][j]) / (2.0 * delta)) 
                                - (dt / 2) * vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2.0 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * delta)) 
                                + dt / 2.0 * D * ((u0[i + 1][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / pow(delta, 2) 
                                + (u1[i + 1][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1]) / pow(delta, 2)));
                    }
                }
            }
        }

        // u0 = u1
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                u0[i][j] = u1[i][j];
            }
        }

        // c xsr
        double c = 0;
        double xsr = 0;

        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                c += u0[i][j];  
                xsr += (i * delta) * u0[i][j];
            }
        }
        c *= pow(delta,2);
        xsr *= pow(delta,2);

        fileCXsr<<(dt*it)<<" "<<c<<" "<<xsr<<endl;

        if( (it) == (nrMap * (IT_MAX / 5))){
            for(int i=0; i<=nx; i++){
                for(int j=0; j<=ny; j++){
                    if(nrMap == 1){
                        fileu1<<i<<" "<<j<<" "<<u0[i][j]<<endl;
                    }
                    else if(nrMap == 2){
                        fileu2<<i<<" "<<j<<" "<<u0[i][j]<<endl;
                    }
                    else if(nrMap == 3){
                        fileu3<<i<<" "<<j<<" "<<u0[i][j]<<endl;
                    }
                    else if(nrMap == 4){
                        fileu4<<i<<" "<<j<<" "<<u0[i][j]<<endl;
                    }
                    else if(nrMap == 5){
                        fileu5<<i<<" "<<j<<" "<<u0[i][j]<<endl;
                    }
                }

                if(nrMap == 1){
                    fileu1<<endl;
                }
                else if(nrMap == 2){
                    fileu2<<endl;
                }
                else if(nrMap == 3){
                    fileu3<<endl;
                }
                else if(nrMap == 4){
                    fileu4<<endl;
                }
                else if(nrMap == 5){
                    fileu5<<endl;
                }
            }
            nrMap++;
        }

    }


    fileCXsr.close();
    fileu1.close();
    fileu2.close(); 
    fileu3.close(); 
    fileu4.close(); 
    fileu5.close(); 
}

int main(){

    double psi[nx+1][ny+1];
    double vx[nx+1][ny+1];
    double vy[nx+1][ny+1];

    getPsiFromFile(psi);
    setV(vx,vy,psi);

    // wyznaczanie dt (kroku czasowego) 
    double vmod, vmax = 0;
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            vmod = sqrt(pow(vx[i][j],2) + (pow(vy[i][j],2)));
            if(vmod > vmax) 
                vmax = vmod;
        }
    }
    double dt = delta / (4 * vmax);
    double D = 0;

    AD(vx, vy, psi, dt, D);
    D = 0.1;
    AD(vx, vy, psi, dt, D);


    return 0;
}
