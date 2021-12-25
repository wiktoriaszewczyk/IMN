#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

const double delta = 0.2;
const int nx = 128;
const int ny = 128;
const double xmax = delta*nx;
const double ymax = delta*ny;
const double TOL = pow(10.0,-8);

double VB1(double y){
    return sin(M_PI * y / ymax);
}

double VB2(double x){
    return -sin(2 * M_PI * x / xmax);
}

double VB3(double y){
    return sin(M_PI * y / ymax);
}

double VB4(double x){
    return sin(2 * M_PI * x / xmax);
}

void relaksacja(double V[nx+1][ny+1], int k, int& it){
    string filename_S = "data/S_k"+to_string(k)+".txt";
    ofstream file_S;
    file_S.open(filename_S);

    double S_old = 0;
    double S_new = 0;

    while(1){
        for(int i=k; i<=(nx-k); i+=k){
            for(int j=k; j<=(ny-k); j+=k){
                V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
            }
        }

        S_old = S_new;
        S_new = 0.0;
        for(int i=0; i<=(nx-k); i+=k){
            for(int j=0; j<=(ny-k); j+=k){
                S_new += (pow((k*delta),2)/2.0) * (pow(( ((V[i+k][j] - V[i][j]) / (2 * k * delta)) + ((V[i+k][j+k] - V[i][j+k]) / (2 * k * delta)) ),2) + pow(( ((V[i][j+k] - V[i][j]) / (2 * k * delta)) + ((V[i+k][j+k] - V[i+k][j]) / (2 * k * delta)) ),2));
            }
        }
        if(S_old != 0){
            if( abs((S_new - S_old) / S_old) < TOL)
                break;
        }
        file_S<<it<<"\t"<<S_new<<endl;
        it ++;
    }
    file_S.close();

    if(k!=1){
        for(int i=0; i<=(nx-k); i+=k){
            for(int j=0; j<=(ny-k); j+=k){
                V[i+k/2][j+k/2] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                if(i != nx-k)
                    V[i+k][j+k/2] = 0.5 * (V[i+k][j] + V[i+k][j+k]);
                if(j != ny-k)
                    V[i+k/2][j+k] = 0.5 * (V[i][j+k] + V[i+k][j+k]);
                if(j != 0)
                    V[i+k/2][j] = 0.5 * (V[i][j] + V[i+k][j]);
                if(i != 0)
                    V[i][j+k/2] = 0.5 * (V[i][j] + V[i][j+k]);
            }
        }
    }

    string filename_V = "data/V_k"+to_string(k)+".txt";
    ofstream file_V;
    file_V.open(filename_V);

    for(int i=0; i<=nx; i+=k){
        for(int j=0; j<=ny; j+=k){
            file_V<<delta*i<<"\t"<<delta*j<<"\t"<<V[i][j]<<endl;
        }
    }

    file_V.close();

}

int main(){
    int k[5] = {16, 8, 4, 2, 1};
    // Wszystkie obliczenia na jednej tablicy potencjału
    double V[nx+1][ny+1] = {0.0};

    // Warunki brzegowe (k=1)
    for(int j=0; j<ny+1; j++){
        V[0][j] = VB1(delta * j);
        V[nx][j] = VB3(delta * j);
    }

    for(int i=0; i<nx+1; i++){
        V[i][ny] = VB2(delta * i);
        V[i][0] = VB4(delta * i);
    }

    //relaksacja(V,k[0]);
    int it = 0;
    for(int i=0; i<5; i++)
        relaksacja(V,k[i],it);


    return 0;
}