#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

// dx = dy = delta
// xi = delta * i,  i=0,1,...,nx
// yi = delta * j,  j=0,1,...,ny
// psi(x,y) = psi(xi,yi) = psi_ij
// dzeta(x,y) = dzeta(xi,yi) = dzeta_ij

const double delta = 0.01;
const int ro = 1;
const int mi = 1;
const int nx = 200;
const int ny = 90;
const int i_1 = 50;
const int j_1 = 55;
const int IT_MAX = 20000;

void WB_psi(double psi[nx+1][ny+1], double Qwe, double Qwy, double y_ny, double y_j1){

    double y;

    // brzeg A (wejście)
    for(int j=j_1; j<=ny; j++){
        y = delta * j;
        psi[0][j] = (Qwe / (2 * mi)) * (pow(y,3) / 3 - (pow(y,2) / 2) * (y_j1 + y_ny) + y * y_j1 * y_ny);
    }

    // brzeg C (wyjście)
    for(int j=0; j<=ny; j++){
        y = delta * j;
        psi[nx][j] = (Qwy / (2 * mi)) * (pow(y,3) / 3 - (pow(y,2) / 2) * y_ny) + (Qwe * pow(y_j1,2) * (-y_j1 + 3 * y_ny)) / (12 * mi);
    }

    // brzeg B
    for(int i=1; i<=nx-1; i++){
        psi[i][ny] = psi[0][ny];
    }

    // breg D
    for(int i=i_1; i<=nx-1; i++){
        psi[i][0] = psi[0][j_1];
    }

    // brzeg E
    for(int j=1; j<=j_1; j++){
        psi[i_1][j] = psi[0][j_1];
    }

    // brzeg F
    for(int i=1; i<=i_1; i++){
        psi[i][j_1] = psi[0][j_1];
    }

}

void WB_dzeta(double dzeta[nx+1][ny+1], double psi[nx+1][ny+1], double Qwe, double Qwy, double y_ny, double y_j1){
    double y;

    // brzeg A (wejście)
    for(int j=j_1; j<=ny; j++){
        y = delta * j;
        dzeta[0][j] = (Qwe / (2 * mi)) * (2 * y - y_j1 - y_ny);
    }

    // brzeg C (wyjście)
    for(int j=0; j<=ny; j++){
        y = delta * j;
        dzeta[nx][j] = (Qwy / (2 * mi)) * (2 * y - y_ny);
    }

    // brzeg B
    for(int i=1; i<=nx-1; i++){
        dzeta[i][ny] = (2 / pow(delta,2)) * (psi[i][ny-1] - psi[i][ny]);
    }

    // brzeg D
    for(int i=i_1+1; i<=nx-1; i++){
        dzeta[i][0] = (2 / pow(delta,2)) * (psi[i][1] - psi[i][0]);
    }

    //brzeg E
    for(int j=1; j<=j_1-1; j++){
        dzeta[i_1][j] = (2 / pow(delta,2)) * (psi[i_1+1][j] - psi[i_1][j]);
    }

    // brzeg F
    for(int i=1; i<=i_1; i++){
        dzeta[i][j_1] = (2 / pow(delta,2)) * (psi[i][j_1+1] - psi[i][j_1]);
    }

    //wierzcholek E/F
    dzeta[i_1][j_1] = 0.5 * (dzeta[i_1-1][j_1] + dzeta[i_1][j_1-1]);

}

void relaksacja(int Qwe, string filename){
    double psi[nx+1][ny+1] = {};
    double dzeta[nx+1][ny+1] = {};
    double u[nx+1][ny+1] = {};
    double v[nx+1][ny+1] = {};

    double y_ny = delta * ny;
    double y_j1 = delta * j_1;
    double Qwy = Qwe * ((pow(y_ny,3) - pow(y_j1,3) - 3 * y_j1 * pow(y_ny,2) + 3 * pow(y_j1,2) * y_ny) / pow(y_ny,3));
    double gamma;
    int j_2 = j_1 +2;

    WB_psi(psi, Qwe, Qwy, y_ny, y_j1);

    bool isBorder;
    bool isInside;
    int omega = 0;
    for(int it=1; it<=IT_MAX; it++){
        if(it < 2000)
            omega = 0;
        else
            omega = 1;
        
        for(int i=1; i<=nx-1; i++){
            for(int j=1; j<=ny-1; j++){
                if(i > i_1 || j > j_1){
                    psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - pow(delta,2) * dzeta[i][j]);
                    dzeta[i][j] = 0.25 * (dzeta[i+1][j] + dzeta[i-1][j] + dzeta[i][j+1] + dzeta[i][j-1])
                        - omega * (ro / (16 * mi)) * ((psi[i][j+1] - psi[i][j-1]) * (dzeta[i+1][j] - dzeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j]) * (dzeta[i][j+1] - dzeta[i][j-1]));
                }
                u[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2 * delta);
                v[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2 * delta);

            }
        }

        WB_dzeta(dzeta,psi,Qwe,Qwy,y_ny,y_j1);

        gamma = 0;

        for(int i=1; i<=nx-1; i++){
            gamma += (psi[i+1][j_2] + psi[i-1][j_2] + psi[i][j_2+1] + psi[i][j_2-1] - 4 * psi[i][j_2] - pow(delta,2) * dzeta[i][j_2]);
        }
        cout<<gamma<<endl;

    }

    ofstream file;
    file.open(filename);

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            file<< delta * i << " " << delta * j << " " << psi[i][j] << " " << dzeta[i][j] << " " << u[i][j] << " " << v[i][j] <<endl;
        }
        file<<endl;
    }

    file.close(); 

}

int main(){

    relaksacja(-1000, "data/Q-1000.txt");
    relaksacja(-4000, "data/Q-4000.txt");
    relaksacja(4000, "data/Q4000.txt");

    return 0;
}


// string filename = "data/S_k"+to_string(k)+".txt";
// ofstream file;
// file.open(filename);
// file<<"\t"<<endl;
// file.close();   