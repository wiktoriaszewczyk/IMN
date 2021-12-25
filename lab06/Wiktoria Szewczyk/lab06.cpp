#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>
#include "mgmres.h"

using namespace std;

const double delta = 0.1; // dx = dy = delta

int epsilon(int epsilon1, int epsilon2, int l, int nx_a){
    int j = floor(l / (nx_a + 1));
    int i = l - j * (nx_a + 1);
    if( i <= (nx_a/2))
        return epsilon1;
    return epsilon2;
}

double p1(double x, double y){
    return 0;
}

double p2(double x, double y){
    double xmax = delta * 100;
    double ymax = delta * 100;
    double sigma = xmax / 10;
    double ro1 =  exp( -pow((x - 0.25 * xmax) / sigma, 2) - pow((y - 0.5 * ymax) / sigma, 2) );
    double ro2 = -exp( -pow((x - 0.75 * xmax) / sigma, 2) - pow((y - 0.5 * ymax) / sigma, 2) );
    return ro1 + ro2;
}

void algorytm(int nx, int ny, int epsilon1, int epsilon2, int V1, int V2, int V3, int V4, double (*p_fun)(double, double), int ifV, string name){
    const int N = (nx+1)*(ny+1);

    double* a = new double[5*N];  // niezerowe wartości el. macierzowych
    int* ja = new int[5*N];    // info. o nr. kolumn
    int* ia = new int[N+1];   // wsk do el. rozpoczynających dany wiersz
    fill_n(ia, N+1, -1);
    double* b = new double[N];
    double* V = new double[N];

    // Algorytm wypełniania macierzy rzadkiej w formacie CSR + WB Dirichleta
    ////////////////////////////////////////////////////////////////////////////
    int k = -1;

    for(int l=0; l<N; l++){
        int brzeg = 0;  // 0 - środek obszaru 1 - brzeg
        double vb = 0;  // potencjał na brzegu

        int j = floor(l / (nx + 1));
        int i = l - j * (nx + 1);
        
        int x = i * delta;
        int y = j * delta;

        if(i==0){   // lewy brzeg
            brzeg = 1;
            vb = V1;
        }
        else if(j==ny){  // górny brzeg
            brzeg = 1;
            vb = V2;
        }
        else if(i==nx){   // prawy brzeg
            brzeg = 1;
            vb = V3;
        }
        else if(j==0){   // dolny brzeg
            brzeg = 1;
            vb = V4;
        }

        if(brzeg == 1){
            b[l] = vb;
        }
        else{
            b[l] = p_fun(x,y);
        }

        ia[l] = -1;

        // lewa skrajna przekątna
        if((l - nx - 1) >= 0 && !brzeg){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsilon(epsilon1, epsilon2, l,nx) / pow(delta,2); // a_(l,l-nx-1)
            ja[k] = l - nx -1;
        }

        // poddiagonala
        if((l - 1) >= 0 && !brzeg){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsilon(epsilon1, epsilon2, l,nx) / pow(delta,2); //a_(l,l-1)
            ja[k] = l - 1;
        }

        // diagonala
        k++;
        if(ia[l] < 0)
            ia[l] = k;
        if(!brzeg){
            a[k] = -(2 * epsilon(epsilon1, epsilon2, l,nx) + epsilon(epsilon1, epsilon2, l+1,nx) + epsilon(epsilon1, epsilon2, l+nx+1,nx)) / pow(delta,2); // a_(l,l)
        }
        else{
            a[k] = 1;
        }
        ja[k] = l;

        // naddiagonala
        if(l < N && !brzeg){
            k++;
            a[k] = epsilon(epsilon1, epsilon2, l+1,nx) / pow(delta,2); // a_(l,l+1)
            ja[k] = l + 1;
        }

        if(l < (N - nx - 1) && !brzeg){
            k++;
            a[k] = epsilon(epsilon1, epsilon2, l+nx+1,nx) / pow(delta,2); //a_(l,l+nx+1)
            ja[k] = l + nx + 1;
        }
    }

    int nz_num = k + 1;
    ia[N] = nz_num;

    if(!ifV){
        string file_a_name = "data/Zadanie_3_A.txt";
        string file_b_name = "data/Zadanie_3_B.txt";
        ofstream file_a, file_b;
        file_a.open(file_a_name);
        file_b.open(file_b_name);

        for(int l=0; l<N; l++){

            int j = floor(l / (nx + 1));
            int i = l - j * (nx + 1);

            file_a<<l<<"\t"<<i<<"\t"<<j<<"\t"<<a[l]<<endl;
            file_b<<l<<"\t"<<i<<"\t"<<j<<"\t"<<b[l]<<endl;
        
        }

        file_a.close();
        file_b.close();

        
    }
    if(ifV){
        pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, 500, 500, 10e-8, 10e-8);
        name = "data/" + name;
        ofstream file;
        file.open(name);

        for(int l=0; l<N; l++){

            double y = floor(l / (nx + 1));
            double x = (l - y * (nx + 1)) * delta;
            y *= delta;

            file<<x<<"\t"<<y<<"\t"<<V[l]<<endl;
        
        }

        file.close();
    }
    cout<<name<<endl;

    delete[] a;
    delete[] ja;
    delete[] ia;
    delete[] b;
    delete[] V;
}

int main(){

    // Zadanie 3
    int nx = 4;
    int ny = 4;
    int epsilon1 = 1;
    int epsilon2 = 1;
    int V1 = 10;
    int V2 = -10;
    int V3 = 10;
    int V4 = -10;
    int ifV = 0;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p1, ifV, "");

    // Zadanie 5a 
    ifV = 1;

    string name = "5a.txt";
    nx = 50;
    ny = 50;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p1, ifV, name);

    // Zadanie 5b

    name = "5b.txt";
    nx = 100;
    ny = 100;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p1, ifV, name);

    // Zadanie 5c

    name = "5c.txt";
    nx = 200;
    ny = 200;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p1, ifV, name);

    // Zadanie 6a

    name = "6a.txt";
    nx = 100;
    ny = 100;
    V1 = 0;
    V2 = 0;
    V3 = 0;
    V4 = 0;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p2, ifV, name);

    // Zadanie 6b

    name = "6b.txt";
    epsilon2 = 2;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p2, ifV, name);

    // Zadanie 6b

    name = "6c.txt";
    epsilon2 = 10;

    algorytm(nx, ny, epsilon1, epsilon2, V1, V2, V3, V4, p2, ifV, name);

    return 0;
}


// string filename = "data/S_k"+to_string(k)+".txt";
// ofstream file;
// file.open(filename);
// file<<"\t"<<endl;
// file.close();   