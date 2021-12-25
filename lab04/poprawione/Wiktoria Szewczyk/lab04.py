import math
import matplotlib.pyplot as plt
import numpy as np

# Dane
epsilon = 1
delta = 0.1
nx = 150
ny = 100
V1 = 10     # WB na dole
V2 = 0      # WB na górze
xmax = delta * nx
ymax = delta * ny
sigmaX = 0.1 * xmax
sigmaY = 0.1 * ymax
TOL = 10**(-8)

xy_array = [[],[]]

for i in range(nx+1):
    for j in range(ny+1):
        xy_array[0].append(i*delta)
        xy_array[1].append(j*delta)

def fun_p1(x,y):
    '''Gęstość 1'''
    return (+1) * math.exp( - ((x - 0.35 * xmax) / sigmaX)**2 - ((y - 0.5 * ymax) / sigmaY)**2 )


def fun_p2(x,y):
    '''Gęstość 2'''
    return (-1) * math.exp( - ((x - 0.65 * xmax) / sigmaX)**2 - ((y - 0.5 * ymax) / sigmaY)**2 )

# xi = delta * i, i=0,1,2....,nx
p = [[fun_p1(delta*i,delta*j) + fun_p2(delta*i,delta*j) for j in range(ny+1)] for i in range(nx+1)]

# Releksacja Lokalna

def metodaRelaksacjiLokalnej(wL):
    V = [[0 for _ in range(ny+1)] for _ in range(nx+1)] # V0 = 0
    for i in range(nx+1):
        V[i][0] = V1

    fun_S = lambda i,j: delta**2 * ( 1/2. * ((V[i+1][j] - V[i][j]) / delta)**2 + 1/2. * ((V[i][j+1] - V[i][j]) / delta)**2 - p[i][j] * V[i][j] )

    S_array = [0]
    while 1:
        for i in range(1,nx):
            for j in range(1,ny):
                V[i][j] = (1-wL)*V[i][j] + (wL/4.)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (delta**2 / epsilon) * p[i][j] )
                


        for j in range(1,ny):
            V[0][j] = V[1][j]
            V[nx][j] = V[nx-1][j]

        S = 0
        for i in range(nx):
            for j in range(ny):
                S += fun_S(i,j)

        S_array.append(S)

        if len(S_array) !=2 and abs((S - S_array[-2]) / S_array[-2]) < TOL:
            break
    
    

    return S_array

def wykres_MRL():
    wL_array = [1.0, 1.4, 1.8, 1.9]
    S_array = []
    for wL in wL_array:
        S_array.append(metodaRelaksacjiLokalnej(wL))
    
    it_array = [[i for i in range(len(S_array[j]))] for j in range(len(S_array))]

    for i in range(len(S_array)):
        plt.plot(it_array[i][1:],S_array[i][1:], label="wL = "+str(wL_array[i]))
    plt.xscale("log")
    plt.title("Metoda relaksacji lokalnej\n Zmiana całki S(it)")
    plt.xlabel("it")
    plt.ylabel("S(it)")
    plt.grid()
    plt.legend()

    plt.show()

# Releksacja Globalna

def metodaRelaksacjiGlobalnej(wG):
    # Vn - nowa, Vs - stara
    Vn = [[0 for _ in range(ny+1)] for _ in range(nx+1)]
    Vs = [[0 for _ in range(ny+1)] for _ in range(nx+1)]

    # fun_S = lambda i,j: delta**2 * ( 1/2. * ((Vn[i+1][j] - Vn[i][j]) / delta)**2 + 1/2. * ((Vn[i][j+1] - Vn[i][j]) / delta)**2 - p[i][j] * Vn[i][j] )

    for i in range(nx+1):
        Vn[i][0] = V1
        Vs[i][0] = V1

    S_array = [0]

    while 1:
        for i in range(1,nx):
            for j in range(1,ny):
                Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + ((delta**2)/epsilon) * p[i][j])

        for j in range(ny+1):
            Vn[0][j] = Vn[1][j]
            Vn[nx][j] = Vn[nx-1][j]

        Vn2 = [] # do wykresu
        bladRozwiazania = []
        xy_array2 = [[],[]] 

        for i in range(nx+1):
            for j in range(ny+1):
                Vs[i][j] = (1. - wG) * Vs[i][j] + wG * Vn[i][j]
                Vn2.append(Vs[i][j])
        

        for i in range(1,nx):
            for j in range(1,ny):
                err = ( (Vs[i+1][j] - 2*Vs[i][j] + Vs[i-1][j]) / (delta**2) + (Vs[i][j+1] - 2*Vs[i][j] + Vs[i][j-1]) / (delta**2) ) + p[i][j] / epsilon
                bladRozwiazania.append(err)
                xy_array2[0].append(i*delta)
                xy_array2[1].append(j*delta)

        S = 0
        for i in range(nx):
            for j in range(ny):
                S += delta**2 * ( 1/2. * ((Vs[i+1][j] - Vs[i][j]) / delta)**2 + 1/2. * ((Vs[i][j+1] - Vs[i][j]) / delta)**2 - p[i][j] * Vs[i][j] )

        S_array.append(S)

        if len(S_array) !=2 and abs((S - S_array[-2]) / S_array[-2]) < TOL:
            break

    # print(len(xy_array2[0]), len(xy_array2[1]), len(bladRozwiazania))

    plt.figure()
    plt.tricontourf(xy_array[0],xy_array[1],Vn2, levels=np.linspace(min(Vn2),max(Vn2),999))
    plt.colorbar(ticks=np.linspace(min(Vn2),max(Vn2),11))

    plt.title("Releksacja globalna\n Zrelaksowany potencjał V(x,y), wG = "+str(wG))
    plt.show()

    plt.figure()
    plt.tricontourf(xy_array2[0],xy_array2[1],bladRozwiazania, levels=np.linspace(min(bladRozwiazania),max(bladRozwiazania),999))
    plt.colorbar(ticks=np.linspace(min(bladRozwiazania),max(bladRozwiazania),11))

    plt.title("Releksacja globalna\n Błąd rozwiązania, wG = "+str(wG))
    plt.show()



    return S_array

if __name__ == "__main__":
    # wykres_MRL()
    
    wG_array = [0.6,1.0]
    S_array = []
    for wG in wG_array:
        S_array.append(metodaRelaksacjiGlobalnej(wG))
        # print("a")

    it_array = [[i for i in range(len(S_array[j]))] for j in range(len(S_array))]

    for i in range(len(S_array)):
        plt.plot(it_array[i][1:],S_array[i][1:], label="wG = "+str(wG_array[i]))
    plt.xscale("log")
    plt.title("Metoda relaksacji globalnej\n Zmiana całki S(it)")
    plt.xlabel("it")
    plt.ylabel("S(it)")
    plt.grid()
    plt.legend()

    plt.show()
