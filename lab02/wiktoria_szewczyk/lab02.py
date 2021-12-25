import math
import numpy as np
import matplotlib.pyplot as plt

def wykres(t, u, z, text):
    '''Generowanie wykresow'''

    plt.plot(t,u, color='r', label='u(t)')
    plt.plot(t,z, color='g', label='z(t)')
    plt.title(text)
    plt.xlabel('t')
    plt.ylabel('u(t),z(t)')
    plt.grid()
    plt.legend()

    plt.show()

def Picard(u_n, dt, alpha, beta, u_n1):
    '''Metoda Picarda'''
    return u_n + (dt/2) * ( (alpha * u_n - beta * u_n * u_n) + (alpha * u_n1 - beta * u_n1 * u_n1) )

def Newton(u_n, dt, alpha, beta, u_n1):
    '''Metoda Newtona'''
    return u_n1 - ( (u_n1 - u_n - (dt/2) * ( (alpha * u_n - beta * u_n * u_n) + (alpha * u_n1 - beta * u_n1 * u_n1) )) / (1 - (dt/2) * (alpha - 2*beta * u_n1)) )

def MetodaTrapezow(beta, N, gamma, t_max, dt, u0, TOL):
    '''Metoda Picarda i metoda Newtona''' 

    u_P = [u0]      # nosiciele metoda Picarda
    u_N = [u0]      # nosiciele metoda Newtona
    alpha = (beta * N) - gamma

    t_array = np.arange(0,t_max+dt,dt)


    for t in t_array[1:]:
        u_mi0 = u_P[-1]  # u w trakcie obliczen mi (stare)
        u_mi1 = u_P[-1]  # u w trakcie obliczen mi + 1 (nowe)
        for mi in range(21):    # mi <= 20
            u_mi1 = Picard(u_P[-1], dt, alpha, beta, u_mi0)
            if abs(u_mi1 - u_mi0) < TOL:
                break
        u_P.append(u_mi1)

    for t in t_array[1:]:
        u_mi0 = u_N[-1]  # u w trakcie obliczen mi (stare)
        u_mi1 = u_N[-1]  # u w trakcie obliczen mi + 1 (nowe)
        for mi in range(21):    # mi <= 20
            u_mi1 = Newton(u_N[-1], dt, alpha, beta, u_mi0)
            if abs(u_mi1 - u_mi0) < TOL:
                break
        u_N.append(u_mi1)

    z_P = [N - u for u in u_P]
    z_N = [N - u for u in u_N]

    wykres(t_array,u_P,z_P,"Metoda Picarda")
    wykres(t_array,u_N,z_N,"Metoda Newtona")


def RK2(beta, N, gamma, t_max, dt, u0, TOL):
    '''Niejawna metoda RK2'''
    alpha = (beta * N) - gamma

    # tablica Butchera
    a11 = 1/4
    a12 = 1/4 - ((3**(1/2)) / 6)
    a21 = 1/4 + ((3**(1/2)) / 6)
    a22 = 1/4
    b1 = 1/2
    b2 = 1/2
    c1 = 1/2 - ((3**(1/2)) / 6)
    c2 = 1/2 + ((3**(1/2)) / 6)

    fun_F1 = lambda U1, U2, u_n : ( U1 - u_n - dt * (a11 * (alpha * U1 - beta * U1 * U1) + a12 * (alpha * U2 - beta* U2 * U2)) )
    fun_F2 = lambda U1, U2, u_n : ( U2 - u_n - dt * (a21 * (alpha * U1 - beta * U1 * U1) + a22 * (alpha * U2 - beta* U2 * U2)) )

    fun_m11 = lambda U1 : 1 - dt * a11 * (alpha - 2 * beta * U1)
    fun_m12 = lambda U2 : -dt * a12 * (alpha - 2 * beta * U2)
    fun_m21 = lambda U1 : -dt * a21 * (alpha - 2 * beta * U1)
    fun_m22 = lambda U2 : 1 - dt * a22 * (alpha - 2 * beta * U2)

    fun_du1 = lambda U1, U2, u_n: (fun_F2(U1, U2, u_n) * fun_m12(U2) - fun_F1(U1, U2, u_n) * fun_m22(U2)) / (fun_m11(U1) * fun_m22(U2) - fun_m12(U2) * fun_m21(U1))
    fun_du2 = lambda U1, U2, u_n: (fun_F1(U1, U2, u_n) * fun_m21(U2) - fun_F2(U1, U2, u_n) * fun_m11(U2)) / (fun_m11(U1) * fun_m22(U2) - fun_m12(U2) * fun_m21(U1))

    fun_uf = lambda x: alpha * x - beta * x * x
    fun_u = lambda u_n, t, u1, u2 : (u_n + dt * ((b1 * fun_uf(u1)) + (b2 * fun_uf(u2)))) 

    u = [u0]

    t_array = np.arange(0,t_max+dt,dt)

    for t in t_array[1:]:
        u10 = u[-1]
        u20 = u[-1]
        for mi in range(21): 
            du1 = fun_du1(u10,u20,u[-1])
            du2 = fun_du2(u10,u20,u[-1])
            u11 = u10 + du1 
            u21 = u20 + du2
            if abs(u11 - u10) < TOL or abs(u21 - u20) < TOL:
                break

        u.append(fun_u(u[-1], t, u11, u21))

    z = [N - i for i in u]

    wykres(t_array,u,z,"Niejawna metoda RK2")


if __name__ == "__main__":

    # Dane:
    beta = 0.001    # czestosc kontaktow zarazonych ze zdrowymi
    N = 500         # populacja
    gamma = 0.1     # sredni czas trwania choroby
    t_max = 100     # czas maksymalny
    dt = 0.1        # delta czasu
    u0 = 1          # co najmniej jeden zarazony na poczatku
    TOL = 10**(-6)  # warunek stopu dodatkowo mi <= 20     

    # MetodaTrapezow(beta, N, gamma, t_max, dt, u0, TOL)
    RK2(beta, N, gamma, t_max, dt, u0, TOL)