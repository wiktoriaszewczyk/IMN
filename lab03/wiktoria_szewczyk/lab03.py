import math
import matplotlib.pyplot as plt

def wykres(x1, y1, x2, y2, txtTitle, txtXlabel, txtYlabel):
    '''Generowanie wykresow'''

    plt.plot(x1,y1, color='r', label='TOL = 10^(-2)')
    plt.plot(x2,y2, color='g', label='TOL = 10^(-5)')
    plt.title(txtTitle)
    plt.xlabel(txtXlabel)
    plt.ylabel(txtYlabel)
    plt.grid()
    plt.legend()

    plt.show()

def funf(v):
    '''f(t,x,v) = v'''
    return v

def fung(x,v,alpha):
    '''g(t,x,v) = alpha(1 - x^2)v - x'''
    return alpha * (1.0 - x**2) * v - x

def metodaTrapezow(xn, vn, dt, alpha):
    ''' ma zwracać x_(n+1) i v_(n+1) '''
    d = 10**(-10)

    xn1 = xn
    vn1 = vn

    while True:
        F = xn1 - xn -(dt / 2.0) * (funf(vn) + funf(vn1))
        G = vn1 - vn -(dt / 2.0) * (fung(xn,vn,alpha) + fung(xn1,vn1,alpha))

        a11 = 1
        a12 = -dt / 2.0
        a21 = (-dt / 2.0) * (-2 * alpha * xn1 * vn1 - 1)
        a22 = 1 - (dt / 2.0) * alpha * (1 - xn1**2)

        dx = (-F * a22 + G * a12) / (a11 * a22 - a12 * a21)
        dv = (a11 * (-G) - a21 * (-F)) / (a11 * a22 - a12 * a21)

        xn1 += dx
        vn1 += dv  

        if (abs(dx) < d and abs(dv) < d):
            break

    return xn1, vn1

def metodaRK2(xn, vn, dt, alpha):
    ''' ma zwracać x_(n+1) i v_(n+1) '''
    k1x = vn
    k1v = alpha * (1.0 - xn**2) * vn - xn

    k2x = vn + dt * k1v
    k2v = alpha * (1.0 - (xn + dt * k1x)**2) * (vn + dt * k1v) - (xn + dt * k1x)

    xn1 = xn + (dt / 2.0) * (k1x + k2x)
    vn1 = vn + (dt / 2.0) * (k1v + k2v)

    return xn1, vn1

def wnetrze_kontrolaKrokuCzasowego(tab_t, tab_dt, tab_x, tab_v, x0, v0, dt0, tmax, p, S, alpha, TOL, metoda, txtTitle):
    t = 0
    dt = dt0
    xn = x0
    vn = v0

    while True:
        # dwa kroki dt
        x2_n1, v2_n1 = metoda(xn, vn, dt, alpha)
        x2_n2, v2_n2 = metoda(x2_n1, v2_n1, dt, alpha)
        # jeden krok 2*dt
        x1_n2, v1_n2 = metoda(xn, vn, 2.0*dt, alpha)
        Ex = (x2_n2 - x1_n2) / ((2**p) - 1.0)
        Ev = (v2_n2 - v1_n2) / ((2**p) - 1.0)
        
        if max( abs(Ex),abs(Ev) ) < TOL:
            t = t + 2*dt
            xn = x2_n2
            vn = v2_n2
            tab_t.append(t)
            tab_dt.append(dt)
            tab_x.append(xn)
            tab_v.append(vn)

        dt = (((S * TOL) / max( abs(Ex),abs(Ev) ))**(1.0 / (p + 1))) * dt
        # print(t, xn, vn, dt)

        if(t >= tmax):
            break

def kontrolaKrokuCzasowego(x0, v0, dt0, tmax, p, S, alpha, TOL, metoda, txtTitle):
    '''algorytm kontroli kroku czasowego'''

    # wyniki TOL = 10^-2
    tab_t1 = []
    tab_dt1 = []
    tab_x1 = []
    tab_v1 = []

    # wyniki TOL = 10^-5
    tab_t2 = []
    tab_dt2 = []
    tab_x2 = []
    tab_v2 = []

    wnetrze_kontrolaKrokuCzasowego(tab_t1, tab_dt1, tab_x1, tab_v1, x0, v0, dt0, tmax, p, S, alpha, TOL[0], metoda, txtTitle)
    wnetrze_kontrolaKrokuCzasowego(tab_t2, tab_dt2, tab_x2, tab_v2, x0, v0, dt0, tmax, p, S, alpha, TOL[1], metoda, txtTitle)

    #wykresy x(t), v(t), dt(t), v(x) dla obu TOL na jednym rysunku
    wykres(tab_t1, tab_x1, tab_t2, tab_x2, txtTitle, "t", "x(t)")
    wykres(tab_t1, tab_v1, tab_t2, tab_v2, txtTitle, "t", "v(t)")
    wykres(tab_t1, tab_dt1, tab_t2, tab_dt2, txtTitle, "t", "dt(t)")
    wykres(tab_x1, tab_v1, tab_x2, tab_v2, txtTitle, "x", "v(x)")

if __name__ == "__main__":
    '''Rozwiązywanie równania różniczkowego 2 rzędu oscylatora Van der Pola'''
    # parametry startowe
    x0 = 0.01
    v0 = 0
    dt0 = 1
    S = 0.75
    p = 2   # rząd dokładności obu metod
    tmax = 40
    alpha = 5
    TOL = [10**(-2), 10**(-5)]

    kontrolaKrokuCzasowego(x0, v0, dt0, tmax, p, S, alpha, TOL, metodaTrapezow, "Metoda trapezow")
    kontrolaKrokuCzasowego(x0, v0, dt0, tmax, p, S, alpha, TOL, metodaRK2, "Metoda RK2")

