import math
import matplotlib.pyplot as plt

l = -1  # lambda

def fun_y_dok(t):
    '''Rozwiazanie: y(t) = e^(lambda*t)'''
    return math.exp(l * t)

# argumenty funkcji
t1 = [i/100 for i in range(0,501,1)]    # delta_t 0.01
t2 = [i/10  for i in range(0,51,1)]     # delta_t 0.1
t3 = [i     for i in range(0,6,1)]      # delta_t 11

# rozwiazania analityczne
y_dok1 = [fun_y_dok(i) for i in t1]
y_dok2 = [fun_y_dok(i) for i in t2]
y_dok3 = [fun_y_dok(i) for i in t3]

def wykresy(y_num1, y_num2, y_num3, bg1, bg2, bg3, name):
    '''Wykresy do zadan 1, 2, 3'''
    figure, axis = plt.subplots(2)

    figure.suptitle(name)

    axis[0].plot(t1,y_dok1, color='k', label='y_dok')
    axis[0].plot(t1,y_num1, color='r', label='y_num1')
    axis[0].plot(t2,y_num2, color='g', label='y_num2')
    axis[0].plot(t3,y_num3, color='b', label='y_num3')
    axis[0].set_title('Rozwiazania')
    axis[0].grid()
    axis[0].legend()

    axis[1].plot(t1,bg1, color='r', label='bg1')
    axis[1].plot(t2,bg2, color='g', label='bg2')
    axis[1].plot(t3,bg3, color='b', label='bg3')
    axis[1].set_title('Zmiany bledu globalnego')
    axis[1].grid()
    axis[1].legend()

    for ax in axis.flat:
        ax.set(xlabel='t', ylabel='y')

    plt.show()

def fun_y_Euler(y, dt):
    '''Metoda jawna Eulera'''
    return y + (dt * l * y)

def fun_y_RK2(y, dt):
    '''Metoda jawna RK2 (trapezow)'''
    k1 = l * y
    k2 = l * (y + (dt * k1))
    return y + ( (dt / 2) * (k1 + k2) )

def fun_y_RK4(y,dt):
    '''Metoda jawna RK4'''
    k1 = l * y
    k2 = l * (y + ( (dt/2) * k1 ))
    k3 = l * (y + ( (dt/2) * k2 ))
    k4 = l * (y + ( dt * k3 ))
    return y + ( (dt/6) * (k1 + 2*k2 + 2*k3 + k4) )

def problem_autonomiczny(fun, name):
    '''
        rozwiazywanie zadan 1, 2, 3
        fun - wybor metody (Eulera, trapezow, RK4)
        name - nazwa metody (do tytulu wykresu)
    '''

    # rozwiazania numeryczne
    y_num1 = [1]
    y_num2 = [1]
    y_num3 = [1]

    for t in t1[1:]:
        y_num1.append( fun(y_num1[-1], 0.01) )

    for t in t2[1:]:
        y_num2.append( fun(y_num2[-1], 0.1) )

    for t in t3[1:]:
        y_num3.append( fun(y_num3[-1], 1) )

    # blad globalny
    bg1 = [y_num1[i] - y_dok1[i] for i in range(len(y_num1))]
    bg2 = [y_num2[i] - y_dok2[i] for i in range(len(y_num2))]
    bg3 = [y_num3[i] - y_dok3[i] for i in range(len(y_num3))]

    wykresy(y_num1, y_num2, y_num3, bg1, bg2, bg3, name)

def fun_V(t, wv):
    return 10 * math.sin( wv * t)

def RK4(Q, I, t_array, dt, R, L, C, wv):
    '''Obliczanie Q i I (jawny schemat RK4) dla czestosci zrodla wv'''
    for t in t_array[1:]:
        k1_Q = I[-1]
        k1_I = (fun_V(t, wv) / L) - ((1 / (L * C)) * Q[-1]) - ((R / L) * I[-1])
        k2_Q = I[-1] + ((dt / 2) * k1_I)
        k2_I = (fun_V((t + (dt / 2)), wv) / L) - ((1 / (L * C)) * (Q[-1] + ((dt / 2) * k1_Q))) - ((R / L) * (I[-1] + ((dt / 2) * k1_I)))
        k3_Q = I[-1] + ((dt /2) * k2_I)
        k3_I = (fun_V((t + (dt / 2)), wv) / L) - ((1 / (L * C)) * (Q[-1] + ((dt / 2) * k2_Q))) - ((R / L) * (I[-1] + ((dt / 2) * k2_I)))
        k4_Q = I[-1] + (dt * k3_I)
        k4_I = (fun_V((t + dt), wv) / L) - ((1 / (L * C)) * (Q[-1] + (dt * k3_Q))) - ((R / L) * (I[-1] + (dt * k3_I)))

        Q.append(Q[-1] + ((dt / 6) * (k1_Q + 2*k2_Q + 2*k3_Q + k4_Q)))
        I.append(I[-1] + ((dt / 6) * (k1_I + 2*k2_I + 2*k3_I + k4_I)))


def RLC():
    '''Rozwiazanie ostatniego zadania'''

    # znane parametry
    dt = 10**(-4)
    R = 100
    L = 0.1
    C = 0.001
    w0 = 1 / (L * C)
    T0 = (2 * math.pi) / w0

    # czestosci zrodla
    wv1 = 0.5 * w0
    wv2 = 0.8 * w0
    wv3 = 1.0 * w0
    wv4 = 1.2 * w0


    # argumenty funkcji [0, 4*T0]
    t_array = [0]
    while(t_array[-1] <= (4 * T0)):
        t_array.append(round(( t_array[-1] + dt), 4))
    t_array.pop()

    # warunki poczatkowe
    Q1 = [0]
    I1 = [0]
    Q2 = [0]
    I2 = [0]
    Q3 = [0]
    I3 = [0]
    Q4 = [0]
    I4 = [0]

    # obliczenia
    RK4(Q1, I1, t_array, dt, R, L, C, wv1)
    RK4(Q2, I2, t_array, dt, R, L, C, wv2)
    RK4(Q3, I3, t_array, dt, R, L, C, wv3)
    RK4(Q4, I4, t_array, dt, R, L, C, wv4)

    # wykresy

    # figure, axis = plt.subplots(2)

    # figure.suptitle('RLC')

    # axis[0].plot(t_array,I1, color='k', label='0.5 omega0')
    # axis[0].plot(t_array,I2, color='r', label='0.8 omega0')
    # axis[0].plot(t_array,I3, color='g', label='1.0 omega0')
    # axis[0].plot(t_array,I4, color='b', label='1.2 omega0')
    # axis[0].set_title('I')
    # axis[0].grid()
    # axis[0].legend()

    # axis[1].plot(t_array,Q1, color='k', label='0.5 omega0')
    # axis[1].plot(t_array,Q2, color='r', label='0.8 omega0')
    # axis[1].plot(t_array,Q3, color='g', label='1.0 omega0')
    # axis[1].plot(t_array,Q4, color='b', label='1.2 omega0')
    # axis[1].set_title('Q')
    # axis[1].grid()
    # axis[1].legend()

    # axis.flat[0].set(xlabel='t', ylabel='I')
    # axis.flat[1].set(xlabel='t', ylabel='Q')


    # nie wie czemu ale chyba zapisało sie na odwrot 
    plt.plot(t_array,I1, color='k', label='0.5 omega0')
    plt.plot(t_array,I2, color='r', label='0.8 omega0')
    plt.plot(t_array,I3, color='g', label='1.0 omega0')
    plt.plot(t_array,I4, color='b', label='1.2 omega0')
    plt.title('RLC')
    plt.xlabel('t')
    plt.ylabel('Q')
    plt.grid()
    plt.legend()


    plt.show()


if __name__ == "__main__":
    problem_autonomiczny(fun_y_Euler, "Metoda jawna Eulera")
    problem_autonomiczny(fun_y_RK2, "Metoda jawna RK2 (trapezow)")
    problem_autonomiczny(fun_y_RK4, "Metoda jawna RK4")
    RLC()
