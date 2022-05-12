import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Model bl21 batch
def bl21_b(t, y, params):
    X, S, A, P = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0 = params
    u1 = 0.2 * (S / (0.1 + S))
    u2 = 0.5 * (S / (0.1 + S))
    u3 = 0.1 * (A / (0.6 + A))
    reac = [u1*X + u2*X + u3*X, -k1*u1*X - k2*u2*X, k3*u2*X - k4*u3*X, k11*u1*X] #X, S, A, P || usando as reacoes de crescimento
    return reac

# Initial conditions
X0= 4 #g/L
S0= 10 #g/L
A0= 0 #g/L
P0= 0 #g/L

# Parameters
k1= 5.164
k2= 27.22
k3= 12.90
k4= 4.382
k5= 2.074
k6= 10.89
k7= 4.098
k8= 2.283
k9= 17.01
k10= 4.576
k11= 12.0
V0= 3 #como é batch o volume não se altera

#lista com os valores iniciais fornecida a func
y0= [X0, S0, A0, P0]

#lista com os parametros fornecida a func
params= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0]

t0= 0 #tempo inicial
t= 5 #tempo final
dt= 0.001 #intervalo de tempo entre reads

# Call the ODE solver
r = ode(bl21_b).set_integrator('lsoda', method='bdf', lband= 0) #lband é o limite inferior -- para nao haver valores negativos
r.set_initial_value(y0, t0).set_f_params(params)


#storing variables
T, x, s, a, p= [], [], [], [], []

while r.successful() and r.t < t:
    time= r.t + dt
    T.append(r.t)
    x.append(r.integrate(time)[0])
    s.append(r.integrate(time)[1])
    a.append(r.integrate(time)[2])
    p.append(r.integrate(time)[3])
    #print(time, r.integrate(time))


# using the storing variables to plot
T = np.array([0] + T)
x = np.array([y0[0]] + x)
s = np.array([y0[1]] + s)
a = np.array([y0[2]] + a)
p = np.array([y0[3]] + p)


#plot
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, x, label='Biomassa', color='blue')
ax.plot(T, s, label='Substrato', color='red')
ax.plot(T, a, label='Acetato', color='green')
ax.plot(T, p, label='Proteína Recombinante', color='pink')
ax.plot(T, [V0] * len(T), label='Volume (L)', color='purple')
ax.legend(loc='best')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo BL21 (Batch)')
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, x, label='Biomassa', color='blue')
ax.plot(T, s, label='Substrato', color='red')
ax.plot(T, a, label='Acetato', color='green')
ax.plot(T, p, label='Proteína Recombinante', color='pink')
ax.plot(T, [V0] * len(T), label='Volume (L)', color='purple')
ax.legend(loc='best')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo BL21 (Batch)')
plt.xlim(0., 0.5)
plt.grid()
plt.show()
