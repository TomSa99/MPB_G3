import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Model bl21 fed batch
def bl21_fb(t, y, params):
    X, S, A, P, V = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se = params
    u1 = 0.2 * (S / (0.1 + S))
    u2 = 0.5 * (S / (0.1 + S))
    u3 = 0.1 * (A / (0.6 + A))
    D= Fe / V
    reac = [u1*X + u2*X + u3*X - D*X, -k1*u1*X - k2*u2*X - D*S + D*Se, k3*u2*X - k4*u3*X - D*A, k11*u1*X - D*P, Fe] #X, S, A, P, V || usando as reacoes de crescimento
    return reac

# Initial conditions
X0= 4 #g/L
S0= 0 #g/L
A0= 0 #g/L
P0= 0 #g/L
V= 5 #L

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
V0= 5 #o volume inicial não se altera
Fe= 0.8 #L/h || caudal de entrada || 350 g/L glucose
Se= 450 #concentração do substrato de entrada g/L

#lista com os valores iniciais fornecida a func
y0= [X0, S0, A0, P0, V]

#lista com os parametros fornecida a func
params= [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, V0, Fe, Se]

t0= 0 #tempo inicial
t= 24 #tempo final
dt= 0.001 #intervalo de tempo entre reads


# Call the ODE solver
r = ode(bl21_fb).set_integrator('lsoda', method='bdf', lband=0) #lband é o limite inferior -- para nao haver valores negativos
r.set_initial_value(y0, t0).set_f_params(params)


#storing variables
T, x, s, a, p, v= [], [], [], [], [], []

while r.successful() and r.t < t:
    time= r.t + dt
    T.append(r.t)
    x.append(r.integrate(time)[0])
    s.append(r.integrate(time)[1])
    a.append(r.integrate(time)[2])
    p.append(r.integrate(time)[3])
    v.append(r.integrate(time)[4])
    #print(time, r.integrate(time))


# Plotting
T = np.array([0] + T)
x = np.array([y0[0]] + x)
s = np.array([y0[1]] + s)
a = np.array([y0[2]] + a)
p = np.array([y0[3]] + p)
v = np.array([y0[4]] + v)

#plot
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, x, label='Biomassa', color='blue')
ax.plot(T, s, label='Substrato', color='red')
ax.plot(T, a, label='Acetato', color='green')
ax.plot(T, p, label='Proteína Recombinante', color='pink')
ax.plot(T, v, label='Volume (L)', color='purple')
ax.legend(loc='best')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo BL21 (Fed-Batch)')
plt.grid()
plt.show()

fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(10)
ax.plot(T, x, label='Biomassa', color='blue')
ax.plot(T, s, label='Substrato', color='red')
ax.plot(T, a, label='Acetato', color='green')
ax.plot(T, p, label='Proteína Recombinante', color='pink')
ax.plot(T, v, label='Volume (L)', color='purple')
ax.legend(loc='best')
ax.set_xlabel('Tempo (h)')
ax.set_ylabel('Concentração (g/L)')
ax.set_title('Modelo BL21 (Fed-Batch)')
plt.xlim(0., 2.5)
plt.ylim(0., 10)
plt.grid()
plt.show()
