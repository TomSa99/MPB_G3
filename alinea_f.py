

# ---------------------correr no ficheiro jupyter--------------------------

import sympy as sp
# Valores iniciais das variáveis
# Y = [X, S, A, P, V]
y0 = [4., 0., 0., 0., 5.]


# Parâmetros minimizantes (aprox.) -> k3, umax3, ks3
params=[13.50945035,0.5754241,1.25010977]
k3,umax3,ks3 = params

# Parâmetros k (exceto k3)
k1,k2,k4,k5,k6,k7,k8,k9,k10,k11 = [5.164,27.22, 4.382, 2.074, 10.89, 4.098, 2.283, 17.01, 4.576, 12]

# Parâmetros umax (exceto umax3)
umax1 = 0.2
umax2 = 0.5

# Parâmetros ks (exceto ks3)
ks2 = 0.1
ks1 = 0.1


# Consider using sympy.symbols to create algebric variables to be used on the derivatives (X, S, k1, ks1, ...)

# Criar símbolos algébricos para as variáveis utilizando a função "symbols" do package SymPy
(X,S,A,P,V,K1,K2,K3,Umax1,Umax2,Umax3,Ks1,Ks2,Ks3,Fe,Se) = sympy.symbols('X,S,A,P,V,k1,k2,k3,umax1,umax2,umax3,Ks1,Ks2,Ks3,Fe,Se', real=True)

# Parâmetros mu
u1 = Umax1 * (S / (Ks1 + S))
u2 = Umax2 * (S / (Ks2 + S))
u3 = Umax3 * (A / (Ks3 + A))

# ODEs relativas a X e S
dxdt = u1*X + u2*X + u3*X - Fe/V*X
dsdt = -K1*u1*X - K2*u2*X - Fe/V*S + Fe/V*Se


# Calcular as derivadas parciais de dX/dt e dS/dt em função de k3, umax3 e ks3, usando a função "diff" do package SymPy
dXk3 = sp.diff(dxdt, k3)
dXks3 = sp.diff(dxdt, Ks3)
dXumax3 = sp.diff(dxdt, Umax3)

dSk3 = sp.diff(dsdt, k3)
dSks3 = sp.diff(dsdt, Ks3)
dSumax3 = sp.diff(dsdt, Umax3)


# Criar funções para substituir incógnitas das funções das derivadas parciais, usando a função "lambdify" do package SymPy
dXk3 = sp.lambdify((), dXk3, "numpy")
dXks3 = sp.lambdify((X, S, Ks3, Umax3), dXks3, "numpy")
dXumax3 = sp.lambdify((X, S, Ks3), dXumax3, "numpy")
dSk3 = sp.lambdify((X, S, Ks3, Umax3), dSk3, "numpy")
dSks3 = sp.lambdify((X, S, Ks3, Umax3, K3), dSks3, "numpy")
dSumax3 = sp.lambdify((X, S, Ks3, K3), dSumax3, "numpy")


t0 = 0    #Tempo inicial
t1 = 24 #Tempo final
dt = 0.01 #Passo

# ode
res = ode(jm109).set_integrator("lsoda", method="bdf", lband=0) # nsteps=10000
res.set_initial_value(y0, t0).set_f_params(params)

# Listas para guardar os valores das sensibilidades ao longo do tempo
T,dxk3,dxks3,dxumax3,dsk3,dsks3,dsumax3 = [],[],[],[],[],[],[]

# Resolução da ode para cada passo de tempo
while res.successful() and res.t < t1:
    result = res.integrate(res.t+dt)
    x,s,a,p,v = result
    
    T.append(res.t)
    
    dxk3.append(dXk3())
    dxks3.append(dXks3(x,s,ks3,umax3))
    dxumax3.append(dXumax3(x,s,ks3))
    
    dsk3.append(dSk3(x,s,ks3,umax3))
    dsks3.append(dSks3(x,s,ks3,umax3,k3))
    dsumax3.append(dSumax3(x,s,ks3,k3))


# Criar gráficos das sensibilidades de X em função de k3, umax3 e ks3
fig, ax = plt.subplots(3,2)
fig.set_figheight(8)
fig.set_figwidth(10)
ax[0,0].plot(T, dxk3, linewidth=2, linestyle='solid', label='dX/dK3', color='lightblue')
ax[0,0].set_title('Sensibilidade do X')
ax[0,0].get_yaxis().set_label_coords(-0.20,0.5)
ax[0,0].legend()
ax[0,0].set_ylabel('Sensibilidade a K1')
ax[1,0].get_yaxis().set_label_coords(-0.2,0.5)
ax[1,0].plot(T, dxumax3, linewidth=2, linestyle='solid', label='dX/dumax3', color='blue')
ax[1,0].legend()
ax[1,0].set_ylabel('Sensibilidade a umax1')
ax[2,0].plot(T, dxks3, linewidth=2, linestyle='solid', label='dX/dKs3', color='darkblue')
ax[2,0].legend()
ax[2,0].set_ylabel('Sensibilidade a Ks1')
ax[2,0].get_yaxis().set_label_coords(-0.2,0.5)
ax[2,0].set_xlabel('Tempo (h)')

# Criar gráficos das sensibilidades de S em função de k3, umax3 e ks3
ax[0,1].plot(T, dsk3, linewidth=2, linestyle='solid', label='dS/dK3', color='salmon')
ax[0,1].set_title('Sensibilidade do S')
ax[0,1].legend()
ax[1,1].plot(T, dsumax3, linewidth=2, linestyle='solid', label='dS/dumax3', color='red')
ax[1,1].legend()
ax[2,1].plot(T, dsks3, linewidth=2, linestyle='solid', label='dS/dKs3', color='darkred')
ax[2,1].legend()
ax[2,1].set_xlabel('Tempo (h)')

fig.subplots_adjust(left=0.2, wspace=0.2)
