# Version: Beta
# @author: Leonel Gerardo González Pérez
# Este codigo esta basado en en el excel de Richard Nakka 'SRM_2013': https://www.nakka-rocketry.net/soft/SRM_2023.zip

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as C
from scipy.interpolate import CubicSpline
from scipy.integrate import odeint

# %%
def mach(P, Po, k):
    return np.sqrt((2/(k-1))*((Po/P)**((k-1)/k)- 1))

def RatArea(M1, M2, k):
    a=(k+1)/(2*(k-1))
    b=(k-1)/2
    return (M2/M1)*(((1+b*(M1**2))/(1+b*(M2**2)))**a)

def A1A2(Po,P,k):
    a = (k-1)/k
    b = (((k+1)/2)**(1/(k-1)))*((P/Po)**(1/k))
    c = np.sqrt((k+1)/(k-1)*(1-(P/Po)**a))
    return b*c

def Cf_(Pe, P0, Pa, At, Ae, k):
    a = (k+1)/(k-1)
    return np.sqrt((2*k**2/(k-1)*(2/(k+1))**a)*(1-(Pe/P0)**((k-1)/k))) + (Pe - Pa)*Ae/(P0*At)

def diametro(A):
    return np.sqrt((4*A)/C.pi)

def ex_vel(R, To, k):
    a = (k+1)/(k-1)
    c = np.sqrt(R*To/(k*(2/(k+1))**a))
    return c

def Ab_(Outer, Core, Ends, p):
    N, Do, do, Lo, = p

    t_web = (Do-do)/2
    x = np.linspace(0,t_web ,1000)
    D = Do-Outer*(2*x)
    d = do + Core*(2*x)
    L = Lo - Ends*(2*x)

    if Outer == 1:
        A_o = N*np.pi*D*L
    else :
        A_o = [0]

    if Core == 1:
        A_c = N*np.pi*d*L
    else:
        A_c = [0]

    if Ends == 1:
        A_e = N*1/2*np.pi*(D**2-d**2)
    else:
        A_e = [0]

    Ab_tot = A_o + A_c + A_e

    return [A_o[0], A_c[0], A_e[0], Ab_tot]

def Kn_(At, Ab):
    #A_o, A_c, A_e = Ab

    Kn = (Ab)/At
    Kno = min(Kn)
    Kn_prom = np.mean(Kn)
    Kn_max = max(Kn)

    return [Kno, Kn_prom, Kn_max]


def Kn_pol(Po):
      Kn_max =  32.954 + 44.108*(Po/145.038) -1.1025*(Po/145.038)**2  # De psi a Mpa
      return Kn_max

def Ve(To, R, k, Pe, Po ):

    Ve = np.sqrt(2*To*R*(k/(k-1))*(1-(Pe/Po)**((k-1)/k)))
    return Ve
    
def mdot(At,Po, k, To, R):
    m = At*Po*np.sqrt((1/(R*To))*k*(2/(k+1))**((k+1)/(k-1)))
    return m

def tail_off(Pab,R, To, At, V0, c, t):
    Pc = Pab*np.exp(-R*To*At*t/(V0*c))
    return Pc


# Variación de la densidad del aire segun la altura.
def densidad_aire(h):
    # Definición de constantes
    rho_0 = 1.225  # kg/m^3
    g = C.g  # m/s^2
    M = 0.02897  # kg/mol
    R = 8.3144598  # J/mol·K
    T0 = 288.15  # K

    rho = rho_0 * np.exp(-g * M * h / (R * T0))
    return rho

# %%
alpha = 15*np.pi/180 # angulo de la divergencia de la tobera 

N_noz = 0.85 # # eficiencia de la tobera
N_com = 0.95 # eficiencia de combustion
N_div = 1/2*(1 + np.cos(alpha)) # factor de corrección por la divergencia
N_po = 0.95 # factor de corrección de la presion de la camara

rho_rat = 0.95


# %%
k = 1.141607     # relación de calores especificos
M = 0.03573117   # kg/mol masa molecular de los gases 
To_T = 1820.426    # K temperatura de combustion teorica
rho_T = 1.8892   # g/cm^3 densidad teorica

To = To_T*N_com # Temperatura "real"
rho = rho_T*rho_rat    # g/cm^3 densidad medida

# Coeficiente y exponente obtenido por Richard Nakka para el propelente KNSU
n = 0.319
a = 8.260

# %%
P0 = 800 #psi, presion de la camara objetivo
Pa = 14.69594878 # psi, presion atmosferica
Pe = Pa #psi, presion de salida de la tobera

# %% [markdown]
# ## Condiciones del Cohete

# %%
theta_0 = np.pi/2 #angulo

# Mach
Mt = 1
Me = mach(Pe,P0,k)

# %%
Cd = 0.333 # coeficiente de arrastre


Dc = 102.26 #mm Diametro de la camara 
Df = 6 # in, diametro fuselaje

Dp =  50 # Diametro del port
In = 0 # ancho de los inhibidores

Af = np.pi*(Df/(2*39.37))**2 # área transversal (6 inch es el diametro del tubo de fuselaje)
Ac = C.pi*(Dc/2)**2 #mm**2 # área de la camara

# 1/7.2 
AtAc_l = 1/7.2
At = Ac*AtAc_l
Dt = diametro(At)

ApAt = np.pi*(Dp/2)**2/At
# 1.6, 1.75


AeAt = RatArea(Me, Mt, k)
Ae = At*AeAt
De = diametro(Ae)

# %%

# Variables del ciclo iterativo

Inc_M = 0.02 # kg, incremento en la masa total por cada iteración
Datos = 5000 # cantidad de datos para la regresión del grano
h0 = 3000 # altura deseada 


# Superficies de quemado

# Si es 1 significa que se toma en cuenta si es 0 no se toma en cuenta, si se pone algun otro número se tendran resultados erroneos
Bs = [0, 1, 1] # Superficies: exterior, nucleo, caras
N = 4 # Número de segmentos del propelente


Cf = Cf_(Pe, P0, Pa, At, Ae, k)*N_noz
Ve = ex_vel(C.R/M, To, k)*Cf # Velocidad de salida de los gases de la tobera

m = 60 # masa del cohete sin propelente
mp0 = m*(np.exp(np.sqrt(2*C.g*h0)/Ve)-1) # masa minima
m_p = [mp0] # se utiliza la formula del método que no se considera la fuerza de arrastre para la primera iteración
M_tot = [m + mp0]



# Ciclo iterativo



i = 0
h = [0]
while h[i] < h0:
    print(f'Iteración = {i + 1}')


    # Área de quemado
    V0 = m_p[i]*1000/rho # cm^3 Volumen con densidad experiemental del grano
    Va0 = (V0)/(1-ApAt*At*4/(np.pi*Dc**2)) # cm^3 Volumen disponible de la camara

    L0 = Va0*10**3/Ac # mm  Longitud del grano con densidad experimental

    p0 = [N, Dc-In, Dp, L0/N] #  N, D, d, L, At, = p
    Ab0 = Ab_(Bs[0],Bs[1],Bs[2], p0) # área de quemado

    Kn_max= Kn_pol(P0)
    Kn = Kn_(At, Ab0[-1])

    At2 = max(Ab0[-1])/Kn_max

    # %%
    
    t_web0  = np.array([(Dc-In-Dp)/2])
    Xinc    = np.linspace(0, t_web0[0]/(1+Bs[0]*Bs[1]), Datos)
    t_web   = np.append(t_web0, t_web0[0]-Xinc)

    D = Dc - In - Bs[0]*2*Xinc
    d = Dp + Bs[1]*2*Xinc
    L = L0 - Bs[2]*N*2*Xinc

    Vc   = np.pi*(Dc/2)**2*L0/(1000**3) # m^3
    V_G  = 1/4*np.pi*(D**2-d**2)*L/(1000**3) # m^3
    V_F  = Vc - V_G

    P_a          = 0.101325

    Po_abs1     = [] # En Pa
    Po_abs2     = [P_a] # En Mpa
    m_grain     = rho*1000*V_G # kg
    m_gen       = [0]
    m_noz       = [0]
    m_sto       = [0]
    mass_sto    = [0]
    rho_prod    = []
    t1          = [0]
    R           = [a*(Po_abs2[0]**n)]

    for l in range(len(Xinc)-1):
        rho_prod.append(mass_sto[l]/V_F[l])

        Po_abs1.append(rho_prod[l]*C.R/M*To+P_a*10**6)
        Po_abs2.append(Po_abs1[l]/(10**6))
        R.append(a*(Po_abs2[l])**n)
        t1.append(Xinc[1]/R[l+1] + t1[l]) # tiempo de quemado 
        m_gen.append((m_grain[l]-m_grain[l+1])/(t1[l+1]-t1[l]))
        
        if m_gen[l+1] < mdot(At/(10**6),(Po_abs2[l+1]-P_a)*10**6, k,To,C.R/M):
            if Po_abs2[l+1] > 0:
                m_noz.append( mdot(At/(10**6),(Po_abs2[l+1]-P_a)*10**6, k,To,C.R/M))
            else :
                m_noz.append(0)
        else :
            m_noz.append(mdot(At/(10**6),(Po_abs2[l+1]-P_a)*10**6, k,To,C.R/M))
        
        m_sto.append(m_gen[l+1]-m_noz[l+1])
        

        if m_sto[l+1]*(t1[l+1]-t1[l]) + mass_sto[l] < 0 :
            mass_sto.append(0)
        else :
            mass_sto.append(m_sto[l+1]*(t1[l+1]-t1[l]) + mass_sto[l])

    Vc = Va0*10**-6 # Volumen de la camára en m^3
    A_t = At*10**-6 # Área de la garganta en m^3
    A_e = Ae*10**-6 # Área de la salida en m^3

    Po_gage =  np.array(Po_abs2) - P_a
    Po_max = max(Po_abs2) - P_a
    Po_final = (0.02/100)*Po_max + P_a

    t_final = - np.log(Po_final/(Po_abs1[-1]*10**-6))*((Vc)*ex_vel(C.R/M, To, k))/(C.R/M*To*(A_t))
    t2 = np.linspace(t1[-1], t1[-1]+ t_final, 100) # tiempo de taill off

    t_thrust = np.append(np.array(t1), t2[1:])

    A1 = (C.R/M)*To*(A_t)*(t2[1:]-t1[-1])
    A2 = (Vc)*ex_vel(C.R/M, To, k)

    Pc = np.array(Po_abs2[-1]*(np.exp(-A1/A2))) # presión de taill off

    Po_gage = np.array(Po_abs1)*10**-6 - P_a
    Presion_camara = np.append(Po_gage, Pc)

    # figura 1
    # %%
    F       = [0]
    C_f     = [N_noz]
    P_e     = []
    Ae_At   = [1]
    I_t     = [0]
    abc     = np.append(np.array(Po_abs1)*10**-6, Pc) # Es solo un paso intermedio
    Po_thrust = np.append(abc, 0)
    P       = lambda Po: Po*(1+(k-1)/2*Me**2)**(-k/(k-1))

    for j in range(len(t_thrust)-1):

        if P(Po_thrust[j]*10**6) < P_a*10**6:
            P_e.append(P_a*10**6)
        else :
            P_e.append(P(Po_thrust[j]*10**6))
        
        C_f.append(N_noz*Cf_(P_e[j], Po_thrust[j+1]*10**6, P_a*10**6, A_t, A_e, k))

        F.append(C_f[j]*Po_thrust[j]*10**6*A_t)

        I_t.append((F[j+1] + F[j])/2*(t_thrust[j+1]-t_thrust[j]))

        if Presion_camara[j] > 0 :
            Ae_At.append(A1A2(Presion_camara[j]*10**6, P_a*10**6, k))
    # figura 4

    # %%

    # Crear una función spline para F en función del tiempo t_thrust
    F_spline = CubicSpline(t_thrust, F)
    m_noz_spline = CubicSpline(t1, m_noz)


    #Graficar para verificar la interpolación

    t1_new = np.linspace(t1[0], t1[-1], 1000)
    m_noz_new = m_noz_spline(t1_new)


    t_new = np.linspace(min(t_thrust), max(t_thrust), 1000)
    F_new = F_spline(t_new)

    # figura 2
    # %%
    def F_(t, t_thrust):
        return np.piecewise(t, [t <= t_thrust, t > t_thrust], [F_spline(t), 0])

    def dm_(t, t_b):
        return np.piecewise(t, [t <= t_b, t > t_b], [m_noz_spline(t), 0])

    # Sistema de ecuaciones
    def Sis(CI, t, p):
        m, Cd, A, t_thrust, t_b, theta_0, fase = p
        x, u, y, v = CI

        F_t = F_(t, t_thrust)
        dm_t = dm_(t, t_b)
        
        dxdt = u
        dudt = (1/(m-fase*(dm_t)*t))*(fase*F_t - 0.5*densidad_aire(y)*A*Cd*(u**2 + v**2))*np.cos(theta_0)

        dydt = v
        dvdt = (1/(m-fase*(dm_t)*t))*(fase*F_t - 0.5*densidad_aire(y)*A*Cd*(u**2 + v**2))*np.sin(theta_0) - C.g
        
        return [dxdt, dudt, dydt, dvdt]


    t_thr = t_thrust[-1] # tiempo de empuje

    # Calculos fase 1
    CI1 = [0, 0, 0, 0] # Condiciones iniciales x0, Vx0, y0, Vy0
    t_1 = np.linspace(0, t_thr, 1000)
    p1 = [M_tot[i], Cd , Af, t_thr, t1[-1], theta_0, 1]
    Sol_1 = odeint(Sis, CI1, t_1, args=(p1, ))

    # Calculos fase 2
    
    CI2 = [Sol_1[-1,0], Sol_1[-1,1], Sol_1[-1,2], Sol_1[-1,3]] # Condiciones iniciales
    t_2 = np.linspace(t_thr,  50, 1000) # Si en la gráfica no se ve todo el descenso cambiar el tiempo final
    p2 = [m, Cd , Af, t_thr, t1[-1], theta_0, 0]
    Sol_2 = odeint(Sis, CI2, t_2, args=(p2, ))
    

    h.append(max(Sol_2[:,2]))

    if h[i+1] < h0:
        M_tot.append(M_tot[i] + Inc_M)
        m_p.append(M_tot[i+1]-m)
    
    print(f'Altura alcanzada = {h[i+1]} m')
    i+=1
    

    # Fd_max = 1/2*(p1[1]*p1[2]*densidad_aire(Sol_1[-1, 2])*(Sol_1[-1,3]**2+Sol_1[-1,1]**2)) # Fuerza de arrastre maxima
    # figura 3

# Figura 1
plt.figure(figsize=(7,7))
plt.plot(t_thrust[0:-1], Presion_camara*1000000/6895, 'b-', markersize='2.5')
plt.ylim(bottom=0)
plt.ylabel('Presión (psi)')
plt.xlabel('tiempo (s)')
plt.title('Presión de la camára')
plt.grid()

# Figura 4
plt.figure(figsize=(7,7))
plt.plot(t_thrust, F, 'b-')
plt.ylim(bottom=0)
plt.ylabel('Empuje (N)')
plt.xlabel('tiempo (s)')
plt.title('Gráfica de empuje')
plt.grid()

# Figura 2
fig1, ax = plt.subplots(2, 1, figsize=(7, 7))

ax[0].plot(t_thrust, F, 'bo', markersize='2.5', label='Datos originales')
ax[0].plot(t_new, F_new, 'r-', label='Interpolación Spline')
ax[0].legend()
ax[0].set_xlabel('Tiempo (s)')
ax[0].set_ylabel('Empuje (N)')
ax[0].set_title('Interpolación de F usando Spline Cúbico')
ax[0].set_ylim(bottom=0)

ax[1].plot(t1, m_noz, 'bo', markersize='2.5', label='Datos originales' )
ax[1].plot(t1_new, m_noz_new, 'r-')
ax[1].set_xlabel('Tiempo (s)')
ax[1].set_ylabel('Flujo masico (kg/s)')
ax[1].set_title(r'Interpolación de $\dot{m}$ usando Spline Cúbico')
ax[1].set_ylim(bottom=0)

# Figura 3
fig2, (ax1, ax2) = plt.subplots(2, figsize=(7, 7))

ax1.plot(t_1, Sol_1[:, 2], 'r-') # posición en y 
ax2.plot(t_1, Sol_1[:, 0], 'r-') # posición en x
ax1.plot(t_2, Sol_2[:, 2], 'b-') # posición en y
ax2.plot(t_2, Sol_2[:, 0], 'b-') # posición en x

ax1.set_xlabel('Tiempo [s]')
ax1.set_ylabel('y [m]')
ax1.set_ylim(bottom=0, top = max(Sol_2[:, 2]) + 5 )

ax2.set_xlabel('Tiempo [s]')
ax2.set_ylabel('x [m]')
ax2.set_ylim(bottom=0)


indice_y = np.where(np.diff(np.sign(Sol_2[:, 2])))
indice_x = np.where(Sol_2[:, 0] == max(Sol_2[:, 0]))
print(indice_x, max(Sol_2[:, 0]), Sol_2[0, 0] )
        

# Establecer el límite en el eje x en función del tiempo en el que y es igual a cero, si se encontró
ax1.set_xlim(left=0, right = t_2[indice_y[0][0]])
ax2.set_xlim(left=0, right = t_2[indice_x[0][0]])


# Datos
print('Datos relevantes')
print(f'Tiempo de quemado = {t1[-1]}', f'Tiempo de empuje = {t_thrust[-1]}')
print(f'Po max = {max(Po_abs2)*1000000/6895}', f'Po prom = {np.mean(Po_abs2)*1000000/6895}' )
print(f'F max = {max(F)}', f'F prom = {np.mean(F)}' )
print(f'It = {sum(I_t)}', f'Isp = {sum(I_t)/(C.g*m_p[-1])}',  f'm prop = {m_p[-1]}', f'm tot = {M_tot[-1]}')
print('Dimensiones de la camára y tobera')
print(f'L = {L0}', f' Dt = {Dt}', f' De = {De}',f'Dp = {Dp}')
print(f'Ap/At = {ApAt}', f' Ae/At = {AeAt}', f' Ae/At max = {1/min(Ae_At)}', f'Ae/At prom = {1/np.mean(Ae_At)}')
print(f'Kn = {Kn}', )
print(f'Apogeo = {max(h)}' )
#print(Ae_At[1], Presion_camara[0])



# Muestra todas las figuras
plt.tight_layout()
plt.show()


#Notas
# Modificando At/Ac y Ap/At se puede obtener una presion menor y un empuje menor igual
# Disminuyendo la masa del cohete sin propelente se reduce mucho la presion maxima y el empuje maximo

# No estoy muy seguro porque me da solo un poco los datos diferentes al excel de Nakka, creo que es porque se tiene que poner una presion deseada, 
# que realmente aqui no se hace eso, mas bien se impone un Ac/At y el Ap/At para así trabajar, ahora una pregunta es si se podria hacer igual con la presión.

# 14/10/2023
# El coeficiente de arrastre del cohete no solo es del cohete unicamente sino que tambien influyen las aletas y la
# el largo del cohete.

# 15/10/2023
# las graficas obtenidas no cuentan con el descenso con paracaidas es una edición que en el futuro se tienen que añadir.
# el angulo theta en el acenso debe cambiar con el tiempo