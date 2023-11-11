# Version: Beta
#   Esta versión no puede calcular correctamente los lanzamiento parabolicos solo cuando el angulo es igual a pi/2 (Sin inclinación)
#   

# @author: Leonel Gerardo González Pérez
# Este codigo esta basado en en el excel de Richard Nakka 'SRM_2013': https://www.nakka-rocketry.net/soft/SRM_2023.zip

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cte
import pandas as pd
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
    return np.sqrt((4*A)/np.pi)

def ex_vel(R, To, k):
    a = (k+1)/(k-1)
    c = np.sqrt(R*To/(k*(2/(k+1))**a))
    return c

def Ab_(Outer, Core, Ends, p):
    N, Do, do, Lo = p
    t_web = (Do - do) / 2
    x = np.linspace(0, t_web, 1000)
    D = Do - Outer * (2 * x)
    d = do + Core * (2 * x)
    L = Lo - Ends * (2 * x)

    A_o = N * np.pi * D * L if Outer == 1 else 0
    A_c = N * np.pi * d * L if Core == 1 else 0
    A_e = N * 1/2 * np.pi * (D**2 - d**2) if Ends == 1 else 0

    Ab_tot = [A_o, A_c, A_e, A_o + A_c + A_e]

    return Ab_tot

def Kn_(At, Ab):

    Kn = (Ab)/At
    Kno = min(Kn)
    Kn_prom = np.mean(Kn)
    Kn_max = max(Kn)

    return [Kno, Kn_prom, Kn_max]

def Ve(To, R, k, Pe, Po ):
    Ve = np.sqrt(2*To*R*(k/(k-1))*(1-(Pe/Po)**((k-1)/k)))
    return Ve
    
def mdot(At,Po, k, To, R):
    m = At*Po*np.sqrt((1/(R*To))*k*(2/(k+1))**((k+1)/(k-1)))
    return m


# Variación de la densidad del aire segun la altura.
def densidad_aire(h):
    # Definición de constantes
    rho_0 = 1.225  # kg/m^3
    g = cte.g  # m/s^2
    M = 0.02897  # kg/mol
    R = 8.3144598  # J/mol·K
    T0 = 288.15  # K

    rho = rho_0 * np.exp(-g * M * h / (R * T0))
    return rho

#  Diccionario con los valores para cada combinación de Propelente y Creador
Datos_propelente = {
    ('KNSU', 'UDEG SPACE'): {'k': 1.141607, 'M': 0.03573117, 'To_T': 1820.426, 'rho_T': 1.8892, 'a': 8.26, 'n': 0.319},
    ('KNSU', 'Nakka'): {'k': 1.133, 'M': 0.04202, 'To_T': 1720, 'rho_T': 1.889, 'a': 8.26, 'n': 0.319},
    ('KNDX', 'Nakka'): {'k': 1.131, 'M': 0.04242, 'To_T': 1710, 'rho_T': 1.879}, # Tiene distintos valores para a y n en distintos rangos de presión
    ('KNSB', 'Nakka'): {'k': 1.137, 'M': 0.0399, 'To_T': 1600, 'rho_T': 1.841, 'a': 5.130, 'n': 0.22},
    ('KNER', 'Nakka'): {'k': 1.140, 'M': 0.03858, 'To_T': 1608, 'rho_T': 1.820, 'a': 2.9, 'n': 0.4},
    ('KNMN', 'Nakka'): {'k': 1.136, 'M': 0.03983, 'To_T': 1616, 'rho_T': 1.854, 'a': 5.130, 'n': 0.22},
    ('KNFR', 'Nakka'): {'k': 1.131, 'M': 0.04242, 'To_T': 1710, 'rho_T': 1.942, 'a': 7.4, 'n': 0.250},
    ('KNPSB', 'Nakka'): {'k': 1.163, 'M': 0.03639, 'To_T': 1858, 'rho_T': 1.923, 'a': 6.5, 'n': 0.628},
    ('ACPC', 'MIT'): {'k': 1.21, 'M': 0.02367, 'To_T': 2821, 'rho_T': 1.680,'a': 3.237295, 'n': 0.327392,} 
}



# %%
alpha = 15*np.pi/180 # angulo de la divergencia de la tobera 

N_noz = 0.85 # # eficiencia de la tobera
N_com = 0.95 # eficiencia de combustion
N_div = 1/2*(1 + np.cos(alpha)) # factor de corrección por la divergencia
N_po = 0.95 # factor de corrección de la presion de la camara


# %%

# Nombre del proyecto
Nombre = 'Cohete chikito'
formato ='excel' # Posibles: 'csv', 'txt', 'excel'
# %%
Propelente = 'KNSU' # Posibles: KNSU, KNDX, KNSB, KNER, KNMN, KNFR, KNPSB, ACPC
Creador = 'Nakka' # Posibles: Nakka, UDEG SPACE, MIT
Altura = 1538 # m, Altura sobre el nivel del mar
# %%

# Obtener los valores correspondientes
Prop = Datos_propelente.get((Propelente, Creador))

if Prop:
    k = Prop['k']         # relación de calores especificos
    M = Prop['M']         # kg/mol masa molecular de los gases
    To_T = Prop['To_T']   # K temperatura de combustion teorica
    print(To_T)
    rho_T = Prop['rho_T'] # g/cm^3 densidad teorica
    a = Prop['a']
    n = Prop['n']

    print(f"Valores para Propelente: {Propelente}, Creador: {Creador}")
    print(f"k = {k}, M = {M}, To_T = {To_T}, rho_T = {rho_T}, a = {a}, n = {n}" )

rho_rat = 0.9
To = To_T*N_com # Temperatura "real"
rho = rho_T*rho_rat    # g/cm^3 densidad medida

# %%
P0 = 250 #psi, presion de la camara objetivo
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
Cd = 0.6 # coeficiente de arrastre

#102.26
Dc = 1*25.4 #mm Diametro de la camara 
Df = 2 # in, diametro fuselaje

Dp =  (8/16)*25.4 # mm, de pulgadas a mm [(in)*(mm/in)] # Diametro del port
In = 0 # ancho de los inhibidores y/o revestimiento del grano

Dt = (4/16)*25.4 # mm, de pulgadas a mm [(in)*(mm/in)]

m = .4 # masa del cohete sin propelente


Af = np.pi*(Df/(2*39.37))**2 # área transversal (6 inch es el diametro del tubo de fuselaje)
Ac = np.pi*(Dc/2)**2 #mm**2 # área de la camara
At = np.pi*(Dt/2)**2
AcAt = Ac/At
ApAt = np.pi*(Dp/2)**2/At
# 1.6, 1.75


AeAt = RatArea(Me, Mt, k)
Ae = At*AeAt
De = diametro(Ae)

# %%

# Variables del ciclo iterativo


Datos = 10000 # cantidad de datos para la regresión del grano
h0 = 300 # altura deseada 


# Superficies de quemado

# Si es 1 significa que se toma en cuenta si es 0 no se toma en cuenta, si se pone algun otro número se tendran resultados erroneos
Bs = [0, 1, 1] # Superficies: exterior, nucleo, caras
N = 1 # Número de segmentos del propelente


Cf = Cf_(Pe, P0, Pa, At, Ae, k)*N_noz
c_ast = ex_vel(cte.R/M, To, k)
Ve = c_ast*Cf # Velocidad de salida de los gases de la tobera


mp0 = m*(np.exp(np.sqrt(2*cte.g*h0)/Ve)-1) # masa minima de propelente
m_p = [mp0] # se utiliza la formula del método que no se considera la fuerza de arrastre para la primera iteración
M_tot = [m + mp0]



# Ciclo iterativo
 
Tolerancia = 0.1 # Error de la altura alcanzada en metros

i = 0
h = [0]
while np.abs(h0-h[i]) > Tolerancia:
#while h[i] < h0:
    print(f'Iteración = {i + 1}')


    # Área de quemado
    V0 = m_p[i]*1000/rho # cm^3 Volumen con densidad experiemental del grano
    Va0 = (V0)/(1-ApAt*At*4/(np.pi*(Dc - 2*In)**2)) # cm^3 Volumen disponible de la camara
    
    Lc = Va0*10**3/Ac # Longitud de la cámara
    L0 = Va0*10**3/(np.pi*((Dc - 2*In)/2)**2) # mm  Longitud del grano con densidad experimental

    p0 = [N, Dc - 2*In, Dp, L0/N] #  N, D, d, L, At, = p
    Ab0 = Ab_(Bs[0],Bs[1],Bs[2], p0) # área de quemado
    Kn = Kn_(At, Ab0[-1])

    # %%
    
    t_web0  = np.array([(Dc-2*In-Dp)/2])
    Xinc    = np.linspace(0, t_web0[0]/(1+Bs[0]*Bs[1]), Datos)
    t_web   = np.append(t_web0, t_web0[0]-Xinc)

    D = Dc - 2*In - Bs[0]*2*Xinc
    d = Dp + Bs[1]*2*Xinc
    L = L0 - Bs[2]*N*2*Xinc

    A_t = At*10**-6 # Área de la garganta en m^3
    A_e = Ae*10**-6 # Área de la salida en m^3

    Vc   = Va0*10**-6 # m^3 # Volumen de la cámara
    V_G  = 1/4*np.pi*(D**2-d**2)*L/(1000**3) # m^3 # Volumen del grano dada la regresión
    V_F  = Vc - V_G # m^3

    P_a          = 0.101325 # Mpa Presión atmosferica

    Po_abs1     = [] # En Pa
    Po_abs2     = [P_a + 0.00001] # En Mpa
    m_grain     = rho*1000*V_G # kg
    m_gen       = [0]
    m_noz       = [0]
    m_sto       = [0]
    mass_sto    = [0]
    rho_prod    = []
    t1          = [0]
    R           = [a*(Po_abs2[0]**n)] # tasa de quemado

    for l in range(len(Xinc)-1):
        rho_prod.append(mass_sto[l]/V_F[l])

        Po_abs1.append(rho_prod[l]*cte.R/M*To+P_a*10**6)
        Po_abs2.append(Po_abs1[l]/(10**6))
        R.append(a*(Po_abs2[l])**n)
        t1.append(Xinc[1]/R[l+1] + t1[l]) # tiempo de quemado 
        m_gen.append((m_grain[l]-m_grain[l+1])/(t1[l+1]-t1[l]))
        
        if m_gen[l+1] < mdot(A_t,(Po_abs2[l+1]-P_a)*10**6, k,To,cte.R/M):
            if Po_abs2[l+1] > 0:
                m_noz.append( mdot(A_t,(Po_abs2[l+1]-P_a)*10**6, k,To,cte.R/M))
            else :
                m_noz.append(0)
        else :
            m_noz.append(mdot(A_t,(Po_abs2[l+1]-P_a)*10**6, k,To,cte.R/M))
        
        m_sto.append(m_gen[l+1]-m_noz[l+1])
        

        if m_sto[l+1]*(t1[l+1]-t1[l]) + mass_sto[l] < 0 :
            mass_sto.append(0)
        else :
            mass_sto.append(m_sto[l+1]*(t1[l+1]- t1[l]) + mass_sto[l])


    Po_gage =  np.array(Po_abs2) - P_a
    Po_max = max(Po_abs2) - P_a
    Po_final = (0.02/100)*Po_max + P_a # Presión final del taill off


    A1 = (cte.R/M)*To*(A_t)
    A2 = (Vc)*c_ast

    t_final = - np.log(Po_final/Po_abs2[-1])*A2/A1
    t2 = np.linspace(0, t_final, 100) # tiempo de taill off
    t_thrust = np.append(np.array(t1), t2[1:]+ t1[-1])

    Pc = np.array((Po_abs2[-1])*(np.exp(-A1*t2[1:]/A2))) # presión de taill off
    Po_gage = np.append(np.array(Po_abs1)*10**-6, Pc) - P_a
    Presion_camara = Po_gage 

    # %%
    F       = [0]
    C_f     = [N_noz]
    P_e     = []
    Ae_At   = [1]
    I_t     = [0]
    Po_thrust = np.append(Po_gage + P_a, 0.01) # presion de taill off
    P       = lambda Po: Po/((1+(k-1)/2*Me**2)**(k/(k-1))) # Funcion definida para hacer mas simple los calculos

    for j in range(len(t_thrust)-1):

        if P(Po_thrust[j]*10**6) < P_a*10**6:
            P_e.append(P_a*10**6)
        else :
            P_e.append(P(Po_thrust[j]*10**6))
        
        C_f.append(N_noz*Cf_(P_e[j], Po_thrust[j+1]*10**6, P_a*10**6, A_t, A_e, k))

        F.append(C_f[j]*Po_thrust[j]*10**6*A_t)

        I_t.append((F[j+1] + F[j])/2*(t_thrust[j+1]-t_thrust[j]))

        if Presion_camara[j] > 0 and P_a/Presion_camara[j] < 1 : # Aqui hay un poblema que resolver, no deberia tener las condiciones
            Ae_At.append(1/(A1A2(Presion_camara[j], P_a, k)))

    # %%

    # Crear una función spline para F en función del tiempo t_thrust
    F_spline = CubicSpline(t_thrust, F)
    m_noz_spline = CubicSpline(t1, m_noz)
    

    #Graficar para verificar la interpolación

    t1_new = np.linspace(t1[0], t1[-1], 1000)
    m_noz_new = m_noz_spline(t1_new)


    t_new = np.linspace(min(t_thrust), max(t_thrust), 1000)
    F_new = F_spline(t_new)

    # %%
    def F_(t, t_thrust):
        return np.piecewise(t, [t <= t_thrust, t > t_thrust], [F_spline(t), 0])

    def dm_(t, t_b):
        return np.piecewise(t, [t <= t_b, t > t_b], [m_noz_spline(t), 0])
    
    #def theta_(t, theta_0, Vy, Vx):
    #    return np.piecewise(t, [t == 0, t > 0], [theta_0, np.arctan(Vy/Vx)])

    # Sistema de ecuaciones
    def Sistema(CI, t, p):
        m, Cd, A, t_thrust, t_b, theta_0, fase = p
        x, u, y, v = CI

        F_t = F_(t, t_thrust)
        dm_t = dm_(t, t_b)
        
        dxdt = u
        dudt = (1/(m-fase*(dm_t)*t))*(fase*F_t - 0.5*densidad_aire(y+Altura)*A*Cd*(u**2 + v**2))*np.cos(theta_0)

        dydt = v
        dvdt = (1/(m-fase*(dm_t)*t))*(fase*F_t - 0.5*densidad_aire(y+Altura)*A*Cd*(u**2 + v**2))*np.sin(theta_0) - cte.g
        
        return [dxdt, dudt, dydt, dvdt]


    t_thr = t_thrust[-1] # tiempo de empuje

    # Calculos fase 1
    CI1 = [0, 0, 0, 0] # Condiciones iniciales x0, Vx0, y0, Vy0
    t_1 = np.linspace(0, t_thr, 1000)
    p1 = [M_tot[i], Cd , Af, t_thr, t1[-1], theta_0, 1]
    Sol_1 = odeint(Sistema, CI1, t_1, args=(p1, ))

    # Calculos fase 2
    
    CI2 = [Sol_1[-1,0], Sol_1[-1,1], Sol_1[-1,2], Sol_1[-1,3]] # Condiciones iniciales
    t_2 = np.linspace(t_thr,  50, 1000) # Si en la gráfica no se ve todo el descenso cambiar el tiempo final
    p2 = [m, Cd , Af, t_thr, t1[-1], theta_0, 0]
    Sol_2 = odeint(Sistema, CI2, t_2, args=(p2, ))
    

    h.append(max(Sol_2[:,2]))

    if np.abs(h0-h[i]) > Tolerancia:
        Inc_M = m_p[i]/h[i+1]*(h0-h[i+1])

        M_tot.append(M_tot[i] + Inc_M)
        m_p.append(M_tot[i+1]-m)
    
    print(f'Altura alcanzada = {h[i+1]} m')
    i+=1
  
# It2 = integrate.quadrature(F_spline, t_thrust[0], t_thrust[-1], maxiter=900) # Sirve para poder tener una referencia con el It calculado

# Figura 1
plt.figure(figsize=(7,7))
plt.plot(t_thrust[:-1], Presion_camara*1000000/6895, 'b-', markersize='2.5')
#plt.ylim(bottom=0)
plt.ylabel('Presión (psi)')
plt.xlabel('Tiempo (s)')
plt.title('Presión de la cámara')
plt.grid()

# Figura 2
plt.figure(figsize=(7,7))
plt.plot(t_thrust, F, 'b-')
#plt.ylim(bottom=0)
plt.ylabel('Empuje (N)')
plt.xlabel('Tiempo (s)')
plt.title('Gráfica de empuje')
plt.grid()

# Figura 3
fig1, ax = plt.subplots(2, 1, figsize=(7, 7))

ax[0].plot(t_thrust, F, 'bo', markersize='2.5', label='Datos originales')
ax[0].plot(t_new, F_new, 'r-', label='Interpolación Spline')
ax[0].legend()
ax[0].set_xlabel('Tiempo (s)')
ax[0].set_ylabel('Empuje (N)')
ax[0].set_title('Interpolación de F usando Spline Cúbico')
ax[0].set_ylim(bottom =0 )
ax[0].grid()

ax[1].plot(t1, m_noz, 'bo', markersize='2.5', label='Datos originales' )
ax[1].plot(t1_new, m_noz_new, 'r-')
ax[1].set_xlabel('Tiempo (s)')
ax[1].set_ylabel('Flujo masico (kg/s)')
ax[1].set_title(r'Interpolación de $\dot{m}$ usando Spline Cúbico')
ax[1].set_ylim(bottom=0)
ax[1].grid()

# Figura 4
plt.figure(figsize=(7,7))

Ymax = np.where(Sol_2[:, 2] == max(Sol_2[:, 2])) # indice donde la altura es maxima

plt.plot(t_1, Sol_1[:, 2], 'r-', label = 'Fase propulsada' ) # posición en y 
plt.plot(t_2[:Ymax[0][0]], Sol_2[:Ymax[0][0] , 2], 'b-', label = 'Fase sin propulsción') # posición en y

plt.xlabel('Tiempo (s)')
plt.ylabel('y (m)')


plt.xlim(left=0, right = t_2[Ymax[0]])
plt.ylim(bottom=0, top = max(Sol_2[:, 2]) + 5 )
plt.legend()
plt.grid()


## Figura 5: Velocidad

fig3, (ax1, ax2) = plt.subplots(2, figsize=(7, 7))

ax1.plot(t_1, Sol_1[:, 3],'r-', markersize='2', label = 'Fase propulsada') # Velocidad en y 
ax2.plot(t_1, Sol_1[:, 1],'r-', markersize='2', label = 'Fase propulsada') # Velocidad en x
ax1.plot(t_2[:Ymax[0][0]], Sol_2[:Ymax[0][0], 3], 'b-', label = 'Fase sin propulsción') # Velocidad en y
ax2.plot(t_2[:Ymax[0][0]], Sol_2[:Ymax[0][0], 1], 'b-', label = 'Fase sin propulsción') # Velocidad en x

ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel(r'$V_{y}$ (m)')
#ax1.set_ylim(bottom=0, top = max(Sol_2[:, 3]) + 5 )
ax1.legend()
ax1.grid()

ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel(r'$V_{x}$ (m)')
ax2.legend()
ax2.grid()

# Figura 6

plt.figure(figsize=(7,7))
plt.plot(t_thrust[0:-1], np.ones(len(t_thrust[0:-1]))*Pa,'b-', label = r'$P_{a}$')
plt.plot(t_thrust[0:-1], P(Presion_camara*1000000/6895 + Pa), color='#FFA500', label = r'$P_{o}$')
plt.plot(t_thrust[0:-1], np.array(P_e)*(6895)**-1, 'k--', linewidth='2', label = r'$P_{e}$')
plt.title(r'Relación de presiones con $P_{e}$')
plt.ylabel('Tiempo (s)')
plt.ylabel('Presión (psi)')
plt.legend()
plt.grid()

t_Brn = t1[-1]
# t_thr = t_thr
Po_max = (Po_max + P_a)*1000000/6895
Po_prom = np.mean(Po_abs2)*1000000/6895
F_max = max(F)
F_prom = np.mean(F)
It = sum(I_t)
Isp = sum(I_t)/(cte.g*m_p[-1])
mp = m_p[-1]
m_tot = M_tot[-1]


# Datos
print('')
print('Datos relevantes')
print('')
print(f'Tiempo de quemado = {t_Brn}', f'Tiempo de empuje = {t_thr}')
print(f'Po max = {Po_max}', f'Po prom = {Po_prom}')
print(f'F max = {F_max}', f'F prom = {F_prom}' )
print(f'It = {It}', f'Isp = {Isp}')
print(f'm = {m}', f'm prop = {mp}', f'm tot = {m_tot}')
#print(f'It2 = {It2}' )
print(f'Me = {Me}', )

print('')
print('Dimensiones de la cámara y tobera')
print('')
print(f'L = {Lc}', f' Dt = {Dt}', f' De = {De}', f'Dp = {Dp}')
print(f'Ap/At = {ApAt}',f'Ac/At = {AcAt}')
print(f'Ae/At a Po(obj) = {AeAt}', f' Ae/At max = {np.nanmax(Ae_At)}', f'Ae/At prom = {np.nanmean(Ae_At)}')

data = {
    'Altura objetivo' : [h0],
    'Tiempo de quemado (s)' : [t_Brn],
    'Tiempo de empuje (s)' : [t_thr],
    'Presión Objetivo (psi)' : [P0],
    'Presión máxima (psi)' : [Po_max],
    'Presión promedio (psi)' : [Po_prom],
    'Empuje máximo (N)' : [F_max],
    'Empuje promedio (N)' : [F_prom],
    'Impulso total (Ns)' : [It],
    'Impulso especifico (s)' : [Isp],
    'Masa del cohete (Kg)' : [m],
    'Masa propelente (Kg)' : [mp],
    'Masa total (Kg)' : [m_tot],
    'Longitud de la cámara (mm)' : [Lc],
    'Diametro cámara (mm)' : [Dc],
    'Diametro garganta (mm)' : [Dt], 
    'Diametro salida (mm)' : [De],
    'Ac/At' : [AcAt],
    'Ae/At, Po(obj)' : [AeAt],
    'Ae/At max' : [np.nanmax(Ae_At)],
    'Ae/At prom' : [np.nanmean(Ae_At)],
    'Ap/At' : [ApAt],
    'Numero de segmentos' : [N],
    'Diametro del grano (mm)' : [Dc-2*In],
    'Diametro nucleo del grano (mm)' : [Dp],
    'Relación de calores especificos' : [k],
    'Peso molecular (kg/mol)' : [M],
    'Temperatura de combustión': [To],
    'Densidad teorica (g/cm^3)' : [rho_T],
    'Densidad medida (g/cm^3)' : [rho],
    'Coeficiente de presión (Mpa^1/n)' : [a],
    'Exponente de presión ' : [n],
    'Angulo de lanzamiento ' : [theta_0*180/np.pi],
    'Diametro fuselaje (mm)' : [Df],
    'Coeficiente de arrastre cohete' : [Cd],
}

df1 = pd.DataFrame(data)


if formato == 'csv':
    df1.to_csv(f'{Nombre}.csv', index = False)
elif formato == 'txt':
    df1.to_csv(f'{Nombre}.txt', sep = '\t', index = False)
elif formato == 'excel':
    df1.to_excel(f'{Nombre}.xlsx' , index = False)
else:
    print("Formato no compatible")


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

# 28/10/23
# El rendimiento de el vuelo y las medidas en general tambien pueden estar afectadas por la presion deseada o P0 en el codigo, esto
# para que se tenga en cuenta en futuras referencias
#
# Si la tobera se sella por dentro de la cámara se puede poner un angulo de convergencia mayor para el flujo de los gases y con eso se deben
# de hacer cambios al codigo, pero seria una mejora, ya que puede reducir el peso de la tobera 
#
#
# 04/11/2023
# Una adición podria ser obtener un archivo .csv donde muestre todos los datos de obtenidos en los ultimos 'print' además de su conversión 
# en pulgadas
#
# 06/11/2023
# Se le acaba de hacer una adición para poder utilizar otros propelentes de manera mas simple sin tener que andar canbiando tediosamente los
# datos del propelente a mano. Ademas se le tiene que hacer una modificación al codigo de igual manera para los coeficientes de presión
# debido a que para cada uno es diferente
#
# 07/11/2023
# Se le tiene que añadir la altura barometrica y la presion barometrica para que a cualquier altura se pueda obtener los datos necesarios, además de agregarle la 
# velocidad del viento
#
# 10/11/2023
# Problema de la grafica de presión de la cámara en tail off arreglado
#
#