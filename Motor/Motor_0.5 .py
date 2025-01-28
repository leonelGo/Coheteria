# Version: Beta
#   Esta versión no puede calcular correctamente los lanzamiento parabolicos solo cuando el angulo es igual a pi/2 (Sin inclinación)
#   

# @author: Leonel Gerardo González Pérez
# Este codigo esta basado en en el excel de Richard Nakka 'SRM_2013': https://www.nakka-rocketry.net/soft/SRM_2023.zip

# %%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as cte
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from scipy.integrate import quad
import time
from scipy.optimize import fsolve, minimize_scalar
from scipy.stats import linregress
from scipy.optimize import curve_fit
import os
from scipy.integrate import simps
#import scipy.integrate as integrate

start_time = time.time()

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
############################################################
#def Cf_(Pe, P0, Pa, At, Ae, k):
#    a = (k+1)/(k-1)
#    return np.sqrt((2*k**2/(k-1)*(2/(k+1))**a)*(1-(Pe/P0)**((k-1)/k))) + (Pe - Pa)*Ae/(P0*At)

def Cf_(Pe, P0, Pa, At, Ae, k):
    # Convertir a numpy array solo si la entrada es una lista
    if isinstance(Pe, list):
        Pe = np.array(Pe)
    if isinstance(P0, list):
        P0 = np.array(P0)
    if isinstance(Pa, list):
        Pa = np.array(Pa)
    if isinstance(At, list):
        At = np.array(Pa)

    term1 = (2 * k**2 / (k - 1) * (2 / (k + 1))**((k + 1) / (k - 1)))
    term2 = (Pe / P0)**((k - 1) / k)
    term2 = np.minimum(term2, 1)  # Asegurar que term2 sea siempre <= 1 para evitar raíces cuadradas de negativos
    term3 = np.sqrt(term1 * (1 - term2))
    term4 = (Pe - Pa) * Ae / (P0 * At)

    Cf = term3 + term4

    # Condición para Pe < Pa
    if isinstance(Pe, np.ndarray):
        condition = Pe < Pa - P0*term3*At/Ae
        Cf[condition] = N_throat * N_div * N_noz * N_skin * 0
    else:
        if Pe < Pa - P0*term3*At/Ae:
            Cf = N_throat * N_div * N_noz * N_skin * 0

    return Cf

########################################################
def diametro(A):
    return np.sqrt((4*A)/np.pi)

def ex_vel(R, To, k):
    a = (k+1)/(k-1)
    c = np.sqrt(R*To/(k*(2/(k+1))**a))
    return c

def Ab_(Outer, Core, Ends, p):
    N, Do, do, Lo = p
    t_web = (Do - do) / 2
    x = np.linspace(0, t_web, 2000)
    D = Do - Outer * (2 * x)
    d = do + Core * (2 * x)
    L = Lo - Ends * (2 * x)

    A_o = N * np.pi * D * L if Outer == 1 else 0
    A_c = N * np.pi * d * L if Core == 1 else 0
    A_e = N * 1/2 * np.pi * (D**2 - d**2) if Ends == 1 else 0

    Ab_tot = [A_o, A_c, A_e, A_o + A_c + A_e]

    return Ab_tot


def Ve(To, R, k, Pe, Po ):
    Ve = np.sqrt(2*To*R*(k/(k-1))*(1-(Pe/Po)**((k-1)/k)))
    return Ve
    
def mdot(At,Po, k, To, R):
    m = At*Po*np.sqrt((1/(R*To))*k*(2/(k+1))**((k+1)/(k-1)))
    return m

def Po_aprox(Kn,k, R, T, rho, a, n):
    C = np.sqrt(R*T/(k*(2/(k+1))**((k+1)/(k-1))))
    Po =((Kn*a*rho*C)*10**(-6))**(1/(1-n))
    return Po*145.03

def Ab_(Dg, In, Dp, Lg0, x, N, Bs):
    D = Dg - 2 * In - Bs[0] * 2 * x # m
    d = Dp + Bs[1] * 2 * x          # m
    L = Lg0 - Bs[2] * N * 2 * x     # m

    Ab_e = Bs[2] * 2 * N * np.pi * (D**2 - d**2) / 4
    Ab_c = Bs[1] * np.pi * d * L
    Ab_s = Bs[0] * np.pi * D * L


    if type(x) == np.ndarray:
        Ab = Ab_e + Ab_c + Ab_s
    elif x <= (Dg - Dp) / 2 and x > 0:
        Ab = Ab_e + Ab_c + Ab_s # mm^2
    else:
        Ab = 0

    return Ab


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

def Propiedades_Atm(h):
    P0 = 101325  # Presión a nivel del mar en Pascales
    T0 = 288.15  # Temperatura a nivel del mar en Kelvin
    g = 9.80665  # Aceleración debido a la gravedad en m/s^2 (valor aproximado)
    M = 0.0289644  # Masa molar del aire en kg/mol
    R = 8.31447  # Constante de los gases ideales en J/(mol·K)

    # Capas de la atmósfera con sus alturas y parámetros
    layers = [
        (0, 288.15, 101325, -0.0065),
        (11000, 216.65, 22632.1, 0),
        (20000, 216.65, 5474.89, 0.001),
        (32000, 228.65, 868.02, 0.0028),
        (47000, 270.65, 110.91, 0),
        (51000, 270.65, 66.94, -0.0028),
        (71000, 214.65, 3.956, -0.0028),
        (84852, 186.87, 0.3734, 0),
        (95000, 186.87, 0.0971, 0.002),
        (105000, 190.65, 0.0135, 0)
    ]
    
    # Comprobar si h es escalar
    if np.isscalar(h):
        h = np.array([h])
        scalar_input = True
    else:
        h = np.array(h)
        scalar_input = False
    
    P = np.zeros_like(h, dtype=float)
    T = np.zeros_like(h, dtype=float)
    rho = np.zeros_like(h, dtype=float)
    
    for i, (h_b, T_b, P_b, L) in enumerate(layers):
        if i == 0:
            mask = h < layers[1][0]
        else:
            mask = (h >= layers[i][0]) & (h < layers[i + 1][0]) if i + 1 < len(layers) else (h >= layers[i][0])
        
        if L == 0:
            T[mask] = T_b
            P[mask] = P_b * np.exp(-g * M * (h[mask] - h_b) / (R * T_b))
        else:
            T[mask] = T_b + L * (h[mask] - h_b)
            P[mask] = P_b * (T[mask] / T_b) ** (-g * M / (R * L))
    
    rho = P / (R / M * T)
    
    if scalar_input:
        return P[0], rho[0], T[0]
    return P, rho, T

def coeficiente_empuje_optimo(k, Me):
    # Definimos la parte de la ecuación antes de la raíz
    term1 = k**2 * (2 / (k + 1))**((k + 1) / (k - 1))
    
    # Calculamos el denominador del segundo término dentro de la raíz
    denominator = 1 + ((k - 1) / 2) * Me**2
    
    # Calculamos el término dentro de la raíz
    term_inside_sqrt = term1 * (Me**2 / denominator)
    
    # Calculamos el primer término del coeficiente de empuje óptimo (dentro de la raíz)
    sqrt_term = np.sqrt(term_inside_sqrt)
    
    
    # Sumar ambos términos
    C_F = sqrt_term 
    
    return C_F

def Cf2_(pe, po, pa, At, Ae, k):
    term1 = np.sqrt((2 * k**2 / (k - 1)) * (2 / (k + 1))**((k + 1) / (k - 1)) * (1 - (pe / po)**((k - 1) / k)))
    term2 = ((pe - pa) / po) * (Ae / At)
    return term1 + term2 



def write_eng_file(DATA, t_thrust, F):
    filename, motor_name, motor_diameter, motor_length, delays, propellant_weight, total_weight, manufacturer = DATA
    with open(filename, 'w') as file:
        # Escribir los comentarios y encabezado
        file.write(f"; {motor_name} data generated\n")
        file.write(f"{motor_name} {motor_diameter} {motor_length} {delays} {propellant_weight:.6f} {total_weight:.6f} {manufacturer}\n")
        
        # Escribir los datos de tiempo y empuje
        for t, thrust in zip(t_thrust, F):
            file.write(f"{t:.3f} {thrust:.3f}\n")
        
        # Asegurarse de que el último punto tenga empuje cero
        if F[-1] != 0:
            file.write(f"{t_thrust[-1]:.3f} 0.000\n")
        file.write(";")

def P_Po(Po, Me, k):
    P       =  Po/((1+(k-1)/2*Me**2)**(k/(k-1)))
    return P

def Ec_Presion(t, CI, p):
    global fase2_activa, tb, Guardado
    a, n, Dc, Dp, Lg0, Bs, In, rho_p, k, M, To, At, N = p
    Po, rho_0, v0, m, mg, x = CI
    
    PresionAtm = Propiedades_Atm(Altura)[0]
    print(x*1e3)
    if x <= (Dc - Dp) / 2 and not fase2_activa == True:
        r = a * (Po + PresionAtm )**n
        Ab = Ab_(Dc, In, Dp, Lg0, x, N, Bs)

    else:
        if x >= (Dc - Dp) / 2 and not fase2_activa == True:
                tb = t
                Guardado = a * (Po + Propiedades_Atm(Altura)[0])**n
                print(f'Valor de la regresión {x*1e3} mm al tiempo tb = {tb} s')
                fase2_activa = True

                if x <= (Dc - Dp) / 2:
                    raise ValueError(f"La regresión lineal es menor que {(Dc-Dp)/2}.")
           
        Ab = 0   
        r = 0
    


    dxdt = r
    dmdt = At * (Po + PresionAtm) * np.sqrt((1 / (cte.R / M * To)) * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
    dmgdt = Ab * r * rho_p
    dP0dt = (cte.R * To / (M * v0)) * (dmgdt - Ab * r * rho_0 - dmdt)
    dv0dt = Ab * r
    drho_0dt = 1 / (cte.R / M * To) * dP0dt

    return [dP0dt, drho_0dt, dv0dt, dmdt, dmgdt, dxdt]


def Ec_Presion2(t, CI, p):
        global fase2_activa, tb
        m1, Cd, A, theta_0,  p_P = p
        a, n, Dc, Dp, Lg0, Bs, In, rho_p, k, M, To, At, Ae, Me, N = p_P

        x, u, y, v, Po, rho_0, v0, m2, mg, x_web = CI
        #print(x)
        
        if PresionCambiante == True:
            P_atm, Rho_atm, _ = Propiedades_Atm(y+Altura)
            #print(f'Pa camabiante {P_atm/6894.76}')
        else: 
            P_atm, Rho_atm, _ = Propiedades_Atm(Altura)
            #print(f'Pa no camabiante {P_atm/6894.76}')

        if x_web <= (Dc - Dp) / 2 and not fase2_activa == True:
            r = a * (Po + P_atm)**n
            
        else:
            if x_web >= (Dc - Dp) / 2 and not fase2_activa == True:
                    tb = t

                    print(f'Valor de la regresión {x*1e3} mm al tiempo tb = {tb} s')
                    fase2_activa = True
                    if x_web <= (Dc - Dp) / 2:
                        raise ValueError(f"La regresión lineal es menor que {(Dc-Dp)/2}.")
            r = 0
        Ab = Ab_(Dc, In, Dp, Lg0, x_web, N, Bs)


        dx_webdt = r
        dmdt = At * (Po + P_atm) * np.sqrt((1 / (cte.R / M * To)) * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        dmgdt = Ab * r * rho_p
        dP0dt = (cte.R * To / (M * v0)) * (dmgdt - Ab * r * rho_0 - dmdt)
        dv0dt = Ab * r
        drho_0dt = 1 / (cte.R / M * To) * dP0dt

        P_e = P_Po(Po + P_atm , Me, k)
        C_f = N_throat*N_div*N_noz*(N_skin*Cf_(P_e, (Po + P_atm), P_atm, At, Ae, k)+(1-N_skin)) # agregar las eficiencias a los parametros

        #print(C_f)

        #print(C_f)
        F = C_f*At*(Po + P_atm)
        dxdt = u
        dudt = (1/(m1-(dmdt)*t))*(F 
                                    - 0.5*(Rho_atm)*A*Cd*(u**2 + v**2))*np.cos(theta_0)

        dydt = v
        dvdt = (1/(m1-(dmdt)*t))*(F 
                                    - 0.5*(Rho_atm)*A*Cd*(u**2 + v**2))*np.sin(theta_0) - cte.g

        return [dxdt, dudt, dydt, dvdt, dP0dt, drho_0dt, dv0dt, dmdt, dmgdt, dx_webdt]


    # Sistema de ecuaciones
def Sistema(t, CI, p):
    m, Cd, A, t_thrust, t_b, theta_0, fase = p
    x, u, y, v = CI

    F_t = F_(t, t_thrust)
        #print(f'F {F_t}')
    dm_t = dm_(t, t_b)

    dxdt = u
    dudt = (1/(m-fase*(dm_t)*t))*(fase*F_t 
                                    - 0.5*(Propiedades_Atm(y+Altura)[1])*A*Cd*(u**2 + v**2))*np.cos(theta_0)

    dydt = v
    dvdt = (1/(m-fase*(dm_t)*t))*(fase*F_t 
                                    - 0.5*(Propiedades_Atm(y+Altura)[1])*A*Cd*(u**2 + v**2))*np.sin(theta_0) - cte.g

    return [dxdt, dudt, dydt, dvdt]

    # Evento para detener la integración cuando y = 0
def evento_y_cero(t, CI, p):
    return CI[2]  # y es el tercer componente de CI

    # El evento debe tener direction=-1 para detectar cuando y cruza cero desde valores positivos
evento_y_cero.terminal = True
evento_y_cero.direction = -1


def event_pressure_zero(t, CI, p):
    Po, rho_0, v0, m, mg, X = CI
    return Po

event_pressure_zero.terminal = True
event_pressure_zero.direction = -1

def event_pressure_zero2(t, CI, p):
    x, u, y, v, Po, rho_0, v0, m, mg, x = CI
    return Po

event_pressure_zero2.terminal = True
event_pressure_zero2.direction = -1

def ClasificacionMotores(impulso_total):
    clases = [
        (1.26, 2.5, "A"), (2.5, 5, "B"), (5, 10, "C"), (10, 20, "D"), 
        (20, 40, "E"), (40, 80, "F"), (80, 160, "G"), (160, 320, "H"), 
        (320, 640, "I"), (640, 1280, "J"), (1280, 2560, "K"), (2560, 5120, "L"), 
        (5120, 10240, "M"), (10240, 20480, "N"), (20480, 40960, "O"), 
        (40960, 81920, "P"), (81920, 163840, "Q"), (163840, 327680, "R"), 
        (327680, 655360, "S"), (655360, 1310720, "T"), (1310720, 2621440, "U"), 
        (2621440, 5242880, "V")
    ]
    
    return next((clase for min_impulso, max_impulso, clase in clases 
                 if min_impulso <= impulso_total <= max_impulso), 
                "NE")


def generar_archivo_eng(manofacturador, mision, diametro, longitud, masa_propelente, masa_motor, tiempos, empuje, ImpulsoTotal):
    Clase = ClasificacionMotores(ImpulsoTotal)
    nombre_motor = Clase + str(round(np.mean(empuje)))
    with open(f'{manofacturador}_{nombre_motor}_{mision}.eng', 'w') as file:
        
        # Escribir la cabecera
        file.write(f"; Archivo generado para el motor {nombre_motor}\n")
        file.write(f"{nombre_motor} {diametro} {longitud} P {masa_propelente} {masa_motor} {manofacturador}\n")
        
        # Escribir los datos de tiempo y empuje
        for t, F in zip(tiempos, empuje):
            file.write(f"{t} {F}\n")
        
        # Escribir el final del archivo
        file.write("0.0 0.0\n")

#  Diccionario con los valores para cada combinación de Propelente y Creador
Datos_propelente = {
    ('KNSU', 'UDEG SPACE'): {'k': 1.141607, 'M': 0.03573117, 'To_T': 1820.426, 
                             'rho_T': 1.8892, 'a': 8.26, 'n': 0.319},
    ('KNSU', 'Nakka'): {'k': 1.133, 'M': 0.04202, 'To_T': 1720, 
                        'rho_T': 1.889, 'a': 8.26, 'n': 0.319},
    ('ASSS', 'Altum'): {'k': 1.131402, 'M': 0.0336351, 'To_T': 2439.8, 
                        'rho_T': 1.7613, 'a': 8.26, 'n': 0.319},
    ('Altum2', 'Altum'): {'k': 1.225644, 'M': 0.03288467, 'To_T': 2108.41, 
                        'rho_T': 1.9527, 'a': 2.57, 'n': 0.371},
    ('KNDX', 'Nakka'): {'k': 1.131, 'M': 0.04242, 'To_T': 1710,
                         'rho_T': 1.879}, # Tienes distintos valores 
    ('KNSB', 'Nakka'): {'k': 1.137, 'M': 0.0399, 'To_T': 1600,
                         'rho_T': 1.841, 'a': 5.130, 'n': 0.22},
    ('KNER', 'Nakka'): {'k': 1.140, 'M': 0.03858, 'To_T': 1608,
                         'rho_T': 1.820, 'a': 2.9, 'n': 0.4},
    ('KNMN', 'Nakka'): {'k': 1.136, 'M': 0.03983, 'To_T': 1616,
                         'rho_T': 1.854, 'a': 5.130, 'n': 0.22},
    ('KNFR', 'Nakka'): {'k': 1.131, 'M': 0.04242, 'To_T': 1710,
                         'rho_T': 1.942, 'a': 7.4, 'n': 0.250},
    ('KNPSB', 'Nakka'): {'k': 1.163, 'M': 0.03639, 'To_T': 1858,
                          'rho_T': 1.923, 'a': 6.5, 'n': 0.628},
    ('RNX-57', 'Nakka'): {'k': 1.128, 'M': 0.04520, 'To_T': 1644,
                          'rho_T': 1.8445, 'a': 1.95, 'n': 0.477},
    ('RNX-71V', 'Nakka'): {'k': 1.14, 'M': 0.04183, 'To_T': 1492,
                          'rho_T': 1.8161, 'a': 2.57, 'n': 0.371},
    ('SNX', 'Altum'): {'k': 1.1405, 'M': 0.034957, 'To_T': 2230.3,
                          'rho_T': 1.932, 'a': 2.57, 'n': 0.371},
    ('ACPC', 'MIT'): {'k': 1.21, 'M': 0.02367, 'To_T': 2821,
                       'rho_T': 1.680,'a': 3.237295, 'n': 0.327392,},
    ('ACPC', 'BPS'): {'k': 1.199, 'M': 0.02378, 'To_T': 2884,
                       'rho_T': 1.682,'a': 2.75225, 'n': 0.456,} 

}

# %%

# Nombre del proyecto
Nombre = 'Tlamatiliztli'
formato ='csv' # Posibles: 'csv', 'txt', 'excel'
Manofacturador = 'Leonel'
# %%
Propelente = 'Altum2' # Posibles: KNSU, KNDX, KNSB, KNER, KNMN, KNFR, KNPSB, RNX-57, RNX-71V, ACPC
Creador = 'Altum' # Posibles: Nakka, UDEG SPACE, MIT, Altum
Altura = 1350 # m, Altura sobre el nivel del mar
# %%

# Obtener los valores correspondientes
Prop = Datos_propelente.get((Propelente, Creador))

if Prop:
    k = Prop['k']         # relación de calores especificos
    M = Prop['M']         # kg/mol masa molecular de los gases
    To_T = Prop['To_T']   # K temperatura de combustion teorica
    rho_T = Prop['rho_T'] # g/cm^3 densidad teorica
    a = Prop['a']   # mm/s/Mpa^n
    n = Prop['n']

    print(f"Valores para Propelente: {Propelente}, Creador: {Creador}")
    print(f"k = {k}, M = {M}, To_T = {To_T}, rho_T = {rho_T}, a = {a}, n = {n}" )





alpha = 12*np.pi/180 # angulo de la divergencia de la tobera 
beta  = 35*np.pi/180 # Angulo de convergencia

# Diametros
D_ext = 57.5 # mm
Dc    = 53.5      # mm Diametro de la camara 
Df    = 4       # in, diametro fuselaje
Dt    = 12      # mm, de pulgadas a mm [(in)*(mm/in)]
Dp    = 26.6 #Dt*np.sqrt(3) #, de pulgadas a mm [(in)*(mm/in)] # Diametro del port
In    = 0       # ancho de los inhibidores y/o revestimiento del grano


# Areas
Af = np.pi*(Df/(2*39.3701))**2 # [m**2] # área  del fuselaje
Ac = np.pi*(Dc/2)**2 # mm**2 # área de la camara
At = np.pi*(Dt/2)**2 # mm**2 # área de la garganta
Ap = np.pi*(Dp/2)**2 # mm**2 # área del port
AcAt = Ac/At
ApAt = Ap/At

L_t    = 4 # mm
L_Dt   = L_t/Dt
Cd     =  0.46  # coeficiente de arrastre

if L_Dt > 0.45:
    N_throat = 0.95
else:
    N_throat = 0.99 - 0.0333*L_Dt
N_noz = 0.85 # # eficiencia de la tobera
N_com = 0.95 # eficiencia de combustion
N_div = 1/2*(1 + np.cos(alpha)) # factor de corrección por la divergencia
N_po = 0.95 # factor de corrección de la presion de la camara
N_skin = 0.99


rho_rat = 0.92
rho = rho_T*rho_rat    # g/cm^3 densidad medida
To = To_T*N_com # Temperatura "real"
RelacionSeparacion = 0.4 # 40 % segun el criterio de Summerfield

# %%


Pa = Propiedades_Atm(Altura)[0]/6894.76 # psi, presion atmosferica
Pe = Pa #P0*((k+1)/2)**(-k/(k-1)) #psi, presion de salida de la tobera
P0 = 400 #Pe*((k+1)/2)**(k/(k-1)) #psi, presion de la camara objetivo

# Mach
Mt = 1
Me = mach(Pe, P0, k)


#%%
# Condiciones del Cohete
theta_0 = np.pi/2 #angulo
MasaCarcasa = .920 # kg
L_motor = 365 # mm
m = 5 # masa del cohete sin propelente

mdot_obj = mdot(At*10**-6, (P0 + Pa)*6894.76, k,To, cte.R/M)

# Medidas de salida
AeAt = RatArea(Me, Mt, k)
Ae = At*AeAt
De = diametro(Ae)


Tb_aprox = (Dc-Dp)/(2*a*(P0/145.038)**n)
print(f'tb aprox = {Tb_aprox}')
# %%



Datos = 100000 # cantidad de datos para la regresión del grano
h0    = 1000 # altura deseada 


Cf = N_throat*N_div*N_noz*(N_skin*Cf2_(Pe, (P0 + Pa), Pa, At, Ae, k)+(1-N_skin))
c_ast = ex_vel(cte.R/M, To, k)
Ve1 = c_ast*Cf # Velocidad de salida de los gases de la tobera

k1 = 1/(2*m)*Propiedades_Atm(Altura*1.05)[1]*Cd*Af 
DeltaV = np.sqrt((cte.g/k1)*(np.exp(2*k1*(h0*0.98))-1))
print(f'Delta V {DeltaV}')
mp0 = m*(np.exp((DeltaV + cte.g*Tb_aprox)/Ve1)-1)

#mp0 = m*(np.exp((np.sqrt(2*cte.g*h0)+cte.g*Tb_aprox)/Ve1)-1) # masa minima de propelente

print(f'mp0 = {mp0}')
m_p = [mp0] # se utiliza la formula del método que no se considera la fuerza de arrastre para la primera iteración
M_tot = [m + mp0]

# Si es 1 significa que se toma en cuenta si es 0 no se toma en cuenta, 
# si se pone algun otro número se tendran resultados erroneos
Bs = [0, 1, 1] # Superficies: exterior, nucleo, caras
L_0 = (3*Dc+Dp)/2
N = 1 #round(4*mp0/L_0/rho/(Dc**2-Dp**2)/np.pi) # Número de segmentos del propelente

# Ciclo iterativo
Tolerancia = 10 # Error de la altura alcanzada en metros

A_t = At*10**-6 # Área de la garganta en m^3
A_e = Ae*10**-6 # Área de la salida en m^3
P_a = Pa/145.038 # Mpa Presión atmosferica








# Condiciones de calculos
PresionCambiante = True
Vuelo = True

if  Vuelo:  
    i = 0
    Iteraciones = []
    h = [0]
    while np.abs(h0-h[i]) > Tolerancia:
        print(f'Iteración = {i + 1}')
        Iteraciones.append(i)


        print(f'N granos {round(4*m_p[i]/(L_0)/(rho*1e-6)/(Dc**2-Dp**2)/np.pi) }')
        fase2_activa = False
        tb = 0
        Guardado = 0

        # Área de quemado
        V0 = m_p[i]/(rho*1000) # m^3 Volumen con densidad experiemental del grano
        Lg0 = V0*10**9/(np.pi/4*((Dc - 2*In)**2-Dp**2)) # mm  Longitud del grano con densidad experimental
        
        Vc =Ac*(Lg0*1.01)*10**-9 #Vc = (V0)/(1-ApAt*At*4/(np.pi*(Dc - 2*In)**2)) # m^3 Volumen disponible de la camara
        Lc = Vc*10**9/Ac # mm Longitud de la cámara
        v0 = Vc - V0 # m^3 volumen libre de la cámara
        

        t_span = (0, 50)
        t_eval = np.linspace(0, 50, Datos)

    #---------------------------------------------------------
        if PresionCambiante == True:
            #print(f'm_tot {M_tot[i]}')
            p_p1 = [a * 10**(-3-6*n), n, Dc* 10**-3, Dp * 10**-3, Lg0*10**-3, Bs, In, rho*1000, k, M, To, A_t, A_e, Me, N]
            p_P = [M_tot[i], Cd , Af, theta_0, p_p1]
            CI1 = [0, 0, 0, 0, 0, 0, v0, 0, 0, 0]

            Solucion_1 = solve_ivp(Ec_Presion2, t_span, CI1, args=(p_P,), events = event_pressure_zero2,
                                    t_eval=t_eval, rtol=1e-9, atol=1e-9)

            t_sol1 = Solucion_1.t
            x_sol1, Vx_sol1, y_sol1, Vy_sol1, Presion_camara,_,_,_,_,_  = Solucion_1.y ### checar (_,_)
            Presion_camara = Solucion_1.y[4]/1e6 # Mpa
            t_thrust = Solucion_1.t

            
            PresionAtmosferica = Propiedades_Atm(y_sol1 + Altura)[0]*1e-6
            PresionAbsoluta = Presion_camara + PresionAtmosferica

            r = a*(PresionAbsoluta)**n # mm/s
            X_num = Solucion_1.y[9]*1e3
            X = r*t_thrust            # mm

            Ab = Ab_(Dc, In, Dp, Lg0, X_num, N, Bs)
            Kn = Ab/At

            P_e = P_Po(PresionAbsoluta , Me, k)

            C_f = N_throat*N_div*N_noz*(N_skin*Cf_(P_e , 
                                                PresionAbsoluta, PresionAtmosferica , A_t, A_e, k)+(1-N_skin))
            F = C_f*At*(PresionAbsoluta)

            # Verificar que las longitudes coincidan
            if len(t_thrust) != len(F):
                print(len(t_thrust), len(F))
                raise ValueError("La longitud de 't_thrust' es diferente de la longitud de 'F'")

            # Integrar la función spline en el intervalo deseado, por ejemplo, de t_min a t_max

            M_e = mach(P_e, PresionAbsoluta, k)
            Ae_At = RatArea(M_e, 1, k)

            mdot_ = mdot(A_t, PresionAbsoluta*1e6, k, To, cte.R/M)

            F_spline = CubicSpline(t_thrust, F)
            m_noz_spline = CubicSpline(t_thrust, mdot_)

            def F_(t, t_thrust):
                return np.piecewise(t, [t <= t_thrust, t > t_thrust], [F_spline(t), 0])

            def dm_(t, t_b):
                return np.piecewise(t, [t <= t_b, t > t_b], [m_noz_spline(t), 0])

            def F_spline_func(t):
                return F_spline(t)


            t_nuevo = np.linspace(t_thrust[0], t_thrust[-1], 1000)
            
            m_noz_nuevo = m_noz_spline(t_nuevo)
            F_nuevo = F_spline(t_nuevo)

            
        else :
            #print(f'm_tot {M_tot[i]}')
            p_p1 = [a * 10**(-3-6*n), n, Dc* 10**-3, Dp * 10**-3, Lg0*10**-3, Bs, In, rho*1000, k, M, To, A_t, A_e, Me, N]
            p_P = [M_tot[i], Cd , Af, theta_0, p_p1]
            CI1 = [0, 0, 0, 0, 0, 0, v0, 0, 0, 0]

            Solucion_1 = solve_ivp(Ec_Presion2, t_span, CI1, args=(p_P,), events = event_pressure_zero2,
                                    t_eval=t_eval, rtol=1e-10, atol=1e-10)

            t_sol1 = Solucion_1.t
            x_sol1, Vx_sol1, y_sol1, Vy_sol1, Presion_camara,_,_,_,_,_  = Solucion_1.y ### checar (_,_)
            Presion_camara = Solucion_1.y[4]/1e6 # Mpa
            t_thrust = Solucion_1.t


            
            PresionAtmosferica = np.ones(len(t_thrust))*P_a
            PresionAbsoluta = Presion_camara + PresionAtmosferica

            r = a*(PresionAbsoluta)**n # mm/s
            X_num = Solucion_1.y[9]*1e3
            X = r*t_thrust            # mm

            Ab = Ab_(Dc, In, Dp, Lg0, X_num, N, Bs)
            Kn = Ab/At

            P_e = P_Po(PresionAbsoluta , Me, k)

            C_f = N_throat*N_div*N_noz*(N_skin*Cf_(P_e, 
                                                PresionAbsoluta, PresionAtmosferica , A_t, A_e, k)+(1-N_skin))
            F = C_f*At*(PresionAbsoluta)

            # Verificar que las longitudes coincidan
            if len(t_thrust) != len(F):
                print(len(t_thrust), len(F))
                raise ValueError("La longitud de 't_thrust' es diferente de la longitud de 'F'")

            # Integrar la función spline en el intervalo deseado, por ejemplo, de t_min a t_max

            M_e = mach(P_e, PresionAbsoluta, k)
            Ae_At = RatArea(M_e, 1, k)

            mdot_ = mdot(A_t, PresionAbsoluta*1e6, k, To, cte.R/M)

            F_spline = CubicSpline(t_thrust, F)
            m_noz_spline = CubicSpline(t_thrust, mdot_)

            def F_(t, t_thrust):
                return np.piecewise(t, [t <= t_thrust, t > t_thrust], [F_spline(t), 0])

            def dm_(t, t_b):
                return np.piecewise(t, [t <= t_b, t > t_b], [m_noz_spline(t), 0])

            def F_spline_func(t):
                return F_spline(t)


            t_nuevo = np.linspace(t_thrust[0], t_thrust[-1], 1000)
            
            m_noz_nuevo = m_noz_spline(t_nuevo)
            F_nuevo = F_spline(t_nuevo)

        Po_max = np.nanmax(Presion_camara) 
        t_thr = t_thrust[-1] # tiempo de empuje
        t_min = t_thrust[0]
        t_max = t_thr
        integral_result, error = quad(F_spline_func, t_min, t_max, limit=100)
        I_t = integral_result

        mp_salida, _ = quad(m_noz_spline, t_min, t_max, limit = 100)

        print(f'mp calculado {m_p[i]}, salida {mp_salida}, relación {mp_salida/m_p[i]*100}')


    
        x_fin1, Vx_fin1, y_fin1, Vy_fin1 = [x_sol1[-1], Vx_sol1[-1], y_sol1[-1], Vy_sol1[-1]]

        # Calculos fase 2
        

        CI2 = [x_fin1, Vx_fin1, y_fin1, Vy_fin1] # Condiciones iniciales
        t_2 = (t_thr, 200) # Si en la gráfica no se ve todo el descenso cambiar el tiempo final
        p2 = [m, Cd , Af, t_thr, tb, theta_0, 0]
        Solucion_2 = solve_ivp(Sistema, t_2, CI2, args=(p2,), t_eval=np.linspace(t_thr, 200, 2000), events=evento_y_cero)


        t_sol2 = Solucion_2.t
        x_sol2, Vx_sol2, y_sol2, Vy_sol2 = Solucion_2.y

        # Si necesitas las soluciones en el mismo formato que odeint
        
        # Metodo de convergencia
        h.append(max(y_sol2))
        xd = True

        if np.abs(h0-h[i+1]) > Tolerancia:
            if i==200:
                    break
            
            if i < 50:

                zeta = m_p[i]/M_tot[i]
                # Función que calcula el integrando para cada valor de t
                def integrand(t, u, rho):
                    return 0.5 * rho * u**2 / (1 - zeta * t / tb)

                # Valores del integrando
                f_values = integrand(t_sol1, Vy_sol1, Propiedades_Atm(y_sol1 + Altura)[1])

                # Calcular la integral usando el método de Simpson
                B = simps(f_values, t_sol1)

                print(f'B { B*Cd*Af/M_tot[i]}, {cte.g * tb }')

                #print('Primera')
                h_fin1 = y_sol1[-1]
                rhoAtm_prom = np.mean(Propiedades_Atm(y_sol1+Altura)[1])
                DeltaT = t_thrust[1]-t_thrust[0]

                CD_n = 1/(2 * (m + m_p[i]- mp_salida)) * rhoAtm_prom * Cd *Af
                C_Star = At* np.sum(Presion_camara)*DeltaT / m_p[i]
                Ve2 = np.mean(C_Star*Cf)
                
                DeltaV2 = np.sqrt((cte.g / CD_n )*(np.exp(2 * CD_n  *(h0 - h_fin1))-1))
                mp_inc = (m)*(np.exp((DeltaV2 + cte.g * tb + B*Cd*Af/M_tot[i]) / Ve2)-1)
                print(f'mp0 nuevo {mp_inc}', f'densida prom {rhoAtm_prom}', f'C* ideal{c_ast} medido {C_Star} ')
                print(f'Ve medido {Ve2}, ideal {Ve1}')

                Inc_M = mp_inc - M_tot[i] + m

            if i > 150:
                if i == 151 and xd == True and np.abs((h0-h[i+1])/h[i+1]) < 1 :
                    print('Ajuste de la masa')
                #slope, intercept, r_value, p_value, std_err = linregress(h[:len(m_p)], m_p)
                #Inc_M = slope*h0 + intercept - M_tot[i] + m
                
                coeficientes = np.polyfit(h[:len(m_p)], m_p, 2)
                polinomio = np.poly1d(coeficientes)

                Inc_M = polinomio(h0) - M_tot[i] + m

                if i==200:
                    break

            else:

                #print('Cambio en la convergencia')
                Inc_M = m_p[i]*0.1*np.sign((h0-h[i+1])/h[i+1])#/h[i+1]*(h0-h[i+1])

                #print(f'Inc_M 2={Inc_M}')

            #if np.abs((h0-h[i+1])/h[i+1]) < 1:
            #    Inc_M = m_p[0]*0.1*np.sign((h0-h[i+1])/h[i+1])
            #    #print(f'Inc_M 1={Inc_M}')
            #else:
            #    Inc_M = m_p[0]*0.1*np.sign((h0-h[i+1])/h[i+1])#/h[i+1]*(h0-h[i+1])
                #print(f'Inc_M 2={Inc_M}')

            M_tot.append(M_tot[i] + Inc_M)
            m_p.append(M_tot[i+1]-m)
        ####
        print(f'Altura alcanzada = {h[i+1]} m')
        i+=1
else:
    print('Vuelo desactivado')
    
    fase2_activa = False
    tb = 0
    Guardado = 0
     
    L_0 = 380 # mm 
    Lg0 = N*L_0
    Lc = Lg0*1.01
    Vc =  Ac*Lc*1e-9
    V0 = Lg0*np.pi/4*((Dc - 2*In)**2-Dp**2)*1e-9
    m_p.append(V0*rho*1000)
    v0 = Vc - V0

    CI_P = [0, 0, v0, 0, 0, 0]

    p_P = [a * 10**(-3-6*n), n, Dc* 10**-3, Dp * 10**-3, Lg0*10**-3, Bs, In, rho*1000, k, M, To, A_t, N]

    #print(f'parametros {p_P}')

    t_span = (0, 50)
    t_eval = np.linspace(0, 50, Datos)
    sol = solve_ivp(Ec_Presion, t_span, CI_P, args=(p_P,), events=event_pressure_zero,
                     t_eval=t_eval, rtol=1e-12, atol=1e-14)

    Presion_camara = sol.y[0]/1e6 # Mpa
    t_thrust = sol.t
    PresionAtmosferica = np.ones(len(t_thrust))*P_a
    PresionAbsoluta = Presion_camara + PresionAtmosferica
    #print(f'Po={Presion_camara}')
    #print(f't={sol.t}')

    r = a*(PresionAbsoluta)**n # mm/s
    X_num = sol.y[5]*1e3
    X = r*t_thrust            # mm

    print(f'Regresion maxima :[{np.nanmax(X_num)}]', f'Espesor de la red : [{(Dc-Dp)/2}]')

    Ab = Ab_(Dc, In, Dp, Lg0, X_num, N, Bs)
    Kn = Ab/At
    Po_max = np.nanmax(Presion_camara) 



    P_e = P_Po(PresionAbsoluta, Me, k)
    C_f = N_throat*N_div*N_noz*(N_skin*Cf_(P_e, PresionAbsoluta, PresionAtmosferica, A_t, A_e, k)+(1-N_skin))
    F = C_f*At*(PresionAbsoluta)

    M_e = mach(P_e, PresionAbsoluta, k)
    Ae_At = RatArea(M_e, 1, k)


    mdot_ = mdot(A_t, PresionAbsoluta*1e6, k, To, cte.R/M)


    if len(t_thrust) < 2 or len(F) < 2:
        raise ValueError("`t_thrust` y `F` deben contener al menos 2 elementos cada uno.")

    # Crear una función spline para F en función del tiempo t_thrust
    F_spline = CubicSpline(t_thrust, F)
    m_noz_spline = CubicSpline(t_thrust, mdot_)

    def F_spline_func(t):
        return F_spline(t)

    # Integrar la función spline en el intervalo deseado, por ejemplo, de t_min a t_max
    t_thr = t_thrust[-1] # tiempo de empuje
    t_min = t_thrust[0]
    t_max = t_thr
    integral_result, error = quad(F_spline_func, t_min, t_max, limit=100)
    I_t = integral_result
    

    #Graficar para verificar la interpolación

    t_nuevo = np.linspace(t_thrust[0], t_thrust[-1], 1000)
    m_noz_nuevo = m_noz_spline(t_nuevo)


    F_nuevo = F_spline(t_nuevo)


# Figura 1
plt.figure(figsize=(7,7))
plt.plot(t_thrust, (Presion_camara)*145.038 , 'b-', markersize='2.5')
#plt.ylim(bottom=0)
plt.ylabel('Presión, Mpa')
plt.xlabel('Tiempo, s')
plt.title('Presión de la cámara')
plt.grid()

# Figura 2
plt.figure(figsize=(7,7))
plt.plot(t_thrust, F, 'b-')
#plt.ylim(bottom=0)#, top= 9000)
plt.ylabel('Empuje, N')
plt.xlabel('Tiempo, s')
plt.title('Gráfica de empuje')
plt.grid()


plt.figure(figsize=(7,7))
plt.plot(t_thrust, C_f, 'b.')
plt.ylabel('C_f')
plt.xlabel('Tiempo, s')
plt.title('Gráfica de Cf')
#plt.legend()plt.grid()



# Figura 3
fig1, ax = plt.subplots(2, 1, figsize=(7, 7))

ax[0].plot(t_thrust, F, 'bo', markersize='2.5', label='Datos originales')
ax[0].plot(t_nuevo, F_nuevo, 'r-', label='Interpolación Spline')
ax[0].legend()
ax[0].set_xlabel('Tiempo, s')
ax[0].set_ylabel('Empuje, N')
ax[0].set_title('Interpolación de F usando Spline Cúbico')
ax[0].set_ylim(bottom = 0)
ax[0].grid()

ax[1].plot(t_thrust, mdot_, 'bo', markersize='2.5', label='Datos originales' )
ax[1].plot(t_nuevo, m_noz_nuevo, 'r-')
ax[1].set_xlabel('Tiempo, s')
ax[1].set_ylabel('Flujo masico, kg/s')
ax[1].set_title(r'Interpolación de $\dot{m}$ usando Spline Cúbico')
ax[1].set_ylim(bottom=0)
ax[1].grid()
fig1.tight_layout()




## Figura 5: Velocidad
if Vuelo:

    # Figura 4
    plt.figure(figsize=(7,7))

    Ymax = np.where(y_sol2 == max(y_sol2)) # indice donde la altura es maxima

    plt.plot(t_sol1, y_sol1, 'r-', label = 'Fase propulsada' ) # posición en y 
    plt.plot(t_sol2, y_sol2, 'b-', label = 'Fase sin propulsión') # posición en y

    plt.xlabel('Tiempo, s')
    plt.ylabel(r'Altura $y$, m')


    plt.xlim(left=0, right = t_sol2[Ymax[0]])
    plt.ylim(bottom=0, top = max(y_sol2) + 5 )
    plt.legend()
    plt.grid()

    fig3, (ax1, ax2) = plt.subplots(2, figsize=(7, 7))

    ax1.plot(t_sol1, Vy_sol1,'r-', markersize='2', label = 'Fase propulsada') # Velocidad en y 
    ax2.plot(t_sol1, Vx_sol1,'r-', markersize='2', label = 'Fase propulsada') # Velocidad en x
    ax1.plot(t_sol2, Vy_sol2, 'b-', label = 'Fase sin propulsión') # Velocidad en y
    ax2.plot(t_sol2, Vx_sol2, 'b-', label = 'Fase sin propulsión') # Velocidad en x

    ax1.set_xlabel('Tiempo, s')
    ax1.set_ylabel(r'$V_{y}$, m/s')
    #ax1.set_ylim(bottom=0, top = max(Sol_2[:, 3]) + 5 )
    ax1.legend()
    ax1.grid()

    ax2.set_xlabel('Tiempo, s')
    ax2.set_ylabel(r'$V_{x}$, m/s')
    ax2.legend()
    ax2.grid()

    coeficientes = np.polyfit(m_p, h[:len(m_p)], 1)
    polinomio = np.poly1d(coeficientes)
    x_mp = np.linspace(min(m_p), max(m_p), 1000)
    y_mp = polinomio(x_mp)


    # Crear figura con GridSpec
    fig8 = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1])  # 3:1 ratio for height

    # Primer subplot (fila 1, columnas 0-2, abarca ambas columnas)
    ax1 = fig8.add_subplot(gs[0, :])
    ax1.plot(m_p, h[1:], 'r.', markersize='2.5' )
    #ax1.plot(x_mp, y_mp, 'm-' )
    ax1.hlines(h0, min(m_p), max(m_p), color="black", linestyle="dashed", label='Apogeo deseado')
    ax1.fill_between(m_p, h0 - Tolerancia, h0 + Tolerancia, color='#6495ED', alpha=0.3, label='Tolerancia')
    ax1.set_xlabel(r'$m_{p}$, kg')
    ax1.set_ylabel('Apogeo, m')
    ax1.grid()
    ax1.legend()
    #ax1.set_title('Figura Grande')

    # Segundo subplot (fila 2, columna 0)
    ax2 = fig8.add_subplot(gs[1, 0])
    ax2.plot(Iteraciones, h[1:], 'g-', markersize='2.5')
    ax2.hlines(h0, min(Iteraciones), max(Iteraciones), color="black", linestyle="dashed", label='Apogeo objetivo')
    ax2.set_xlabel('Iteración')
    ax2.set_ylabel('Apogeo, m')
    ax2.grid()
    ax2.legend()
    ax2.set_title('Figura Pequeña 1')

    # Tercer subplot (fila 2, columna 1)
    ax3 = fig8.add_subplot(gs[1, 1])
    ax3.plot(Iteraciones, m_p, 'm-', markersize='2.5', )
    ax3.hlines(m_p[-1], min(Iteraciones), max(Iteraciones), color="black", linestyle="dashed", label=r'$m_{p}$ a $h_{0}$')
    ax3.set_xlabel('Iteración')
    ax3.set_ylabel(r'$m_{p}$, kg')
    ax3.grid()
    ax3.legend()
    ax3.set_title('Figura Pequeña 2')


else: next

# Figura 6 Relaciones de presión atmosferica y presión de salida



#P_ro = (P_Po(Presion_camara*145.038 + Pa, Me, k)) # kp
P_re = (np.array(P_e + PresionAtmosferica))*145.038 # Mpa
#Pe_opt = (PresionAtmosferica - Presion_camara*coeficiente_empuje_optimo(k, Me)/AeAt)*145.038


plt.figure(figsize=(7,7))
plt.plot(t_thrust, Presion_camara*RelacionSeparacion*145.038/Pa, color='#DC143C', label=f'{RelacionSeparacion}'r'$P_0/P_{a0}$')
plt.plot(t_thrust, PresionAtmosferica*145.038 ,'k--', label = r'$P_{a}$')
#plt.plot(t_thrust, P_re, '#FFA500', linewidth='2', label = r'$P_{e}$')
plt.plot(t_thrust, (PresionAtmosferica- Presion_camara*Cf2_(P_e,PresionAbsoluta, P_e, At, Ae, k)/AeAt)*145.038, '#FFA500', linewidth='2', label = r'$P_{a}-P_{o}Cf_{opt}/e$')
plt.plot(t_thrust, P_e*145.038, linewidth='2', label = r'$P_{e}$')
#plt.title(r'Relación de presiones con $P_{e}$')
plt.xlabel('Tiempo, s')
plt.ylabel('Presión, psi')
plt.legend()
plt.grid()

plt.figure(figsize=(7,7))
plt.plot(X, Kn, 'b-')
plt.xlabel('Web regression, mm')
plt.ylim(bottom = 0, top =max(Kn))
#plt.legend()
plt.grid()

plt.figure(figsize=(7,7))
plt.plot(t_thrust, X_num, '.r', label= 'dxdt=r', markersize= 2)
plt.plot(t_thrust, X, '.b', label='x=rt', markersize= 2)
plt.vlines(tb, 0, np.max(X_num)*1.1,color = "black", linestyle = "dashed" )
plt.hlines((Dc-Dp)/2, 0, np.max(t_thrust),color = "green", linestyle = "dashed" )
plt.xlabel('tiempo, s')
plt.ylabel('Web regression, mm')
plt.legend()
plt.ylim(bottom = 0, top =np.max(X_num)*1.1)
plt.legend()
plt.grid()






t_Brn = tb
Po_max = (Po_max + P_a)*1000000/6895
Po_prom = np.mean(Presion_camara)*1000000/6895
F_max = max(F)
F_prom = np.mean(F)
Cf_max = np.nanmax(C_f)
Cf_prom = np.nanmean(C_f)
Cf_min = np.nanmin(C_f)
It = I_t
Isp = It/(cte.g*m_p[-1])
mp = m_p[-1]
m_tot = M_tot[-1]
mdot_prom =np.nanmean(m_noz_nuevo)

tic = (Dc-Dp)/2
Ab_final = N*np.pi*Dc*(L_0-2*tic)

Kn_max = np.nanmax(Kn)
Kn_prom = np.nanmean(Kn)
Kn_min = np.nanmin(Kn)
Kn_in = Kn[0]
Kn_final = Ab_final/At

Po_prom_aprox = Po_aprox(Kn_prom,k, cte.R/M, To,rho,a, n)
Po_max_aprox = Po_aprox(Kn_max,k, cte.R/M, To,rho,a, n)
Po_min_aprox = Po_aprox(Kn_min,k, cte.R/M, To,rho,a, n)
Po_fin_aprox = Po_aprox(Kn_final,k, cte.R/M, To,rho,a, n)


Pe_Pa = np.nanmean(P_re/Pa)
Pa_min = np.nanmin(PresionAtmosferica*145.038)

Pe_max = np.nanmax(P_re)

N_efi= N_throat*N_div*N_noz*N_skin

# Dimensiones de la tobera con presión promedio
AeAt_p = np.nanmean(Ae_At)
De_p = 2*np.sqrt(AeAt_p*At/np.pi)

L_noz = (De_p-Dt)/(2*np.tan(alpha)) + (Dc-Dt)/(2*np.tan(beta)) + L_t
DeltaT = t_thrust[1]- t_thrust[0]

# Datos
print('')
print('Datos relevantes')
print('')
print(f'Tiempo de quemado = {t_Brn:.4f}', f'Tiempo de empuje = {t_thr:.4f}')
print(f'Tiempo de quemado aprox = {Tb_aprox:.4f}', f'Delta t = {DeltaT}')
print(f'Po max = {Po_max:.3f}', f'Po prom = {Po_prom:.3f}')
print(f'Po max aprox = {Po_max_aprox:.3f}', f'Po prom aprox = {Po_prom_aprox:.3f}', f'Po min aprox = {Po_min_aprox:.3f}', f'Po final aprox = {Po_fin_aprox:.3f}')
print(f'Pe_max {Pe_max:.3f}', f'Rat Pe_Pa {Pe_Pa:.4f}')
print(f'F max = {F_max:.2f}', f'F prom = {F_prom:.2f}' )
print(f'Cf max = {Cf_max:.4f}', f'Cf prom = {Cf_prom:.4f}', f'Cf min = {Cf_min:.4f}' )
print(f'It = {It:.2f}', f'Isp = {Isp:.2f}')
print(f'm prop = {mp:.4f}')
print(f'Me = {Me:.3f}', f'mdot ={mdot_prom:.3f}', f'mdot obj ={mdot_obj:.3f}', f'R/M = {cte.R/M:.4f}')
print(f'Kn max = {Kn_max:.2f}', f'Kn in = {Kn_in:.2f}', f'Kn min = {Kn_min:.2f}', f'Kn prom = {Kn_prom:.2f}', f'No. Granos {N}')
print(f'Lg0 = {Lg0:.3f}', f'L0 = {L_0:.3f}', f'L0 neutro {(3*Dc+Dp)/2}')
print(f' Vg = {V0}, Vc = {Vc}')


print('')
print('Dimensiones de la cámara y tobera')
print('')
print(f'L_c = {Lc:.2f}',f'L_noz = {L_noz:.3f}',f' Dc = {Dc:.4f}', f' Dt = {Dt:.4f}', f' De = {De:.4f}', f'De a Po prom = {De_p:.4f}', f'Dp = {Dp:.4f}')
print(f'Ap/At = {ApAt:.3f}',f'Ac/At = {AcAt:.3f}')
print(f'Ae/At a Po(obj) = {AeAt:.3f}', f' Ae/At max = {np.nanmax(Ae_At):.3f}', f'Ae/At prom = {AeAt_p:.3f}')
print(f'Eficiencia de tobera {N_efi:.3f}')


if Vuelo:
    h_PaMin = y_sol1[-1]
    print('')
    print('Vuelo y cohete')
    print('')
    print( f'Pa min(h ={h_PaMin:.2f}): {Pa_min:.3f}')
    print(f'm = {m}', f'm tot = {m_tot:.4f}')

    data = {
    'Altura objetivo' : [h0],
    'Tiempo de quemado (s)' : [t_Brn],
    'Tiempo de empuje (s)' : [t_thr],
    'Presión Objetivo (psi)' : [P0],
    'Presión máxima (psi)' : [Po_max],
    'Presión promedio (psi)' : [Po_prom],
    'Empuje máximo (N)' : [F_max],
    'Empuje promedio (N)' : [F_prom],
    'Coeficiente de empuje maximo' : [Cf_max],
    'Coeficiente de empuje promedio' : [Cf_prom],
    'Impulso total (Ns)' : [It],
    'Impulso especifico (s)' : [Isp],
    'Masa del cohete (Kg)' : [m],
    'Masa propelente (Kg)' : [mp],
    'Masa total (Kg)' : [m_tot],
    'Longitud de la cámara (mm)' : [Lc],
    'Longitud de la garganta (mm)': [L_t],
    'Longitud del grano (mm)' : [Lg0],
    'Diametro cámara (mm)' : [Dc],
    'Diametro garganta (mm)' : [Dt], 
    'Diametro salida (mm)' : [De_p],
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
    'Velocidad caracteristica (m/s)' : [c_ast],
    'Coeficiente de presión (Mpa^1/n)' : [a],
    'Exponente de presión ' : [n],
    'Angulo de lanzamiento ' : [theta_0*180/np.pi],
    'Diametro fuselaje (mm)' : [Df],
    'Coeficiente de arrastre cohete' : [Cd],
}

else: 
    data = {
        'Tiempo de quemado (s)' : [t_Brn],
        'Tiempo de empuje (s)' : [t_thr],
        'Presión Objetivo (psi)' : [P0],
        'Presión máxima (psi)' : [Po_max],
        'Presión promedio (psi)' : [Po_prom],
        'Empuje máximo (N)' : [F_max],
        'Empuje promedio (N)' : [F_prom],
        'Coeficiente de empuje maximo' : [Cf_max],
        'Coeficiente de empuje promedio' : [Cf_prom],
        'Impulso total (Ns)' : [It],
        'Impulso especifico (s)' : [Isp],
        'Masa del cohete (Kg)' : [m],
        'Masa propelente (Kg)' : [mp],
        'Masa total (Kg)' : [m_tot],
        'Longitud de la cámara (mm)' : [Lc],
        'Longitud de la garganta (mm)': [L_t],
        'Longitud del grano  (mm)' : [Lg0],
        'Diametro cámara (mm)' : [Dc],
        'Diametro garganta (mm)' : [Dt], 
        'Diametro salida (mm)' : [De_p],
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
        'Velocidad caracteristica (m/s)' : [c_ast],
        'Coeficiente de presión (Mpa^1/n)' : [a],
        'Exponente de presión ' : [n],
    }


data_tobera = {
    'Alpha': [alpha],
    'Beta': [beta],
    'Longitud tobera': [L_noz],
    'Longitud de la garganta (mm)': [L_t],
    'Presión promedio': [Po_prom],
    'Relación de expansión': [AeAt_p],
    'Ac/At': [AcAt],
    'De' : [De_p],
    'Dt' : [Dt],
    'Dc' : [Dc],

}

data3 = {
    '  Tiempo  ': t_thrust,
    '  Presión  ': Presion_camara,
    '  Empuje  ': F,
    '  Área de quemado  ': Ab,
    '  Kn  ': Kn,
    '  x  ': X_num,

}

# Crear una carpeta
carpeta = Nombre
if not os.path.exists(carpeta):
    os.makedirs(carpeta)

# Crear DataFrames
df1 = pd.DataFrame(data)
df2 = pd.DataFrame(data_tobera)
df3 = pd.DataFrame(data3)

# Guardar los DataFrames en la carpeta creada
if formato == 'csv':
    df1.to_csv(f'{carpeta}/{Nombre}.csv', index=False)
elif formato == 'txt':
    df1.to_csv(f'{carpeta}/{Nombre}.txt', sep='\t', index=False)
elif formato == 'excel':
    df1.to_excel(f'{carpeta}/{Nombre}.xlsx', index=False)
else:
    print("Formato no compatible")

df2.to_csv(f'{carpeta}/Tobera_{Nombre}.csv', index=False)
df3.to_excel(f'{carpeta}/Datos_{Nombre}.xlsx', index=False)

if L_motor == False:
    L_motor = round(Lc*1.2)
else :
    next

#generar_archivo_eng('Motor.eng', 'L1341', 75, round(Lc), m_p[-1], MasaCarcasa + m_p[-1] , t_thrust, F)
generar_archivo_eng(Manofacturador, Nombre, D_ext, L_motor, m_p[-1], MasaCarcasa + m_p[-1], t_thrust, F, It)

end_time = time.time()

# Calcular el tiempo transcurrido
elapsed_time = end_time - start_time
print(f"Tiempo transcurrido del programa: {elapsed_time:.2f} segundos")

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
# 20/04/2024
# -La forma en incrementar la masa (Inc_M) no hacia que convergia la solución para diferentes tipos de propelente
# ya ha sido solucionado
# - Se agrego una lista de datos de para el diseño de la tobera la tobera
#
# 30/06/2024
# - Se cambio el tipo de solución de las ecuaciones de la presión de un ciclo a un metodo numerico utilizando solve_ivp para mayor presición,
# bajo mucho el tiempo de calculo 
# - Comparandolo con el excel de Nakka y OpenMotor hay mas correlación con el Excel pienso que puede ser debido a que aun sigo utilizando la 
# forma de calcular la fuerza con el metodo de Nakka debido a que aqui en el empuje es donde mas diferencia hay entre el Excel y OpenMotor
#
#
# 01/07/2024
# Ya se modifico la forma en como se calcula en empuje, y se agregaron las eficiencias N_throat y N_div multiplicadas a Cf 
# esto para dar mas fiabilidad a nuestros calculos
#
#
# 07/07/2024
# Se le agrego como varia la presión, la densidad y la temperatura de la atmosfera con la función Propiedades_Atm(h) para en un futuro poder calcular como es que se adapta 
# la tobera a diferentes alturas con diferentes presiones, talvez sea la version 0.4, pero si lo agrego ya no podria compararlo con OpenMotor o el excel de Nakka
# Por lo que deberia agregar una variable para activar o desactivar esta función.
#
#
# 11/08/2024
# Ya se puede obtener el archivo .eng de los datos de motor, se agrego la funcion de poder activar o desactivar como cambia la presión respecto a la altura que
# este metodo tambien lo implementa OpenRocket, por lo cual puede ser mas fiable, de igual maner al desactivarlo se puede comparar los datos directamente con OpenMotor
# Tambien se puede desactivar el vuelo, o sea se puede tener una simulación estatica como OpenMotor. Solo hay que tener cuidado con el coeficiente de arraste
# se tiene que elegir el valor promedio de este pero aun así no es lo ideal, seria meter como es que varia pero eso talvez es muy complicado aún
#