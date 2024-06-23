import matplotlib.pyplot as plt
import numpy as np

def draw_fenocyl_grain(ax, Ri, Re, h, N, w, d, rb, steps):
    """
    Dibuja el grano de propelente tipo fenocyl con contornos que representan la quema a lo largo del tiempo.

    Parámetros:
    ax: Axes object de matplotlib para dibujar
    Ri (float): Radio interno inicial del grano (cm)
    Re (float): Radio externo inicial del grano (cm)
    h (float): Altura del grano (cm)
    N (int): Número de aletas
    w (float): Ancho de la aleta (cm)
    d (float): Profundidad de la aleta (cm)
    rb (float): Tasa de quemado (cm/s)
    steps (int): Número de pasos para dibujar los contornos
    """
    
    theta = np.linspace(0, 2 * np.pi, 100)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    
    for t in range(steps):
        tiempo = t * (Re - Ri) / (rb * steps)
        Ri_t = Ri + rb * tiempo
        Re_t = Re - rb * tiempo
        
        if Ri_t >= Re_t:
            break
        
        x_outer = Re_t * np.cos(theta)
        y_outer = Re_t * np.sin(theta)
        
        x_inner = Ri_t * np.cos(theta)
        y_inner = Ri_t * np.sin(theta)
        
        ax.plot(x_outer, y_outer, 'k')
        ax.plot(x_inner, y_inner, 'k')
        
        for angle in angles:
            # Coordenadas del inicio y fin de la aleta en el ángulo dado
            x_start = Ri_t * np.cos(angle)
            y_start = Ri_t * np.sin(angle)
            x_end = (Ri_t + d - rb * tiempo) * np.cos(angle)
            y_end = (Ri_t + d - rb * tiempo) * np.sin(angle)
            
            # Puntos para formar la aleta rectangular centrada
            x_rect = [x_start - w/2 * np.sin(angle), 
                      x_start + w/2 * np.sin(angle), 
                      x_end + w/2 * np.sin(angle), 
                      x_end - w/2 * np.sin(angle)]
            
            y_rect = [y_start + w/2 * np.cos(angle), 
                      y_start - w/2 * np.cos(angle), 
                      y_end - w/2 * np.cos(angle), 
                      y_end + w/2 * np.cos(angle)]
            
            # Dibujar las aletas
            ax.plot([x_rect[0], x_rect[1]], [y_rect[0], y_rect[1]], 'k')
            ax.plot([x_rect[1], x_rect[2]], [y_rect[1], y_rect[2]], 'k')
            ax.plot([x_rect[2], x_rect[3]], [y_rect[2], y_rect[3]], 'k')
            ax.plot([x_rect[3], x_rect[0]], [y_rect[3], y_rect[0]], 'k')

fig, ax = plt.subplots(figsize=(6, 6))
draw_fenocyl_grain(ax, Ri=1.0, Re=5.0, h=10.0, N=6, w=.5, d=2.0, rb=0.1, steps=50)
ax.set_aspect('equal')
ax.set_axis_off()
plt.show()






