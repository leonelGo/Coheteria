import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np

# Función para graficar
def graficar():
    # Obtener los valores ingresados
    a = float(entrada_a.get())
    b = float(entrada_b.get())
    
    # Crear datos para la gráfica (puedes personalizarlo)
    x = np.linspace(0, 10, 100)
    y = a * np.sin(b * x)
    
    # Limpiar el gráfico anterior
    ax.clear()
    
    # Graficar los nuevos datos
    ax.plot(x, y, label=f'{a} * sin({b} * x)')
    ax.legend()
    
    # Actualizar el canvas
    canvas.draw()

# Crear la ventana principal
ventana = tk.Tk()
ventana.title("Interfaz de Python con Gráficas y Múltiples Entradas")
ventana.geometry("800x600")

# Crear el frame para la entrada de datos
frame_entrada = ttk.Frame(ventana)
frame_entrada.pack(pady=20)

# Etiqueta y campo de entrada para el valor 'a'
ttk.Label(frame_entrada, text="Ingrese el valor de 'a':").pack(side=tk.LEFT, padx=5)
entrada_a = ttk.Entry(frame_entrada, width=10)
entrada_a.pack(side=tk.LEFT, padx=5)

# Etiqueta y campo de entrada para el valor 'b'
ttk.Label(frame_entrada, text="Ingrese el valor de 'b':").pack(side=tk.LEFT, padx=5)
entrada_b = ttk.Entry(frame_entrada, width=10)
entrada_b.pack(side=tk.LEFT, padx=5)

# Botón para graficar
boton_graficar = ttk.Button(frame_entrada, text="Graficar", command=graficar)
boton_graficar.pack(side=tk.LEFT, padx=15)

# Crear el frame para la gráfica
frame_grafica = ttk.Frame(ventana)
frame_grafica.pack()

# Crear la figura de matplotlib
fig, ax = plt.subplots()

# Crear el canvas de tkinter para matplotlib
canvas = FigureCanvasTkAgg(fig, master=frame_grafica)
canvas.get_tk_widget().pack()

# Iniciar el bucle principal de la interfaz
ventana.mainloop()
