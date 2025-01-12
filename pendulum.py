import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tkinter as tk
from tkinter import ttk, messagebox

# Ecuaciones de movimiento
def coupled_pendulum_derivatives(state, L, m, k, g):
    theta1, theta2, p1, p2 = state
    delta_theta = theta2 - theta1

    # Aceleraciones angulares
    dtheta1_dt = p1 / (m * L**2)
    dtheta2_dt = p2 / (m * L**2)
    dp1_dt = -m * g * L * np.sin(theta1) - k * delta_theta
    dp2_dt = -m * g * L * np.sin(theta2) + k * delta_theta

    return np.array([dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt])

# Método de integración (RK4)
def integrate(state, dt, L, m, k, g):
    k1 = coupled_pendulum_derivatives(state, L, m, k, g)
    k2 = coupled_pendulum_derivatives(state + dt * k1 / 2, L, m, k, g)
    k3 = coupled_pendulum_derivatives(state + dt * k2 / 2, L, m, k, g)
    k4 = coupled_pendulum_derivatives(state + dt * k3, L, m, k, g)
    return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Simulación principal
def start_simulation():
    try:
        # Leer entradas del usuario
        L = float(length_entry.get())
        m = float(mass_entry.get())
        k = float(spring_entry.get())
        dt = float(timestep_entry.get())
        theta1_0 = float(theta1_entry.get())
        theta2_0 = float(theta2_entry.get())
        T = float(time_entry.get())

        # Constantes
        g = 9.81

        # Validación básica de parámetros
        if L <= 0 or m <= 0 or k < 0 or dt <= 0 or T <= 0:
            messagebox.showerror("Error", "Los valores deben ser positivos y mayores a cero.")
            return

        # Condiciones iniciales
        state = np.array([theta1_0, theta2_0, 0.0, 0.0])
        steps = int(T / dt)

        # Guardado de datos
        theta1_vals, theta2_vals = [], []
        for _ in range(steps):
            state = integrate(state, dt, L, m, k, g)
            theta1_vals.append(state[0])
            theta2_vals.append(state[1])

        # Configuración del gráfico
        fig, ax = plt.subplots()
        ax.set_xlim(-2.5 * L, 2.5 * L)
        ax.set_ylim(-2.5 * L, 2.5 * L)
        ax.set_aspect('equal')
        ax.grid(True)
        ax.set_title("Simulación de Péndulos Acoplados")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        line1, = ax.plot([], [], 'o-', lw=2, color="orange", label="Péndulo 1")
        line2, = ax.plot([], [], 'o-', lw=2, color="blue", label="Péndulo 2")
        spring_line, = ax.plot([], [], 'r--', lw=1, label="Resorte")
        ax.legend()

        # Funciones de animación
        def init():
            line1.set_data([], [])
            line2.set_data([], [])
            spring_line.set_data([], [])
            return line1, line2, spring_line

        def update(frame):
            theta1 = theta1_vals[frame]
            theta2 = theta2_vals[frame]

            # Posición del péndulo 1
            x1, y1 = L * np.sin(theta1), -L * np.cos(theta1)

            # Posición del péndulo 2
            x2, y2 = L * np.sin(theta2) + 2 * L, -L * np.cos(theta2)

            # Actualizar los péndulos
            line1.set_data([0, x1], [0, y1])
            line2.set_data([2 * L, x2], [0, y2])

            # Actualizar resorte
            spring_line.set_data([x1, x2], [y1, y2])
            return line1, line2, spring_line

        ani = FuncAnimation(fig, update, frames=range(steps), init_func=init, blit=True, interval=dt * 1000)
        plt.show()

    except ValueError:
        messagebox.showerror("Error", "Por favor, ingrese valores numéricos válidos.")
    except Exception as e:
        messagebox.showerror("Error inesperado", f"Ocurrió un error: {e}")

# Interfaz gráfica (Tkinter)
root = tk.Tk()
root.title("Simulador de Péndulos Acoplados")
root.geometry("400x400")

# Estilo
ttk.Style().theme_use('clam')
style = ttk.Style()
style.configure("TButton", font=("Helvetica", 12), padding=6)
style.configure("TLabel", font=("Helvetica", 10))
style.configure("TEntry", font=("Helvetica", 10))

# Marco para las entradas
frame = ttk.LabelFrame(root, text="Parámetros de Simulación", padding=(10, 10))
frame.pack(padx=20, pady=20, fill="both", expand=True)

# Espacios para los inputs
inputs = [
    ("Longitud de los péndulos (m):", "1.0"),
    ("Masa de los péndulos (kg):", "1.0"),
    ("Constante del resorte (N/m):", "0.5"),
    ("Paso de tiempo (s):", "0.01"),
    ("Ángulo inicial 1 (θ1) (rad):", "0.2"),
    ("Ángulo inicial 2 (θ2) (rad):", "0.1"),
    ("Tiempo total de la simulación (s):", "10"),
]

entries = []
for i, (label, default) in enumerate(inputs):
    ttk.Label(frame, text=label).grid(row=i, column=0, padx=5, pady=5, sticky="W")
    entry = ttk.Entry(frame)
    entry.insert(0, default)
    entry.grid(row=i, column=1, padx=5, pady=5)
    entries.append(entry)

length_entry, mass_entry, spring_entry, timestep_entry, theta1_entry, theta2_entry, time_entry = entries

# Botón para iniciar la simulación
start_button = ttk.Button(root, text="Iniciar Simulación", command=start_simulation)
start_button.pack(pady=20)

# Ejecutar interfaz gráfica
root.mainloop()