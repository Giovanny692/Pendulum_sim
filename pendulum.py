import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tkinter as tk
from tkinter import ttk

# Function to start the simulation
def start_simulation():
    # Retrieve user inputs
    L = float(length_entry.get())
    m = float(mass_entry.get())
    k = float(spring_entry.get())
    dt = float(timestep_entry.get())
    theta1_0 = float(theta1_entry.get())
    theta2_0 = float(theta2_entry.get())
    T = float(time_entry.get())

    # Constants
    g = 9.81  # Gravitational acceleration (m/s^2)

    # Equations of Motion
    def coupled_pendulum_derivatives(state):
        theta1, theta2, p1, p2 = state
        delta_theta = theta2 - theta1

        # Angular accelerations
        dtheta1_dt = p1 / (m * L**2)
        dtheta2_dt = p2 / (m * L**2)
        dp1_dt = -m * g * L * np.sin(theta1) - k * delta_theta
        dp2_dt = -m * g * L * np.sin(theta2) + k * delta_theta

        return np.array([dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt])

    # Numerical Integration (Runge-Kutta 4th order)
    def integrate(state, dt):
        k1 = coupled_pendulum_derivatives(state)
        k2 = coupled_pendulum_derivatives(state + dt * k1 / 2)
        k3 = coupled_pendulum_derivatives(state + dt * k2 / 2)
        k4 = coupled_pendulum_derivatives(state + dt * k3)
        return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Initial Conditions
    state = np.array([theta1_0, theta2_0, 0.0, 0.0])
    steps = int(T / dt)

    # Data Storage
    theta1_vals, theta2_vals = [], []
    for _ in range(steps):
        state = integrate(state, dt)
        theta1_vals.append(state[0])
        theta2_vals.append(state[1])

    # Visualization with Matplotlib
    fig, ax = plt.subplots()
    ax.set_xlim(-2 * L, 2 * L)
    ax.set_ylim(-2 * L, 2 * L)
    ax.set_aspect('equal')
    line1, = ax.plot([], [], 'o-', lw=2, color="orange")
    line2, = ax.plot([], [], 'o-', lw=2, color="orange")
    spring_line, = ax.plot([], [], 'r-', lw=1)

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        spring_line.set_data([], [])
        return line1, line2, spring_line

    def update(frame):
        theta1 = theta1_vals[frame]
        theta2 = theta2_vals[frame]

        # Pendulum 1 position
        x1, y1 = L * np.sin(theta1), -L * np.cos(theta1)

        # Pendulum 2 position
        x2, y2 = L * np.sin(theta2) + 2 * L, -L * np.cos(theta2)

        # Update pendulums
        line1.set_data([0, x1], [0, y1])
        line2.set_data([2 * L, x2], [0, y2])

        # Update spring
        spring_line.set_data([x1, x2], [y1, y2])
        return line1, line2, spring_line

    ani = FuncAnimation(fig, update, frames=range(steps), init_func=init, blit=True, interval=dt*1000)
    plt.show()

# Create the Tkinter window
root = tk.Tk()
root.title("Coupled Pendulum Simulator")

# Create input fields
ttk.Label(root, text="Length of the pendulums (m):").grid(row=0, column=0, padx=5, pady=5)
length_entry = ttk.Entry(root)
length_entry.insert(0, "1.0")
length_entry.grid(row=0, column=1, padx=5, pady=5)

ttk.Label(root, text="Mass of the pendulums (kg):").grid(row=1, column=0, padx=5, pady=5)
mass_entry = ttk.Entry(root)
mass_entry.insert(0, "1.0")
mass_entry.grid(row=1, column=1, padx=5, pady=5)

ttk.Label(root, text="Spring constant (N/m):").grid(row=2, column=0, padx=5, pady=5)
spring_entry = ttk.Entry(root)
spring_entry.insert(0, "0.5")
spring_entry.grid(row=2, column=1, padx=5, pady=5)

ttk.Label(root, text="Time step (s):").grid(row=3, column=0, padx=5, pady=5)
timestep_entry = ttk.Entry(root)
timestep_entry.insert(0, "0.01")
timestep_entry.grid(row=3, column=1, padx=5, pady=5)

ttk.Label(root, text="Initial angle (θ1) (rad):").grid(row=4, column=0, padx=5, pady=5)
theta1_entry = ttk.Entry(root)
theta1_entry.insert(0, "0.2")
theta1_entry.grid(row=4, column=1, padx=5, pady=5)

ttk.Label(root, text="Initial angle (θ2) (rad):").grid(row=5, column=0, padx=5, pady=5)
theta2_entry = ttk.Entry(root)
theta2_entry.insert(0, "0.1")
theta2_entry.grid(row=5, column=1, padx=5, pady=5)

ttk.Label(root, text="Total simulation time (s):").grid(row=6, column=0, padx=5, pady=5)
time_entry = ttk.Entry(root)
time_entry.insert(0, "10")
time_entry.grid(row=6, column=1, padx=5, pady=5)

# Create start button
start_button = ttk.Button(root, text="Start Simulation", command=start_simulation)
start_button.grid(row=7, column=0, columnspan=2, pady=10)

# Run the Tkinter event loop
root.mainloop()
