import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from matplotlib.patches import Circle  # Import Circle patch

# Load data from a .dat file using numpy
datafile = "output.dat"  # Replace with your actual file path
data = np.loadtxt(datafile, skiprows=10)

# Create a figure and axis
fig, ax = plt.subplots()

# Set x and y axis limits based on the data
x_min, x_max = data[:, 1].min(), data[:, 1].max()
y_min, y_max = data[:, 2].min(), data[:, 2].max()

tolerance = 0.1
max_size = max(x_max + tolerance * (x_max - x_min), y_max + tolerance * (y_max - y_min))
min_size = min(x_min - tolerance * (x_max - x_min), y_min - tolerance * (y_max - y_min))

ax.set_xlim(min_size , max_size)
ax.set_ylim(min_size, max_size)

# Add a circle centered at (0, 0) with radius M
M = 1  # Set the desired radius
circle = Circle((0, 0), radius=2*M, edgecolor='none', facecolor='black')  # Red circle with no fill
ax.add_patch(circle)

# Line that will be updated in the animation
E = 0.1
L_z = 2.4
M = 1.0
line, = ax.plot([], [], "bo", linewidth=None, ms=0.2)
ax.set_title(f"Schwaszchild solver for a massive particle with E = {E}, L_z = {L_z} and M = {M}")

# Function to initialize the plot (for the first frame)
def init():
    line.set_data([], [])
    return line,

# Function to update the plot for each frame
def update(frame):
    # Extract x and y values up to the current frame (row)
    x = data[:frame, 1]  # X-values (first column)
    y = data[:frame, 2]  # Y-values (second column)
    
    # Set the data for the line
    line.set_data(x, y)
    
    if frame - 1 > 0:
        line.set_label(f'Delta proper time : {str((data[frame - 1, 0] - data[frame-2, 0])*100)[:5]}')  # Update the label of the line
        ax.legend(loc='upper right')
    
    return line,

# Create the animation
step = int(len(data) / 500)
ani = FuncAnimation(fig, update, frames=range(1, len(data) + 1, step),
                    init_func=init, interval=20, repeat=True)

# Save the animation as a GIF
writer = PillowWriter(fps=50)
ani.save("massive_falling.gif", writer=writer)

plt.show()
