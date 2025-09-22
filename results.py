import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

cpp_file = "q1_main.cpp"  
exe_file = "text"          

compile_result = subprocess.run(["g++", "-o", exe_file, cpp_file])

if compile_result.returncode != 0:
    print("Compilation failed.")
    sys.exit(1)

run_result = subprocess.run([f"./{exe_file}"], capture_output=True, text=True)

if run_result.returncode != 0:
    print("Execution failed.")
    sys.exit(1)

try:
    output_data = run_result.stdout

except Exception as e:
    print(f"try again ")
    sys.exit(1)

blocks = [block for block in output_data.strip().split('\n\n') if block]
if not blocks:
    print("No output")
    sys.exit(1)


try:
    convergence_blocks = blocks[:3]
    # plot1_data = blocks[0]
    # converged_points = {'x': [], 'y': []}
    # not_converged_points = {'x': [], 'y': []}
    # lines=plot1_data.strip().split('\n')

    # all_x_values=[]
    # all_y_values=[]

    # for line in lines:
    #     parts = line.split()
    #     omega = float(parts[0])
    #     iterations = int(parts[1])
    #     status = int(parts[2])
    #     all_x_values.append(omega)
    #     all_y_values.append(iterations)
    #     if status == 0:
    #         converged_points['x'].append(omega)
    #         converged_points['y'].append(iterations)
    #     else:
    #         not_converged_points['x'].append(omega)
    #         not_converged_points['y'].append(iterations)

    # fig1, ax1 = plt.subplots()
    # ax1.plot(all_x_values, all_y_values,label='Point GS')
    # ax1.plot(converged_points['x'], converged_points['y'],color='green', marker='o', linestyle='None', label='Converged')
    # ax1.plot(not_converged_points['x'], not_converged_points['y'],color='red', marker='d', linestyle='None', label='Not Converged')
        
    fig1, ax1 = plt.subplots(figsize=(10, 7))
    colors = ['blue', 'pink', 'black'] # Colors for Point GS, Line GS, ADI

    for i, block in enumerate(convergence_blocks):
        lines = block.strip().split('\n')
        method_name = lines[0].strip()
        print(method_name)
        data_lines = lines[1:]

        # Prepare lists to hold data for the current method
        all_x_values = []
        all_y_values = []
        converged_points = {'x': [], 'y': []}
        not_converged_points = {'x': [], 'y': []}

        for line in data_lines:
            parts = line.split()
            omega = float(parts[0])
            iterations = int(parts[1])
            status = int(parts[2])

            all_x_values.append(omega)
            all_y_values.append(iterations)

            if status == 0: # Converged
                converged_points['x'].append(omega)
                converged_points['y'].append(iterations)
            else: # Not Converged
                not_converged_points['x'].append(omega)
                not_converged_points['y'].append(iterations)

        # Plotting for the current method
        color = colors[i % len(colors)]
        ax1.plot(all_x_values, all_y_values, color=color, linestyle='-', label=method_name)
        ax1.plot(converged_points['x'], converged_points['y'], color='green', marker='o', linestyle='None')
        ax1.plot(not_converged_points['x'], not_converged_points['y'], color='red', marker='x', linestyle='None', markersize=8)

    # Add dummy plots for a clean legend explaining the markers
    ax1.plot([], [], color='green', marker='o', linestyle='None', label='Converged')
    ax1.plot([], [], color='red', marker='x', linestyle='None', label='Not Converged')

    ax1.set_xlabel('Overrelaxation Factor (ω)')
    ax1.set_ylabel('Total Number of Iterations')
    ax1.set_title('Convergence Analysis Comparison')
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax1.legend()

    # ax1.set_xlabel('Overrelaxation value')
    # ax1.set_ylabel('Total Number of Iterations')
    # ax1.set_title('Convergence Analysis in point GS SOR')
    # ax1.grid()
    # ax1.legend()

except (ValueError, IndexError) as e:
    print(f"Error: {e}")

# this is code for drawing the countors for tempsss
levels=50
# try:
#     data_grid = []
#     output_data=blocks[1]
#     lines = output_data.strip().split('\n')
#     for line in lines:
#         row = [float(value) for value in line.split()]
#         data_grid.append(row)

#     temperature_data = np.array(data_grid)

# except (ValueError, IndexError) as e:
#     sys.exit(1)

# if temperature_data.size == 0:
#     sys.exit(1)

# fig, ax = plt.subplots()

# contour_fill = ax.contourf(temperature_data, levels=levels, cmap='plasma')
# fig.colorbar(contour_fill, ax=ax, label='Temperature (°C)')
# ax.set_title('2D Temperature Distribution')
# ax.set_xlabel('X-coordinate')
# ax.set_ylabel('Y-coordinate')
# ax.set_aspect('equal', adjustable='box')
# fig.savefig('temperature_contour.png', dpi=300)
# plt.show()

residual_blocks = blocks[3:]
if residual_blocks:
    try:
        fig2, ax2 = plt.subplots()
        colors = ["blue","orange","pink"]
        markers = ['o', 'x', 's']

        for i, block in enumerate(residual_blocks):
            lines = block.strip().split('\n')
            legend_label = lines[0].strip()
            # print(legend_label)
            data_lines = lines[1:]

            iters = [int(line.split()[0]) for line in data_lines]
            residuals = [float(line.split()[1]) for line in data_lines]
            
            ax2.plot(iters, residuals, label=legend_label,color=colors[i],marker=markers[i])

        ax2.set_yscale('log')
        ax2.set_xlabel('Iteration Number')
        ax2.set_ylabel('Residual (Log Scale)')
        ax2.set_title('Residual Reduction Comparison')
        ax2.grid()
        ax2.legend()

    except (ValueError, IndexError) as e:
        print(f"Error: {e}")

plt.tight_layout()
fig1.savefig('omega_vs_iterations_non_sym.png', dpi=300)
fig2.savefig('residual_comparison_non_sym_omega=1.2.png', dpi=300)
plt.show()