import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import griddata

import sys
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

# Step 1: Read the CSV file
data = pd.read_csv("ORC\\Results\\tab_ORC_T_tap_30_glide_cd_10.csv")

# Remove leading and trailing whitespace from column names
data.columns = data.columns.str.strip()

# Remove rows where W_dot_exp <= 0
# data = data[data['W_dot_exp [W]'] > 0]
data = data[data['eta_ORC [-]'] > 0]
# data = data[data['T_sto_HT [K]']  80]

# Step 2: Extract required columns
T_sto_HT = data['T_sto_HT [K]']
T_tap_in = data['T_tap_in [K]']
glide_ev = data['glide_ev [K]']
W_dot_exp = data['W_dot_exp [W]']
eta = data['eta_ORC [-]']
m_dot_r = data['m_dot_r [kg/s]']


"--------------------------------------------------------------------------------------------------------------"
# Define grid
W_dot_exp_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(W_dot_exp), max(W_dot_exp), 100),
                                             np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# Interpolate eta_ORC onto the grid
eta_grid = griddata((W_dot_exp, T_sto_HT),  eta*100, (W_dot_exp_grid, T_sto_HT_grid), method='linear')

# Define levels for contour lines
eta_levels = np.arange(0, 0.1, 0.01)*100

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(W_dot_exp_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(W_dot_exp_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, colors='black')

# Label the contour lines with inline text
ax.clabel(contour_lines, inline=True, fmt='%d')

# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$\dot{W}_{exp}$ [W]', fontsize=20)
ax.set_ylabel(r'$T_{sto,HT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label(r'$\eta_{ORC}$ [%]', fontsize=16)

# Save the plot as SVG
plt.savefig('ORC\\Graphs\\W_dot_VS_T_sto_HT_VS_eta_T_tap_30_glide_cd_10.svg', format='svg')

# Show plot
plt.show()

"--------------------------------------------------------------------------------------------------------------"

# Define grid
glide_ev_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(glide_ev), max(glide_ev), 100),
                                             np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# Interpolate eta_ORC onto the grid
eta_grid = griddata((glide_ev, T_sto_HT),  eta*100, (glide_ev_grid, T_sto_HT_grid), method='linear')

# Define levels for contour lines
eta_levels = np.arange(0, 0.1, 0.01)*100

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(glide_ev_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(glide_ev_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, colors='black')

# Label the contour lines with inline text
ax.clabel(contour_lines, inline=True, fmt='%d')

# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$glide_{ev}$ [K]', fontsize=20)
ax.set_ylabel(r'$T_{sto,HT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label(r'$\eta_{ORC}$ [%]', fontsize=16)

# Save the plot as SVG
plt.savefig('ORC\\Graphs\\glide_ev_VS_T_sto_HT_VS_eta_T_tap_30_glide_cd_10.svg', format='svg')

# Show plot
plt.show()

"--------------------------------------------------------------------------------------------------------------"

# Define grid
glide_ev_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(glide_ev), max(glide_ev), 100),
                                             np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# Interpolate eta_ORC onto the grid
W_dot_exp_grid = griddata((glide_ev, T_sto_HT),  W_dot_exp, (glide_ev_grid, T_sto_HT_grid), method='linear')

# Define levels for contour lines
W_dot_exp_levels = np.arange(0, 10000, 1000)

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(glide_ev_grid, T_sto_HT_grid, W_dot_exp_grid, levels=W_dot_exp_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(glide_ev_grid, T_sto_HT_grid, W_dot_exp_grid, levels=W_dot_exp_levels, colors='black')

# Label the contour lines with inline text
ax.clabel(contour_lines, inline=True, fmt='%d')

# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$glide_{ev}$ [K]', fontsize=20)
ax.set_ylabel(r'$T_{sto,HT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label(r'$\dot{W}_{exp}$ [W]', fontsize=16)

# Save the plot as SVG
plt.savefig('ORC\\Graphs\\glide_ev_VS_T_sto_HT_VS_W_dot_exp_T_tap_30_glide_cd_10.svg', format='svg')

# Show plot
plt.show()

"--------------------------------------------------------------------------------------------------------------"
# # Define grid
# W_dot_exp_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(W_dot_exp), max(W_dot_exp), 100),
#                                              np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# # Interpolate eta_ORC onto the grid
# eta_grid = griddata((W_dot_exp, T_sto_HT),  eta, (W_dot_exp_grid, T_sto_HT_grid), method='linear')

# # Define levels for contour lines
# eta_levels = np.arange(0, 0.1, 0.01)

# # Create a new figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Create a 2D contour plot with specified levels
# contour = ax.contourf(W_dot_exp_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, cmap='viridis')

# # Add contour lines
# contour_lines = ax.contour(W_dot_exp_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, colors='black')

# # Label the contour lines with inline text
# ax.clabel(contour_lines, inline=True, fmt='%.2f')

# # Set axis labels with LaTeX formatting and bigger font
# ax.set_xlabel(r'$\dot{W}_{exp}$ [W]', fontsize=20)
# ax.set_ylabel(r'$T_{sto,HT}$ [°C]', fontsize=20)

# # Remove graph title
# ax.set_title('')

# # Add grid
# ax.grid(True)

# # Save the plot as SVG
# plt.savefig('ORC\\Graphs\\W_dot_VS_T_sto_HT_VS_eta_glide_15.svg', format='svg')

# # Show plot
# plt.show()


# #---------------------------------------------------------------------------------------------------

# T_tap_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(T_tap_in), max(T_tap_in), 100),
#                                              np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# # Interpolate eta_ORC onto the grid
# eta_grid = griddata((T_tap_in, T_sto_HT),  eta, (T_tap_grid, T_sto_HT_grid), method='linear')

# # Define levels for contour lines
# eta_levels = np.arange(0, 0.1, 0.01)

# # Create a new figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Create a 2D contour plot with specified levels
# contour = ax.contourf(T_tap_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, cmap='viridis')

# # Add contour lines
# contour_lines = ax.contour(T_tap_grid, T_sto_HT_grid, eta_grid, levels=eta_levels, colors='black')

# # Label the contour lines with inline text
# ax.clabel(contour_lines, inline=True, fmt='%.2f')

# # Set axis labels with LaTeX formatting and bigger font
# ax.set_xlabel(r'$T_{tap}$ [°C]', fontsize=20)
# ax.set_ylabel(r'$T_{sto,HT}$ [°C]', fontsize=20)

# # Remove graph title
# ax.set_title('')

# # Add grid
# ax.grid(True)

# # Save the plot as SVG
# plt.savefig('ORC\\Graphs\\T_tap_VS_T_sto_HT_VS_eta_glide_15.svg', format='svg')

# # Show plot
# plt.show()