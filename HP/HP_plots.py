import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import griddata

import sys
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

# Step 1: Read the CSV file
data = pd.read_csv("HP\\Results\\tab_HP_TDH_40_glide_ev_10.csv")

# Remove leading and trailing whitespace from column names
data.columns = data.columns.str.strip()

# Remove rows where W_dot_exp <= 0
# data = data[data['W_dot_exp [W]'] > 0]
data = data[data['COP [-]'] > 0]
# data = data[data['T_DH+ [K]'] == 50]

# Step 2: Extract required columns
T_sto_HT = data['T_sto_HT [K]']
T_sto_LT = data['T_sto_LT [K]']
glide_cd = data['glide_cd [K]']
T_su_DH = data['T_DH+ [K]']
W_dot_comp = data['W_dot_cp [W]']
COP = data['COP [-]']
m_dot_r = data['m_dot_r [kg/s]']

"--------------------------------------------------------------------------------------------------"
"Plot for T_su_DH imposed"
# Assuming glide_cd, T_sto_LT, and COP are already defined arrays

# Define grid
glide_cd_grid, T_sto_LT_grid = np.meshgrid(np.linspace(min(glide_cd), max(glide_cd), 100),
                                           np.linspace(min(T_sto_LT), max(T_sto_LT), 100))

# Interpolate COP onto the grid
COP_grid = griddata((glide_cd, T_sto_LT), COP, (glide_cd_grid, T_sto_LT_grid), method='linear')

# Define levels for contour lines
COP_levels = np.arange(0, 20, 1)

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(glide_cd_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(glide_cd_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, colors='black')

# Label the contour lines with inline text, without decimal places
ax.clabel(contour_lines, inline=True, fmt='%d')

# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$glide_{cd}$ [K]', fontsize=20)
ax.set_ylabel(r'$T_{sto,LT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label('COP', fontsize=16)

# Save the plot as SVG
plt.savefig('HP/plot/glide_cd_VS_T_sto_LT_VS_COP_TDH_40_glide_ev_10.svg', format='svg')

# Show plot
plt.show()


#---------------------------------------------------------------------------------------------------

# Define grid
glide_cd_grid, T_sto_LT_grid = np.meshgrid(np.linspace(min(glide_cd), max(glide_cd), 100),
                                             np.linspace(min(T_sto_LT), max(T_sto_LT), 100))

# Interpolate eta_ORC onto the grid
W_dot_comp_grid = griddata((glide_cd, T_sto_LT),  W_dot_comp, (glide_cd_grid, T_sto_LT_grid), method='linear')

# Define levels for contour lines
W_dot_comp_levels = np.arange(0, 18000, 1000)

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(glide_cd_grid, T_sto_LT_grid, W_dot_comp_grid, levels=W_dot_comp_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(glide_cd_grid, T_sto_LT_grid, W_dot_comp_grid, levels=W_dot_comp_levels, colors='black')

# Label the contour lines with inline text, without decimal places
ax.clabel(contour_lines, inline=True, fmt='%d')


# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$glide_{cd}$ [K]', fontsize=20)
ax.set_ylabel(r'$T_{sto,LT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label(r'$\dot{W}_{comp}$ [W]', fontsize=16)

# Save the plot as SVG
plt.savefig('HP\\plot\\glide_cd_VS_T_sto_LT_VS_W_TDH_40_glide_ev_10.svg', format='svg')

# Show plot
plt.show()

#-------------------------------------------------------------------------------------------
# Define grid
W_dot_comp_grid, T_sto_LT_grid = np.meshgrid(np.linspace(min(W_dot_comp), max(W_dot_comp), 100),
                                             np.linspace(min(T_sto_LT), max(T_sto_LT), 100))

# Interpolate eta_ORC onto the grid
COP_grid = griddata((W_dot_comp, T_sto_HT),  COP, (W_dot_comp_grid, T_sto_LT_grid), method='linear')

# Define levels for contour lines
COP_levels = np.arange(0, 20, 1)

# Create a new figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create a 2D contour plot with specified levels
contour = ax.contourf(W_dot_comp_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, cmap='viridis')

# Add contour lines
contour_lines = ax.contour(W_dot_comp_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, colors='black')

# Label the contour lines with inline text
ax.clabel(contour_lines, inline=True, fmt='%d')

# Set axis labels with LaTeX formatting and bigger font
ax.set_xlabel(r'$\dot{W}_{comp}$ [W]', fontsize=20)
ax.set_ylabel(r'$T_{sto,LT}$ [°C]', fontsize=20)

# Remove graph title
ax.set_title('')

# Add grid
ax.grid(True)

# Add a color bar to the right side of the graph
colorbar = plt.colorbar(contour, ax=ax)
colorbar.set_label('COP', fontsize=16)

# Save the plot as SVG
plt.savefig('HP\\plot\\W_dot_VS_T_sto_LT_VS_COP_TDH_40_glide_ev_10.svg', format='svg')

# Show plot
plt.show()
# "--------------------------------------------------------------------------------------------------"

# # Define grid
# W_dot_comp_grid, T_su_DH_grid = np.meshgrid(np.linspace(min(W_dot_comp), max(W_dot_comp), 100),
#                                              np.linspace(min(T_su_DH), max(T_su_DH), 100))

# # Interpolate eta_ORC onto the grid
# COP_grid = griddata((W_dot_comp, T_su_DH),  COP, (W_dot_comp_grid, T_su_DH_grid), method='linear')

# # Define levels for contour lines
# COP_levels = np.arange(0, 17, 2)

# # Create a new figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Create a 2D contour plot with specified levels
# contour = ax.contourf(W_dot_comp_grid, T_su_DH_grid, COP_grid, levels=COP_levels, cmap='viridis')

# # Add contour lines
# contour_lines = ax.contour(W_dot_comp_grid, T_su_DH_grid, COP_grid, levels=COP_levels, colors='black')

# # Label the contour lines with inline text
# ax.clabel(contour_lines, inline=True, fmt='%.2f')

# # Set axis labels with LaTeX formatting and bigger font
# ax.set_xlabel(r'$\dot{W}_{comp}$ [W]', fontsize=20)
# ax.set_ylabel(r'$T_{DH+}$', fontsize=20)

# # Remove graph title
# ax.set_title('')

# # Add grid
# ax.grid(True)

# # Save the plot as SVG
# plt.savefig('HP\\plot\\W_dot_VS_T_DH_VS_COP_glide_15.svg', format='svg')

# # Show plot
# plt.show()

# "--------------------------------------------------------------------------------------"

# # Define grid
# W_dot_comp_grid, T_sto_LT_grid = np.meshgrid(np.linspace(min(W_dot_comp), max(W_dot_comp), 100),
#                                              np.linspace(min(T_sto_LT), max(T_sto_LT), 100))

# # Interpolate eta_ORC onto the grid
# COP_grid = griddata((W_dot_comp, T_sto_LT),  COP, (W_dot_comp_grid, T_sto_LT_grid), method='linear')

# # Define levels for contour lines
# COP_levels = np.arange(0, 17, 2)

# # Create a new figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Create a 2D contour plot with specified levels
# contour = ax.contourf(W_dot_comp_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, cmap='viridis')

# # Add contour lines
# contour_lines = ax.contour(W_dot_comp_grid, T_sto_LT_grid, COP_grid, levels=COP_levels, colors='black')

# # Label the contour lines with inline text
# ax.clabel(contour_lines, inline=True, fmt='%.2f')

# # Set axis labels with LaTeX formatting and bigger font
# ax.set_xlabel(r'$\dot{W}_{comp}$ [W]', fontsize=20)
# ax.set_ylabel(r'$T_{sto,LT}$', fontsize=20)

# # Remove graph title
# ax.set_title('')

# # Add grid
# ax.grid(True)

# # Save the plot as SVG
# plt.savefig('HP\\plot\\W_dot_VS_T_sto_LT_VS_COP_glide_15.svg', format='svg')

# # Show plot
# plt.show()

# #---------------------------------------------------------------------------------------------------

# T_su_DH_grid, T_sto_HT_grid = np.meshgrid(np.linspace(min(T_su_DH), max(T_su_DH), 100),
#                                              np.linspace(min(T_sto_HT), max(T_sto_HT), 100))

# # Interpolate eta_ORC onto the grid
# COP_grid = griddata((T_su_DH, T_sto_HT),  COP, (T_su_DH_grid, T_sto_HT_grid), method='linear')

# # Define levels for contour lines
# COP_levels = np.arange(0, 17, 1)

# # Create a new figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Create a 2D contour plot with specified levels
# contour = ax.contourf(T_su_DH_grid, T_sto_HT_grid, COP_grid, levels=COP_levels, cmap='viridis')

# # Add contour lines
# contour_lines = ax.contour(T_su_DH_grid, T_sto_HT_grid, COP_grid, levels=COP_levels, colors='black')

# # Label the contour lines with inline text
# ax.clabel(contour_lines, inline=True, fmt='%.2f')

# # Set axis labels with LaTeX formatting and bigger font
# ax.set_xlabel(r'$T_{DH,+}$', fontsize=20)
# ax.set_ylabel(r'$T_{sto,HT}$', fontsize=20)

# # Remove graph title
# ax.set_title('')

# # Add grid
# ax.grid(True)

# # Save the plot as SVG
# plt.savefig('HP\\plot\\T_DH_VS_T_sto_HT_VS_COP_glide_15.svg', format='svg')

# # Show plot
# plt.show()