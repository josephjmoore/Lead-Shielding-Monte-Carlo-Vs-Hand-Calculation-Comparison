import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from params import (mu, K_dot_no_shield, gp_buildup, Gy_to_uGy, r_detector,
                    Shield_Material, results_dir)

#Results for the chosen material (set Shield_Material in params.py)
csv_path = results_dir / f"{Shield_Material}_openmc_sweep.csv"
csv = pd.read_csv(csv_path)

#Infinite-medium results (only plotted if that sweep has been run for this material)
csv_inf_path = results_dir / f"{Shield_Material}_openmc_infinite_sweep.csv"
csv_inf = pd.read_csv(csv_inf_path) if csv_inf_path.exists() else None

material = Shield_Material.lower()
print(f'Shield material: {Shield_Material}')

shield_thickness = csv["thickness_cm"]
K_total_openmc_uGy_h = csv["K_total_openmc_uGy_h"]
K_total_openmc_err_uGy_h = csv["K_total_openmc_err_uGy_h"]
K_uncollided_openmc_uGy_h = csv["K_uncollided_openmc_uGy_h"]
K_uncollided_openmc_err_uGy_h = csv["K_uncollided_openmc_err_uGy_h"]
B_openmc = csv["B_openmc"]
B_openmc_err = csv["B_openmc_err"]

#Depths in mean free paths and hand-calc values, recalculated from params.py
#(so they always use this material's current mu and GP coefficients, not values saved in the CSV)
mfp = mu * shield_thickness
B_hand = gp_buildup(mfp)
K_uncollided_hand_uGy_h = K_dot_no_shield * np.exp(-mfp) * Gy_to_uGy
K_total_hand_uGy_h = K_uncollided_hand_uGy_h * B_hand

#Hand calculation as a smooth curve (it's a formula, so evaluate it on a fine grid)
t_max = shield_thickness.max()
if csv_inf is not None:
    t_max = max(t_max, csv_inf["thickness_cm"].max())
t_fine = np.linspace(0, t_max, 200) #cm
mfp_fine = mu * t_fine
B_hand_fine = gp_buildup(mfp_fine)
#Air kerma curve only over the thicknesses simulated with the finite shell
t_fine_shell = t_fine[t_fine <= shield_thickness.max()]
K_hand_fine = K_dot_no_shield * np.exp(-mu * t_fine_shell) * gp_buildup(mu * t_fine_shell) * Gy_to_uGy

#Ratio of OpenMC to hand calculation at each simulated thickness
ratio_total = K_total_openmc_uGy_h / K_total_hand_uGy_h
ratio_total_err = K_total_openmc_err_uGy_h / K_total_hand_uGy_h
ratio_uncollided = K_uncollided_openmc_uGy_h / K_uncollided_hand_uGy_h
ratio_uncollided_err = K_uncollided_openmc_err_uGy_h / K_uncollided_hand_uGy_h

#Colours
colour_openmc = "#83d62a"
colour_hand = '#eb6834'
colour_total = "#ff00fb"
colour_uncollided = '#3d3d3a'
colour_reference = '#8a8a85'
colour_grid = '#d9d9d6'
colour_text = '#3d3d3a'

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.5))

#Left: air kerma rate against shield thickness
ax1.plot(t_fine_shell, K_hand_fine,
         color=colour_hand, linewidth=2, linestyle='--', label='Hand Calculation')
ax1.errorbar(shield_thickness, K_total_openmc_uGy_h, yerr=K_total_openmc_err_uGy_h,
             fmt='x', color=colour_openmc, markersize=6, capsize=3, label='OpenMC')
ax1.set_yscale('log')
ax1.set_xlabel(f'{Shield_Material} thickness (cm)')
ax1.set_ylabel(f'Air kerma rate at {r_detector / 100:g} m (µGy/h)')
ax1.set_title(f'Cs-137 (1 GBq) behind {material}', loc='left')
ax1.legend(frameon=False)

#Middle: buildup factor against depth in mean free paths
ax2.plot(mfp_fine, B_hand_fine,
         color=colour_hand, linewidth=2, linestyle='--', label='Hand Calculation (GP)')
ax2.errorbar(mfp, B_openmc, yerr=B_openmc_err,
             fmt='x', color=colour_openmc, markersize=6, capsize=3, label=f'OpenMC ({material} shell)')
if csv_inf is not None:
    #Depth recalculated from mu, as for the shell results
    mfp_inf = mu * csv_inf["thickness_cm"]
    ax2.errorbar(mfp_inf, csv_inf["B_openmc"], yerr=csv_inf["B_openmc_err"],
                 fmt='x', color=colour_text, markersize=6, capsize=3, label=f'OpenMC (infinite {material})')
ax2.set_xlabel(f'{Shield_Material} thickness (mean free paths)')
ax2.set_ylabel('Exposure buildup factor B')
ax2.set_title(f'Buildup factor, {material}, 0.662 MeV', loc='left')
ax2.set_ylim(bottom=1.0)
ax2.legend(frameon=False)

#Right: ratio of OpenMC to hand calculation
ax3.axhline(1.0, color=colour_reference, linewidth=1, linestyle='--')
ax3.errorbar(shield_thickness, ratio_uncollided, yerr=ratio_uncollided_err,
             fmt='x', color=colour_uncollided, markerfacecolor='none', markersize=6,
             capsize=3, label='Uncollided')
ax3.errorbar(shield_thickness, ratio_total, yerr=ratio_total_err,
             fmt='x', color=colour_total, markersize=6, capsize=3, label='Total')
ax3.set_xlabel(f'{Shield_Material} thickness (cm)')
ax3.set_ylabel('OpenMC / Hand Calculation')
ax3.set_title('OpenMC relative to hand calculation', loc='left')
ax3.legend(frameon=False)

#Keep gridlines and axes recessive so the data stands out
for ax in (ax1, ax2, ax3):
    ax.grid(True, which='major', color=colour_grid, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(colour_grid)
    ax.tick_params(colors=colour_text)
    ax.xaxis.label.set_color(colour_text)
    ax.yaxis.label.set_color(colour_text)
    ax.title.set_color(colour_text)

fig.tight_layout()

#Save the figure in this material's Results folder
fig_path = results_dir / f'{Shield_Material}_openmc_vs_hand_calc.png'
fig.savefig(fig_path, dpi=200)
print(f'Plot saved to {fig_path}')

plt.show()