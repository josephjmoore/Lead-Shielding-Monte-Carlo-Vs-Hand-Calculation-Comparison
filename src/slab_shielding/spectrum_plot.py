from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from params import E, air_energy, air_mu_en_rho, loglog_interp

#Lead thicknesses to plot (cm); any without a spectrum file are skipped
plot_thickness = [1.0, 3.0, 5.0, 7.0]

spectrum_dir = Path(__file__).resolve().parent / 'Results' / 'spectra'

#Features to mark on the plot, in keV
electron_rest_energy = 0.511 #MeV
E_backscatter = E / (1 + 2 * E / electron_rest_energy) * 1e3 #keV, single 180 degree Compton scatter
#(energy, label); features closer together than the labels can fit share one label
features = [
    (E * 1e3, 'Uncollided (662 keV)'),
    (E_backscatter, f'Backscatter edge ({E_backscatter:.0f} keV)'),
    (88.0, 'Pb K-edge (88 keV)'),
    (84.9, None),  #Pb K-beta X-ray
    (75.0, 'Pb K X-rays (73, 75, 85 keV)'),  #Pb K-alpha1 X-ray
    (72.8, None),  #Pb K-alpha2 X-ray
]

#Colours (fixed order, colour-blind safe)
colours = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
colour_feature = '#8a8a85'
colour_grid = '#d9d9d6'
colour_text = '#3d3d3a'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

print('Thickness (cm)   flux ratio (total/peak bin)   air kerma ratio (total/peak bin)')

i = 0
for t in plot_thickness:
    spectrum_file = spectrum_dir / f'spectrum_t_{t}cm.csv'
    if not spectrum_file.exists():
        print(f'{t} cm: no spectrum file ({spectrum_file.name}) - skipped')
        continue
    spectrum = pd.read_csv(spectrum_file)
    colour = colours[i % len(colours)]
    i += 1

    E_low = spectrum['E_low_eV'].to_numpy()
    E_high = spectrum['E_high_eV'].to_numpy()
    flux_per_bin = spectrum['flux_per_bin'].to_numpy() #photons/cm^2/s in each bin
    edges_keV = np.append(E_low, E_high[-1]) / 1e3

    #Flux per unit lethargy: divide by ln(E_high/E_low) so bins of different widths are comparable
    lethargy_width = np.log(E_high / E_low)
    flux_per_lethargy = flux_per_bin / lethargy_width

    #Air kerma delivered by each bin: flux x E x mu_en/rho(air), at the bin's geometric-mean energy
    E_mid_MeV = np.sqrt(E_low * E_high) / 1e6
    kerma_per_bin = flux_per_bin * E_mid_MeV * loglog_interp(E_mid_MeV, air_energy, air_mu_en_rho)
    #As a share of the total air kerma, per unit lethargy (so different thicknesses can be compared)
    kerma_share_per_lethargy = kerma_per_bin / kerma_per_bin.sum() / lethargy_width

    #Empty bins can't be drawn on a log axis
    flux_plot = np.where(flux_per_lethargy > 0, flux_per_lethargy, np.nan)
    kerma_plot = np.where(kerma_share_per_lethargy > 0, kerma_share_per_lethargy, np.nan)

    ax1.stairs(flux_plot, edges_keV, color=colour, linewidth=1.5, label=f'{t:g} cm lead')
    ax2.stairs(kerma_plot, edges_keV, color=colour, linewidth=1.5, label=f'{t:g} cm lead')

    #Total vs peak bin: should be close to the flux ratio and the buildup factor from the sweep.
    #Slightly lower, because the ~2% wide peak bin also catches very small-angle scatters.
    peak = np.searchsorted(E_high, E * 1e6)
    print(f'{t:10.1f} {flux_per_bin.sum() / flux_per_bin[peak]:22.3f} {kerma_per_bin.sum() / kerma_per_bin[peak]:30.3f}')

#Mark the physical features on both panels
for ax in (ax1, ax2):
    for energy, label in features:
        ax.axvline(energy, color=colour_feature, linestyle='--', linewidth=0.8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(50, 800)
    ax.set_xlabel('Photon energy (keV)')

#Label the features once, along the top of the left panel
for ax in (ax1, ax2):
    for energy, label in features:
        if label is not None:
            ax.text(energy, 1.02, label, transform=ax.get_xaxis_transform(), rotation=90,
                    ha='center', va='bottom', fontsize=8, color=colour_feature)

ax1.set_ylabel('Flux per unit lethargy (photons cm$^{-2}$ s$^{-1}$)')
ax1.set_title('Photon spectrum at the detector (1 m, 1 GBq Cs-137)', loc='left', pad=110)
ax2.set_ylabel('Share of total air kerma per unit lethargy')
ax2.set_title('Where the air kerma comes from', loc='left', pad=110)

#Keep gridlines and axes recessive so the data stands out
for ax in (ax1, ax2):
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

#One legend for both panels, below the plots so it doesn't cover the data
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(labels), frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1))

#Save the figure in the Results folder next to this script
results_dir = Path(__file__).resolve().parent / 'Results'
results_dir.mkdir(exist_ok=True)
fig.savefig(results_dir / 'spectrum.png', dpi=200)
print(f'\nPlot saved to {results_dir / "spectrum.png"}')

plt.show()