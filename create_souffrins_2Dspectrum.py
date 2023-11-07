# setting up necessary modules

import cmasher as cmr
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
from matplotlib import rc
from tqdm import tqdm
from astropy.io import ascii
import souffrins_dispersion_equations
import spectral_analysis
import isothermal_dispersion_equations

tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Times New Roman",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 14,
    "font.size": 14,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": True,
    "legend.framealpha": 0.8,
    "axes.formatter.use_mathtext": True,
    # "figure.dpi": 300,
    # "figure.autolayout": True,
    "lines.linewidth": 2,
    "ytick.minor.visible": True,
    "xtick.minor.visible": True,
    "figure.facecolor": "white",
    # "savefig.bbox": "tight",
    # "savefig.dpi": 300,
    "pcolor.shading": "auto",
}
matplotlib.rcParams.update(tex_fonts)

# Set the font used for MathJax - more on this later
rc("mathtext", **{"default": "regular"})

# The following %config line changes the inline figures to have a higher DPI.
# You can comment out (#) this line if you don't have a high-DPI (~220) display.
plt.style.use("classic")
plt.rcParams["figure.facecolor"] = "white"

###########################################################################################333

distance = 147.1e6  # Earth-Sun distance in km

conversion_arcseconds_to_Mm = spectral_analysis.conversion_arcseconds_to_Mm(distance)

# spatial scaling/sampling [arcseconds/pixel]
dx = 0.6

# cadence of data
dt = 11.88  # s

# Nyquist wavenumber [1/arcseconds]
kx = np.pi / dx
print("Nyquist wavnumber is %.3f 1/arcseconds" % kx)
print("Nyquist wavnumber is %.3f 1/Mm" % (kx / conversion_arcseconds_to_Mm))
print("Nyquist wavnumber is %.3e 1/km" % (kx / conversion_arcseconds_to_Mm / 1000))

# Nyquist frequency/angular frequency [rad/s]
omega = np.pi / dt
# Nyquist frequency/cyclic frequency [Hz]
v = omega / (2 * np.pi)
frq = v * 1000  # [mHz]
print("Nyquist frequency is %.2f mHz" % frq)

# defining the meshgrid for the data to be displayed on
# The starting time
start = 0
# the duration of the data
end_time = 840
# the entire length of the data
end_space = 160
mid_time = end_time // 2
# The real floor division operator is “//”. It returns floor value for both
# integer and floating point arguments.
mid_space = end_space // 2
print("The middle time is %.1f s and the end time is %.1f s" % (mid_time, end_time))

# horizontal wavenumber [1/arcseconds]
horizontal_array = np.linspace(0.0, kx, int(mid_space), endpoint=True)  # 1/Mm
kx_km = horizontal_array / conversion_arcseconds_to_Mm / 1000  # 1/km

# frequency [mHz]
# freq_array = (np.fft.rfftfreq(end_time, d=dt)) * 1000
freq_array2 = np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz
print("The three first frequencies are ", freq_array2[0:4])
print("The length of the frequency array is", len(freq_array2))
# np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz
omega_array = np.linspace(-omega, omega, 840)  # (rad*hz)


# creates meshgrid for the axes to be able to plot on imshow/pcolormesh
KH, NU = np.meshgrid(horizontal_array, freq_array2)

print("The meshgrid has the following shape: ", KH.shape)


####################### best "known" guesses #######################

wac_val = 5.4  # mHz  ~photospheric value
tau_val = 200  # s   ~ Mihalas & Toomre, 1982 Fig. 3
cs_val = 7.0  # km/s  ~ average photospheric value
z_val = 150  # km  ~ 7699: 470 km / 5434: 570 km. del z ~ 100 km
# convert mHz to angular frequency
mHz_to_Hz_factor = (2 * np.pi) / 1000

wavedisp2d = np.empty([420, len(kx_km)])


for kx in tqdm(range(0, len(kx_km))):
    kx_val = kx_km[
        kx
    ]  # 1/km.  # should be a relatively small number (1e-4 to 1e-1) Miihalas & Toomre. 1982 Fig. 4

    model_vals = np.array([wac_val * mHz_to_Hz_factor, cs_val, tau_val, z_val, kx_val])
    model = souffrins_dispersion_equations.merged_model_fmode(
        omega_array[420:], *model_vals
    )
    wavedisp2d[:, kx] = model

####################### Dispersion Relations #######################

# loads in model atmosphere of the Sun in order to calculate the
# dispersion relation

# col1 = height (Mm); col2 = gravity (cm/s^2); col3 = dP/dz; col4 = density(g/cm^3); col5 = dRho_dz; col6 = sound speed (cm/s)
data = ascii.read("/media/oana/Data1/AGWs/IBISdata/CSM_A.dat")
cs_cgs = data["col6"]  # cm/s
grav = data["col2"]  # cm/s²
height = data["col1"]  # Mm
density_cgs = data["col4"]  # g/cm^3
dRho_dz = data["col5"]  # ?? should be in cgs...
dP_dz = data["col3"]  # ?? should be in cgs...
density_Mm = density_cgs / (1e-8) ** 3  # g/Mm^3
gravity_cgs = 27400  # cm/s²
gravity_Mm = 274 / 1e6  # Mm/s²
gravity_km = 0.274  # km/s²
cs_Mm = cs_cgs * 1e-8  # Mm/s
cs_km = cs_cgs * 1e-5  # km/s


height_index = 1005  # # for a height of 0.25 Mm or 250 km

# fundamental mode (f-mode) curve; w = sqrt(k_h*g)
surface_fmode = isothermal_dispersion_equations.surface_fmodes(
    gravity_km, kx_km
)  # angular frequency [Hz]

# Lamb frequency; w = c_s*k_h
lamb_curve = isothermal_dispersion_equations.lamb_frequency(
    cs_km[height_index], kx_km
)  # angular frequency [Hz]

# Density Scale Height (cm)
density_scale_height_value = isothermal_dispersion_equations.densityscaleheight(
    dRho_dz[height_index], density_cgs[height_index]
)

# acoustic cut-off frequency squared
acoustic_cutoff_frequency_squared = (
    isothermal_dispersion_equations.acoustic_cutoff_frequency(
        cs_cgs[height_index], density_scale_height_value
    )
)  # angular frequency squared  [Hz^2]


Brunt_Vaisala_frequency_Squared_value = (
    isothermal_dispersion_equations.Brunt_Vaisala_frequency_Squared(
        gravity_cgs,
        dP_dz[height_index],
        cs_cgs[height_index],
        density_cgs[height_index],
        dRho_dz[height_index],
    )
)  # angular frequency squared [Hz^2]


# solving the isothermal dispersion relation in the absence of a
# magnetic field

a = 1
b = -(acoustic_cutoff_frequency_squared + cs_km[height_index] ** 2 * kx_km**2)
c = Brunt_Vaisala_frequency_Squared_value * cs_km[height_index] ** 2 * kx_km**2


qval1, qval2 = isothermal_dispersion_equations.quadraticforumala(a, b, c)
# badq1, badq2 = quadraticforumala(a,b,c)

####################### plot #######################

spectrum = np.rad2deg(wavedisp2d)
cmap = cmr.fusion

color_fmode = "gray"
color_lamb = "gray"
color_regime = "k"

style.use("classic")
fig = plt.figure(facecolor="white", figsize=[6, 5])

ax1 = plt.subplot(111)
im1 = ax1.pcolormesh(
    KH / conversion_arcseconds_to_Mm,
    NU,
    spectrum[:, :],
    cmap=cmap,
    shading="auto",
    norm=mcolors.TwoSlopeNorm(vmin=-80.0, vcenter=0.0, vmax=80),
    rasterized=True,
)

# ticks=np.arange(-20, 81, 20)
cbar = plt.colorbar(im1, ax=ax1, pad=0.08)  # , format=ticker.FuncFormatter(fmt_cmpers))
cbar.ax.set_yscale("linear")
cbar.set_label("Phase Difference [deg]", labelpad=10, fontsize=13)
ax1.set_ylabel(r"$\nu \ [{\rm mHz}]$", fontsize=13)
ax1.set_xlabel(r"$k_{\rm h} \ [{\rm Mm^{-1}}]$", fontsize=13)
ax1.set_xlim(0, horizontal_array.max() / conversion_arcseconds_to_Mm)
ax1.set_ylim(0, 10)
ax1.tick_params(which="both")
ax1.minorticks_on()
ax1.tick_params(axis="both", which="major", width=1)
ax1.tick_params(axis="both", which="minor", width=1)


ax1.plot(
    horizontal_array / conversion_arcseconds_to_Mm,
    np.sqrt(qval1) * 1000 / (2 * np.pi),
    color=color_regime,
    linewidth=2.5,
)
ax1.plot(
    horizontal_array / conversion_arcseconds_to_Mm,
    np.sqrt(qval2) * 1000 / (2 * np.pi),
    color=color_regime,
    linewidth=2.5,
)
ax1.plot(
    horizontal_array / conversion_arcseconds_to_Mm,
    surface_fmode * 1000 / (2 * np.pi),
    color=color_fmode,
    linewidth=2.5,
)
ax1.plot(
    horizontal_array / conversion_arcseconds_to_Mm,
    lamb_curve * 1000 / (2 * np.pi),
    color=color_lamb,
    linewidth=2.5,
    linestyle="--",
)

# ####################### Label Plot #######################

# t = ax1.text(
#     4.5,
#     2,
#     "A",
#     ha="center",
#     va="center",
#     size=15,
#     bbox=dict(boxstyle="circle", fc="w", ec="k", lw=2),
# )

# t2 = ax1.text(
#     2.5,
#     8,
#     "C",
#     va="center",
#     ha="center",
#     size=15,
#     bbox=dict(boxstyle="circle,pad=0.3", fc="w", ec="k", lw=2),
# )


# # t2 = plt.annotate(
# #     "C",
# #     xy=(2.3, 7),
# #     xytext=(2.5, 7),
# #     size=15,
# #     bbox=dict(boxstyle="circle,pad=0.3", fc="w", ec="k", lw=2),
# # )


# t3 = ax1.text(
#     0.5,
#     1.4,
#     "B",
#     ha="center",
#     va="center",
#     size=15,
#     bbox=dict(boxstyle="circle", fc="w", ec="k", lw=2),
# )

# t4 = ax1.text(
#     2.5,
#     4.2,
#     "F",
#     ha="center",
#     va="center",
#     size=15,
#     bbox=dict(boxstyle="circle,pad=0.3", fc="w", ec="k", lw=2),
#     rotation_mode="anchor",
#     transform_rotates_text=True,
# )

# t5 = ax1.text(
#     0.8,
#     4,
#     "E",
#     ha="center",
#     va="center",
#     size=15,
#     bbox=dict(boxstyle="circle,pad=0.3", fc="w", ec="k", lw=2),
# )

# # slope y=mx+b

# ax1.plot([0, 1.6], [2, 2], c="darkslategrey", linewidth=2.2)
# ax1.plot([0, 1.6], [0, 2], c="darkslategrey", linewidth=2.2)


# ax1.set_title(r"$k_{\rm h} - \nu$ Diagnostic Diagram")
# # text
# ax1.text(
#     4.8,
#     1.5,
#     "Atmospheric Gravity Waves",
#     fontsize=12,
#     ha="center",
#     bbox=dict(facecolor="white", alpha=1, pad=6),
# )
# # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},


# ax1.text(
#     2.3,
#     7.8,
#     "Acoustic Waves",
#     ha="center",
#     fontsize=12,
#     bbox=dict(facecolor="white", alpha=1, pad=6),
# )

# ax1.text(
#     7.4,
#     5.0,
#     r"$\frac{N}{2 \pi}$",
#     fontsize=16,
#     # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},
# )

# ax1.text(
#     -0.55,
#     5.2,
#     r"$\nu_{\rm ac}$",
#     fontsize=16
#     # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},
# )

# ax1.text(
#     0.7,
#     2.7,
#     r"f-mode",
#     fontsize=15,
#     rotation=38, )

#     # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},


# ax1.text(
#     1.8,
#     4.0,
#     "Evanescent Waves",
#     ha="center",
#     fontsize=12,
#     bbox=dict(facecolor="white", alpha=1, pad=4),
# )
# #    bbox=dict(facecolor="white", alpha=0.9, pad=0.4, boxstyle="round"),

# # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},


# ax1.text(
#     3.4,
#     4.0,
#     r"Lamb Waves",
#     fontsize=15,
#     rotation=35
#     # bbox={"facecolor": "white", "alpha": 0.5, "pad": 10},
# )

plt.tight_layout()
# fig.savefig(
#     "Labelled_Souffrin_Diagnostic_Diagram.pdf",
#     dpi=400,
#     bbox_inches="tight",
# )
plt.show()
