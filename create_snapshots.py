### Program to create quick snapshots for a given time for IBIS datasets

# setting up necessary modules
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from configparser import ConfigParser
import glob
import data_functions
import spectral_analysis
import os
import cmasher as cmr
import sunpy.visualization.colormaps as cm

plt.style.use("classic")

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font
# of your choice!)
# rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':10})
# These parameters can also be put into the style or matplotlibrc.
# This is the dynamic approach of changing parameters.
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
    "savefig.dpi": 300,
    "pcolor.shading": "auto",
}
matplotlib.rcParams.update(tex_fonts)

# Set the font used for MathJax - more on this later
rc("mathtext", **{"default": "regular"})

# The following %config line changes the inline figures to have a higher DPI.
# You can comment out (#) this line if you don't have a high-DPI
# (~220) display.
# plt.rcParams["figure.facecolor"] = "white"

############################### Directories ###################################3


# Contains information for the date

config_object = ConfigParser()
config_object.read("AGWs_config.ini")
date = "04262019"
written_date = "26Apr2019"
# Grab the location and plate scale information associated with chosen date
date_image_information = config_object[date]


data_dir = "/media/oana/Data1/AGWs/" + written_date + "/Data/"
path_to_save_figures = "/media/oana/Data1/AGWs/" + written_date + "/Figures/Snapshots/"

# If folder doesn't exist, create it
if not os.path.isdir(path_to_save_figures):
    os.makedirs(path_to_save_figures)

magnetogram_file = data_dir + written_date + ".HMImag.ibis.aligned.fits"
dopplergram_file = data_dir + written_date + ".HMIdop.ibis.aligned.fits"
continuum_file = data_dir + written_date + ".HMIcont.ibis.aligned.fits"
aia1600_file = data_dir + written_date + ".AIA1600.ibis.aligned.fits"
aia1700_file = data_dir + written_date + ".AIA1700.ibis.aligned.fits"


# determine which FITs file type to grab
fname_type_dict = {0: "int", 1: "vel"}
fname_type = fname_type_dict.get(1)


hmi_dop_file = data_dir + written_date + ".HMIdop.ibis.aligned.fits"

spectral_line = "fe7090"
vfe7090_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]


spectral_line = "k7699"
vk7699_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

spectral_line = "fe5434"
vfe5434_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

spectral_line = "ca8542"
vca8542_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

# determine which FITs file type to grab
fname_type_dict = {0: "int", 1: "vel"}
fname_type = fname_type_dict.get(0)


hmi_cont_file = data_dir + written_date + ".HMIcont.ibis.aligned.fits"

spectral_line = "fe7090"
ife7090_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]


spectral_line = "k7699"
ik7699_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

spectral_line = "fe5434"
ife5434_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

spectral_line = "ca8542"
ica8542_file = glob.glob(
    data_dir + written_date + "*" + fname_type + "*" + spectral_line + ".fits"
)[0]

############################### Read in Data #####################################

vfe7090_file = data_functions.read_file(vfe7090_file, 0, fits_file=True)
vfe7090_file -= np.mean(vfe7090_file)

vk7699_file = data_functions.read_file(vk7699_file, 0, fits_file=True)
vk7699_file -= np.mean(vk7699_file)

vfe5434_file = data_functions.read_file(vfe5434_file, 0, fits_file=True)
vfe5434_file -= np.mean(vfe5434_file)

vca8542_file = data_functions.read_file(vca8542_file, 0, fits_file=True)
vca8542_file -= np.median(vca8542_file[20, 80, :])


ife7090_file = data_functions.read_file(ife7090_file, 0, fits_file=True)
ife7090_file -= np.mean(ife7090_file)

ik7699_file = data_functions.read_file(ik7699_file, 0, fits_file=True)
ik7699_file -= np.mean(ik7699_file)

ife5434_file = data_functions.read_file(ife5434_file, 0, fits_file=True)
ife5434_file -= np.mean(ife5434_file)

ica8542_file = data_functions.read_file(ica8542_file, 0, fits_file=True)
ica8542_file -= np.mean(ica8542_file)


hmi_cont = data_functions.read_file(continuum_file, 0, fits_file=True)
hmi_cont -= np.mean(hmi_cont)

hmi_dop = data_functions.read_file(dopplergram_file, 0, fits_file=True)
hmi_dop -= np.mean(hmi_dop)

hmi_mag = data_functions.read_file(magnetogram_file, 0, fits_file=True)
hmi_mag -= np.mean(hmi_mag)

aia1600 = data_functions.read_file(aia1600_file, 0, fits_file=True)
aia1600 -= np.mean(aia1600)

aia1700 = data_functions.read_file(aia1700_file, 0, fits_file=True)
aia1700 -= np.mean(aia1700)


############################### Grid Information #################################

distance = 147.1e6  # Earth-Sun distance in km
conversion_arcseconds_to_Mm = spectral_analysis.conversion_arcseconds_to_Mm(distance)

# spatial scaling/sampling [arcseconds/pixel]
dx = float(date_image_information["dx"])

# cadence of data
dt = float(date_image_information["dt"])  # s

# Nyquist wavenumber [1/arcseconds]
kx = np.pi / dx
print(
    "Nyquist wavenumber is %.2f 1/arcsecs or %.2f 1/Mm or %.2f 1/km"
    % (
        kx,
        (kx / conversion_arcseconds_to_Mm),
        ((kx / conversion_arcseconds_to_Mm) / 1000),
    )
)

# Nyquist frequency/angular frequency [rad/s]
omega = np.pi / dt
# Nyquist frequency/cyclic frequency [Hz]
v = omega / (2 * np.pi)
frq = v * 1000  # [mHz]
print("Nyquist frequency is %s mHz" % frq)

# defining the meshgrid for the data to be displayed on
# The starting time
start = 0
# the duration of the data
end_time = vfe7090_file.shape[-1]
# the entire length of the data
end_space = vfe7090_file.shape[0]

mid_time = end_time // 2
# The real floor division operator is “//”. It returns floor value for both
# integer and floating point arguments.
mid_space = end_space // 2
print("The middle time is %s and the end time is %s " % (mid_time, end_time))
fov = dx * end_space  # arcseconds; field of view
print("FOV is %s x %s arcseconds" % (fov, fov))

# horizontal wavenumber [1/arcseconds]
horizontal_array = np.linspace(0.0, kx, int(mid_space), endpoint=True)  # 1/arcsec

# length of cube
arcsecond_length_array = np.linspace(0.0, fov, int(end_space), endpoint=True)
# creates meshgrid for the lengths axes
arcx, arcy = np.meshgrid(arcsecond_length_array, arcsecond_length_array)

print("The length meshgrid has the following shape: ", arcx.shape)
# frequency [mHz]
# freq_array = (np.fft.rfftfreq(end_time, d=dt)) * 1000
freq_array = np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz
print("The three first frequencies are ", freq_array[0:4])
print("The length of the frequency array is", len(freq_array))
# np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz

# creates meshgrid for the axes to be able to plot on imshow/pcolormesh
KH, NU = np.meshgrid(horizontal_array, freq_array)

print("The meshgrid has the following shape: ", KH.shape)

# compute grid to azimuthally average over
x = np.linspace(-mid_space, mid_space - 1, end_space)  # size of kx
y = np.linspace(-mid_space, mid_space - 1, end_space)  # sixe of ky
X, Y = np.meshgrid(x, y)
radial_dist = np.hypot(X, Y)
print(
    "The shape of the grid to compute the azimuthal averaging over is ",
    radial_dist.shape,
)

########################### Create Doppler Velocity Snapshots ###############################

cmap = "bwr"

color_regimes = "k"
color_fmode = "dimgray"
color_lamb = "k"
vminp = -1
vmaxp = 1

fig, axs = plt.subplots(
    1,
    4,
    figsize=[18, 5],
    sharex=True,
    sharey=True,
    constrained_layout=True,
)
fig.suptitle(written_date + " IBIS Velocity", fontsize=18)

im1 = axs[0].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    vfe7090_file[:, :, 10],
    cmap=cmap,
    shading="auto",
    vmin=vminp,
    vmax=vmaxp,
    rasterized=True,
)

axs[0].set_title("Fe\,{\Large I}\,7090")
axs[0].set_xlabel(r"[Mm]")
axs[0].set_ylabel(r"[Mm]")


im2 = axs[1].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    vk7699_file[:, :, 10],
    cmap=cmap,
    shading="auto",
    vmin=vminp,
    vmax=vmaxp,
    rasterized=True,
)
axs[1].set_xlabel(r"[Mm]")

axs[1].set_title(r"K\,{\Large I}\,7699")

im3 = axs[2].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    vfe5434_file[:, :, 10],
    cmap=cmap,
    shading="auto",
    vmin=vminp,
    vmax=vmaxp,
    rasterized=True,
)
axs[2].set_title(r"Fe\,{\Large I}\,5434")
axs[2].set_xlabel(r"[Mm]")


im4 = axs[3].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    vca8542_file[:, :, 10],
    cmap=cmap,
    shading="auto",
    vmin=vminp,
    vmax=vmaxp,
    rasterized=True,
)
axs[3].set_title("Ca\,{\Large II}\,8542")
axs[3].set_xlabel(r"[Mm]")

cbar = fig.colorbar(im4, ax=axs[3], aspect=40)

cbar.set_label(label=r"Doppler Velocity [km\,s$^{-1}$]", fontsize=15)


# label customization
for ax in axs.flat:
    ax.minorticks_on()
    ax.xaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_major_locator(mticker.AutoLocator())

    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())

    ax.tick_params(axis="both", which="major", width=1)
    ax.tick_params(axis="both", which="minor", width=1)

    ax.set_xlim(0, arcx.max() * conversion_arcseconds_to_Mm)
    ax.set_ylim(0, arcy.max() * conversion_arcseconds_to_Mm)


fig.set_constrained_layout_pads(wspace=-0.0001, hspace=-0.0001)


fig.savefig(
    path_to_save_figures + written_date + "_IBIS_lc_vel_snapshots.pdf",
    dpi=400,
    bbox_inches="tight",
)
plt.show()

########################### Create Intensity Snapshots ###############################

cmap = cmr.neutral

color_regimes = "k"
color_fmode = "dimgray"
color_lamb = "k"

fig, axs = plt.subplots(
    1,
    4,
    figsize=[18, 5],
    sharex=True,
    sharey=True,
    constrained_layout=True,
)
fig.suptitle(written_date + " IBIS Intensity", fontsize=18)


im1 = axs[0].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(ife7090_file[:, :, 10], 0, 1),
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
axs[0].set_title("Fe\,{\Large I}\,7090")
axs[0].set_xlabel(r"[Mm]")
axs[0].set_ylabel(r"[Mm]")


im2 = axs[1].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(ik7699_file[:, :, 10], 0, 1),
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
axs[1].set_xlabel(r"[Mm]")
axs[1].set_title("K\,{\Large I}\,7699")

im3 = axs[2].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(ife5434_file[:, :, 10], 0, 1),
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
axs[2].set_xlabel(r"[Mm]")
axs[2].set_title("Fe\,{\Large I}\,5434")


im4 = axs[3].pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(ica8542_file[:, :, 10], 0, 1),
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
axs[3].set_title("Ca\,{\Large II}\,8542")
axs[3].set_xlabel(r"[Mm]")

cbar = fig.colorbar(im4, ax=axs[3], aspect=40)

cbar.set_label(label=r"Intensity [counts]", fontsize=15)


# label customization
for ax in axs.flat:
    ax.xaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_major_locator(mticker.AutoLocator())

    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())

    ax.tick_params(axis="both", which="major", width=1)
    ax.tick_params(axis="both", which="minor", width=1)

    ax.set_xlim(0, arcx.max() * conversion_arcseconds_to_Mm)
    ax.set_ylim(0, arcy.max() * conversion_arcseconds_to_Mm)


fig.set_constrained_layout_pads(wspace=-0.0001, hspace=-0.0001)


fig.savefig(
    path_to_save_figures + written_date + "_IBIS_lc_int_snapshots.pdf",
    dpi=400,
    bbox_inches="tight",
)
plt.show()


########################### Create HMI-AIA Snapshots ###############################

############################### Grid Information #################################

distance = 147.1e6  # Earth-Sun distance in km
conversion_arcseconds_to_Mm = spectral_analysis.conversion_arcseconds_to_Mm(distance)

# spatial scaling/sampling [arcseconds/pixel]
dx = float(date_image_information["dx"])

# cadence of data
dt = float(date_image_information["dt"])  # s

# Nyquist wavenumber [1/arcseconds]
kx = np.pi / dx
print("Nyquist wavnumber is %s 1/arcsecs" % kx)
print("Nyquist wavnumber is %s 1/Mm" % (kx / conversion_arcseconds_to_Mm))
print("Nyquist wavnumber is %s 1/km" % ((kx / conversion_arcseconds_to_Mm) / 1000))

# Nyquist frequency/angular frequency [rad/s]
omega = np.pi / dt
# Nyquist frequency/cyclic frequency [Hz]
v = omega / (2 * np.pi)
frq = v * 1000  # [mHz]
print("Nyquist frequency is %s mHz" % frq)

# defining the meshgrid for the data to be displayed on
# The starting time
start = 0
# the duration of the data
end_time = hmi_cont.shape[-1]
# the entire length of the data
end_space = hmi_cont.shape[0]

mid_time = end_time // 2
# The real floor division operator is “//”. It returns floor value for both
# integer and floating point arguments.
mid_space = end_space // 2
print("The middle time is %s and the end time is %s " % (mid_time, end_time))
fov = dx * end_space  # arcseconds; field of view
print("FOV is %s arcseconds" % (fov))

# horizontal wavenumber [1/arcseconds]
horizontal_array = np.linspace(0.0, kx, int(mid_space), endpoint=True)  # 1/arcsec

# length of cube
arcsecond_length_array = np.linspace(0.0, fov, int(end_space), endpoint=True)
# creates meshgrid for the lengths axes
arcx, arcy = np.meshgrid(arcsecond_length_array, arcsecond_length_array)

print("The length meshgrid has the following shape: ", arcx.shape)
# frequency [mHz]
# freq_array = (np.fft.rfftfreq(end_time, d=dt)) * 1000
freq_array = np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz
print("The three first frequencies are ", freq_array[0:4])
print("The length of the frequency array is", len(freq_array))
# np.linspace(0.0, frq, int(mid_time), endpoint=True)  # mHz

# creates meshgrid for the axes to be able to plot on imshow/pcolormesh
KH, NU = np.meshgrid(horizontal_array, freq_array)

print("The meshgrid has the following shape: ", KH.shape)

# compute grid to azimuthally average over
x = np.linspace(-mid_space, mid_space - 1, end_space)  # size of kx
y = np.linspace(-mid_space, mid_space - 1, end_space)  # sixe of ky
X, Y = np.meshgrid(x, y)
radial_dist = np.hypot(X, Y)
print(
    "The shape of the grid to compute the azimuthal averaging over is ",
    radial_dist.shape,
)

############################### Plotting #################################


fig = plt.figure(figsize=[16, 8], facecolor="white")
fig.suptitle(written_date + " SDO", fontsize=18)

ax1 = fig.add_subplot(2, 3, 1)
im1 = ax1.pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(hmi_cont[:, :, 10], 0, 1),
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
ax1.set_title("HMI Continuum")
cbar = fig.colorbar(im1, ax=ax1, aspect=40)
cbar.set_label(label=r"Intensity", fontsize=15)


ax2 = fig.add_subplot(2, 3, 2)
im2 = ax2.pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    hmi_mag[:, :, 10],
    cmap=cmap,
    shading="auto",
    rasterized=True,
)
ax2.set_title("HMI Magnetogram")
cbar = fig.colorbar(im2, ax=ax2, aspect=40)
cbar.set_label(label=r"LOS B [G]", fontsize=15)

ax3 = fig.add_subplot(2, 3, 3)
im3 = ax3.pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    hmi_dop[:, :, 10] * 1e-4,
    cmap="bwr",
    vmin=-1,
    vmax=1,
    shading="auto",
    rasterized=True,
)
ax3.set_title("HMI Dopplergram")
cbar = fig.colorbar(im3, ax=ax3, aspect=40)
cbar.set_label(label=r"Velocity [km/s]", fontsize=15)

ax4 = fig.add_subplot(2, 3, 4)
im4 = ax4.pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(aia1600[:, :, 10], 0, 1),
    cmap=matplotlib.colormaps["sdoaia1600"],
    shading="auto",
    rasterized=True,
)
ax4.set_title("AIA 1600")
cbar = fig.colorbar(im4, ax=ax4, aspect=40)
cbar.set_label(label=r"Intensity", fontsize=15)

ax5 = fig.add_subplot(2, 3, 5)
im5 = ax5.pcolormesh(
    arcx * conversion_arcseconds_to_Mm,
    arcy * conversion_arcseconds_to_Mm,
    data_functions.data_normalization(aia1700[:, :, 10], 0, 1),
    cmap=matplotlib.colormaps["sdoaia1700"],
    shading="auto",
    rasterized=True,
)
ax5.set_title("AIA 1700")
cbar = fig.colorbar(im5, ax=ax5, aspect=40)
cbar.set_label(label=r"Intensity", fontsize=15)

axes_list = [ax1, ax2, ax3, ax4, ax5]
nom_subplots = len(axes_list)

for ax in axes_list:
    ax.set_xlabel(r"[Mm]")
    ax.set_ylabel(r"[Mm]")
    ax.set_xlim(0, arcx.max() * conversion_arcseconds_to_Mm)
    ax.set_ylim(0, arcy.max() * conversion_arcseconds_to_Mm)
    ax.grid(False)

fig.tight_layout()

fig.savefig(
    path_to_save_figures + written_date + "SDO_snapshots.pdf",
    dpi=400,
    bbox_inches="tight",
)
plt.show()
