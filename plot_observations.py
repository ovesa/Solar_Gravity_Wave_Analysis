## Routine to plot observations on a blank Sun

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from matplotlib import rc
from sunpy.coordinates import frames

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
# rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':10})
# These parameters can also be put into the style or matplotlibrc.
# This is the dynamic approach of changing parameters.
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Times New Roman",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 13,
    "font.size": 13,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 12,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.titlesize": 14,
}
matplotlib.rcParams.update(tex_fonts)

# Set the font used for MathJax - more on this later
rc("mathtext", **{"default": "regular"})

# The following %config line changes the inline figures to have a higher DPI.
# You can comment out (#) this line if you don't have a high-DPI (~220) display.
plt.style.use("classic")
plt.rcParams["figure.facecolor"] = "white"

############################### Define Fake Sun ###################################

data = np.full((10, 10), np.nan)

# Define a reference coordinate and create a header using sunpy.map.make_fitswcs_header
skycoord = SkyCoord(
    0 * u.arcsec,
    0 * u.arcsec,
    obstime="2019-04-25",
    observer="earth",
    frame=frames.Helioprojective,
)

# Scale set to the following for solar limb to be in the field of view
header = sunpy.map.make_fitswcs_header(
    data, skycoord, scale=[220, 220] * u.arcsec / u.pixel
)

# Use sunpy.map.Map to create the blank map
blank_map = sunpy.map.Map(data, header)

############################### Define Observation Coordinates ###################################


def plot_observation_coordinates(lon, lat, date_string, time_string, nframe):
    obs_date = date_string + "T" + time_string
    coord = SkyCoord(
        lon * u.deg,
        lat * u.deg,
        obstime=obs_date,
        observer="earth",
        frame=frames.HeliographicStonyhurst,
    )
    coordvb = coord.transform_to(frame=frames.Helioprojective)

    coord = SkyCoord(
        coordvb.Tx.value * u.arcsec,
        coordvb.Ty.value * u.arcsec,
        frame=nframe.coordinate_frame,
    )
    return coord


coord_25Apr2019 = plot_observation_coordinates(
    6.8315904, -3.27290550, "2019-04-25", "14:15:17.542", blank_map
)
coord_26Apr2019 = plot_observation_coordinates(
    5.684587800000e01, 8.886901100000e-01, "2019-04-26", "15:24:07.049", blank_map
)
coord_06May2019 = plot_observation_coordinates(
    -4.888134700000e01, 8.590126200000e00, "2019-05-06", "14:34:40.85", blank_map
)
coord_09May2019 = plot_observation_coordinates(
    -1.048208400000e01, 9.395851600000e00, "2019-05-09", "15:37:08.584", blank_map
)

coord_08May2019 = plot_observation_coordinates(
    -2.464746500000e01, 1.085355400000e01, "2019-05-08", "13:58:30.96", blank_map
)

############################### Plot Figure ###################################

fig = plt.figure()
ax = plt.subplot(projection=blank_map)

blank_map.plot(axes=ax)

blank_map.draw_limb(axes=ax, color="k")
blank_map.draw_grid(axes=ax, color="k")

ax.plot_coord(
    coord_25Apr2019,
    color="k",
    marker="o",
    linewidth=10,
    markersize=12,
    linestyle="None",
    label="25Apr2019",
)
ax.plot_coord(
    coord_26Apr2019,
    color="r",
    marker="o",
    linewidth=10,
    markersize=12,
    linestyle="None",
    label="26Apr2019",
)
ax.plot_coord(
    coord_06May2019,
    color="blue",
    marker="o",
    linewidth=10,
    markersize=12,
    linestyle="None",
    label="06May2019",
)
ax.plot_coord(
    coord_08May2019,
    color="magenta",
    marker="o",
    linewidth=10,
    markersize=12,
    linestyle="None",
    label="08May2019",
)
ax.plot_coord(
    coord_09May2019,
    color="orange",
    marker="o",
    linewidth=10,
    markersize=12,
    linestyle="None",
    label="09May2019",
)

ax.set_title("2019 AGW Observations")
ax.legend(
    numpoints=1,
    ncol=1,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fancybox=True,
    shadow=True,
)
plt.show()
