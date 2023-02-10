import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap, LogNorm
import plotly.express as px
from astropy.modeling import fitting, polynomial
from scipy.optimize import minimize
import plotly.graph_objects as go


# --------------------------------------------- Plot Color Specifications ---------------------------------------------
uc_colors = {'red': '#C61A27', 'blue': '#3891A6', 'orange': '#F79D65', 'yellow': '#FDE74C', 'lightblue': '#C1FFF2',
             'green': '#9BC53D', 'white': '#FFFFFF', 'black': '#000000', 'gray': '#999999'}
top = cm.get_cmap('Reds', 128 * 4)
bottom = cm.get_cmap('BrBG_r', 128 * 4)
top_r = cm.get_cmap('Reds_r', 128 * 4)
bottom_r = cm.get_cmap('BrBG', 128 * 4)
newcolors = np.vstack((bottom(np.linspace(0, 0.5, 128 * 4)), top(np.linspace(0, 1, 128 * 4))))
newcolors_r = np.vstack((top_r(np.linspace(0, 1, 128 * 4)), bottom_r(np.linspace(0.5, 1, 128 * 4))))
newcmp = ListedColormap(newcolors, name='RedBlueish')

# ----------------------------------------------- Physical Parameters -----------------------------------------------
wavelength = 671e-9  # [m]
numerical_aperture = 0.656  # Special Optics: 0.656, Mitutoyo: 0.3
pixel_size = 2.4e-6  # [m], Blackfly: 2.4µm Chameleon: 3.7µm
camera_x_pixels, camera_y_pixels = 5472, 3648
focus_error = 1e-6  # [m] Error of the nanohole position along the optical axis
focal_length_obj = 18.8e-3  # [m] Focal length of objective
focal_length_lens = 300e-3  # [m] The lens determining the magnification used to focus light on the camera
effective_focal_length = 1 / (1 / focal_length_obj + 1 / focal_length_lens)
magnification = focal_length_lens / focal_length_obj
magnification_error = effective_focal_length / (effective_focal_length - focal_length_obj) ** 2 * focus_error
resolution_theo = 0.61 * wavelength / numerical_aperture
tweezer_position = np.array([2757, 1847]) * pixel_size / magnification

# ------------------------------------------------- Data Directory -------------------------------------------------
# Folders required in the working directory: "Images", "Results"
working_directory = "files/"

# ---------------------------------------------- Reading Results Files -----------------------------------------------
data = {}

for filename in sorted(os.listdir(working_directory + "Results")):
    if filename.endswith("txt"):
        if filename.startswith("psf_filenames"):
            data[filename[:-4]] = np.loadtxt(working_directory + "Results/" + filename, dtype='U')
        else:
            data[filename[:-4]] = np.loadtxt(working_directory + "Results/" + filename)

# Coordinates: image plane [px] to object plane [µm] conversion:
data["all_coordinates"] = data["all_coordinates"] * pixel_size / magnification

# ---------------------------------------------- Find Best Resolution  ---------------------------------------------
rmin, rsecondmin = np.partition(data["all_resolutions"], 1)[0:2]

# Use either rmin or rsecondmin, whichever makes sense
i, = np.where(np.isclose(data["all_resolutions"], rsecondmin, atol=1e-10))
print("Best Image", data["psf_filenames"][i])
print("Best Res. ", data["all_resolutions"][i])
print("Best Res. Err. ", data["all_resolutions_error"][i])
print("x" + str(data["all_coordinates"][0][i] * 1e6) + "µm")
print("y" + str(data["all_coordinates"][1][i] * 1e6) + "µm")

tweezer_position = np.array([2757, 1847]) * pixel_size / magnification  # Pixel coordinates given by add. measurement...
tweezer_error = 0.5e-3 / 35e-2 * focal_length_obj  # Beam overlap accuary / distance from objective * focallength => [m]

# --------------------------------------------------- 2D Fitting ----------------------------------------------------
fit_func = fitting.LinearLSQFitter()

model = polynomial.Polynomial2D(3)

fit = fit_func(model, data["all_coordinates"][0], data["all_coordinates"][1], data["all_resolutions"])

x = np.linspace(0, camera_x_pixels * pixel_size / magnification, 500)
y = np.linspace(0, camera_y_pixels * pixel_size / magnification, 500)
x_i, y_i = np.meshgrid(x, y)


def polynom2d_third_deg(x, c):
    return c["c0_0"] + c["c1_0"] * x[0] + c["c2_0"] * x[0] ** 2 + c["c3_0"] * x[0] ** 3 \
           + c["c0_1"] * x[1] + c["c0_2"] * x[1] ** 2 + c["c0_3"] * x[1] ** 3 \
           + c["c1_1"] * x[0] * x[1] + c["c1_2"] * x[0] * x[1] ** 2 \
           + c["c2_1"] * x[0] ** 2 * x[1]


def polynom2d_fourth_deg(x, c):
    return c["c0_0"] + c["c1_0"] * x[0] + c["c2_0"] * x[0] ** 2 + c["c3_0"] * x[0] ** 3 + c["c4_0"] * x[0] ** 4 \
           + c["c0_1"] * x[1] + c["c0_2"] * x[1] ** 2 + c["c0_3"] * x[1] ** 3 + c["c0_4"] * x[1] ** 4 \
           + c["c1_1"] * x[0] * x[1] + c["c2_2"] * x[0] ** 2 * x[1] ** 2 \
           + c["c1_2"] * x[0] * x[1] ** 2 + c["c1_3"] * x[0] * x[1] ** 3 \
           + c["c2_1"] * x[0] ** 2 * x[1] + c["c3_1"] * x[0] ** 3 * x[1]


c_values = dict(zip(fit.param_names, fit.parameters))


def fit_poly(x):
    # return polynom2d_fourth_deg(x, c_values)  # Use this for the fourth degree polynom
    return polynom2d_third_deg(x, c_values)  # Use this for the third degree polynom


# Minimize
minimum = minimize(fit_poly, x0=[0.0004, 0.0003])

# ---------------------------- Calculate the distance between the optical and mechanical axes -----------------------
# Distance between the axes intersections in the focal plane
distance = np.sqrt((tweezer_position[0] - minimum.x[0]) ** 2 + (tweezer_position[1] - minimum.x[1]) ** 2)
distance_error = abs((tweezer_position[0] - minimum.x[0]) / distance * tweezer_error) \
                 + abs((tweezer_position[1] - minimum.x[1]) / distance * tweezer_error)

# Conversion to the axes' angle difference
distance_to_angle = np.arctan((distance / focal_length_obj))
distance_to_angle_error = focal_length_obj / (
        distance ** 2 + focal_length_obj ** 2) * distance_error  # Deriv.: d/ddist arctan(dist/f)

# Test Plot:
# plt.imshow(polynom2d_third_deg([x_i, y_i],c_values),
#            extent=[0, camera_x_pixels*pixel_size/magnification, camera_y_pixels*pixel_size/magnification, 0])
# plt.scatter(minimum.x[0], minimum.x[1])
# plt.show()


# %% ------------------------------------- Calc. Mean Resolution around Minimum ----------------------------------------
fov_radius = 100e-6  # [m]
i, = np.where((data["all_coordinates"][0] < minimum.x + fov_radius) &
              (data["all_coordinates"][0] > minimum.x - fov_radius))
j, = np.where((data["all_coordinates"][1] < minimum.y + fov_radius) &
              (data["all_coordinates"][1] > minimum.y - fov_radius))
k = [x for x in i if x in j]
resolutions_in_FOV = data["all_resolutions"][k]
resolutions_in_FOV_mean = np.mean(resolutions_in_FOV)

mean_error = np.sqrt((np.std(resolutions_in_FOV) / np.sqrt(len(resolutions_in_FOV))) ** 2
                     + np.mean(data["all_resolutions_error"][k]) ** 2 / len(resolutions_in_FOV))
print("FOV radius", fov_radius)
print("r_mean", resolutions_in_FOV_mean)
print("error", mean_error)

# %% ---------------------------------------- PLOT FITTED FOV (INTERACTIVE) ----------------------------------------
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 30})

fig, main_ax = plt.subplots(figsize=(27, 15))
divider = make_axes_locatable(main_ax)
top_ax = divider.append_axes("top", 3, pad=0.3, sharex=main_ax)
right_ax = divider.append_axes("right", 3, pad=0.3, sharey=main_ax)

# make some labels invisible
top_ax.xaxis.set_tick_params(labelbottom=False)
right_ax.yaxis.set_tick_params(labelleft=False)

# Format Labels to µm
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x * 1e6))
ticks_cb = ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x * 1e6))
main_ax.xaxis.set_major_formatter(ticks)
main_ax.yaxis.set_major_formatter(ticks)
top_ax.yaxis.set_major_formatter(ticks)
right_ax.xaxis.set_major_formatter(ticks)

main_ax.set_xlabel('x Atomplane [µm]', labelpad=15)
main_ax.set_ylabel('y Atomplane [µm]', labelpad=15)
top_ax.set_ylabel('Resolution [µm]', labelpad=15)
right_ax.set_xlabel('Resolution [µm]', labelpad=15)

coord_conversion_x = 1 / (magnification / pixel_size / camera_x_pixels * 500)
coord_conversion_y = 1 / (magnification / pixel_size / camera_y_pixels * 500)

x, y = x_i, y_i
pos = np.empty(x.shape + (2,))
pos[:, :, 0] = x
pos[:, :, 1] = y
z = polynom2d_third_deg([x_i, y_i], c_values)
z_max = polynom2d_third_deg([x_i, y_i], c_values).max()

# Fit Plot
im = main_ax.imshow(polynom2d_third_deg([x_i, y_i], c_values),
                    extent=[0, camera_x_pixels * pixel_size / magnification,
                            camera_y_pixels * pixel_size / magnification, 0],
                    cmap="RdBu_r",
                    # vmin=0.775e-6,
                    # vmax=2.85e-6,
                    alpha=0.9,
                    norm=LogNorm(0.7e-6, 2.85e-6)
                    )
# Contour Plot
cont = main_ax.contour(z, levels=np.arange(rsecondmin * 1.1, 5e-6, 0.2e-6), origin="upper", cmap="Greys_r",
                       linewidths=2,
                       alpha=0.6,
                       extent=[0, camera_x_pixels * pixel_size / magnification,
                               camera_y_pixels * pixel_size / magnification, 0])

main_ax.autoscale(enable=False)
right_ax.autoscale(enable=False)
right_ax.grid()
top_ax.autoscale(enable=False)
top_ax.grid()
right_ax.set_xlim(right=z_max)
top_ax.set_ylim(top=z_max)
v_line = main_ax.axvline(np.nan, color=uc_colors["black"], lw=2, ls="--")  # 3891A6
h_line = main_ax.axhline(np.nan, color=uc_colors["black"], lw=2, ls="--")  # F79D65
v_prof, = right_ax.plot(np.zeros(x.shape[1]), np.arange(x.shape[1]) * coord_conversion_y, uc_colors["black"], lw=4)
h_prof, = top_ax.plot(np.arange(x.shape[0]) * coord_conversion_x, np.zeros(x.shape[0]), uc_colors["black"], lw=4)

top_ax.set_ylim([0, 3.5e-6])
top_ax.set_yticks([0, 1e-6, 2e-6, 3e-6])
yticks = top_ax.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)

right_ax.set_xlim([0, 3.5e-6])
right_ax.set_xticks([0, 1e-6, 2e-6, 3e-6])


def on_move(event):
    if event.inaxes is main_ax:
        # UNCOMMENT FOR FIXED POSITION
        # cur_x = minimum.x[0] * coord_conversion_x ** -1
        # cur_y = minimum.x[1] * coord_conversion_y ** -1

        # UNCOMMENT FOR INTERACTIVE POSITIONING:
        cur_x = event.xdata * coord_conversion_x ** -1
        cur_y = event.ydata * coord_conversion_y ** -1

        v_line.set_xdata([cur_x * coord_conversion_x, cur_x * coord_conversion_x])
        h_line.set_ydata([cur_y * coord_conversion_y, cur_y * coord_conversion_y])
        v_prof.set_xdata(z[:, int(cur_x)])
        h_prof.set_ydata(z[int(cur_y), :])

        fig.canvas.draw_idle()


# Optical Axis Plot:
main_ax.errorbar(tweezer_position[0], tweezer_position[1], tweezer_error, tweezer_error, capsize=10, capthick=4, lw=4,
                 color=uc_colors["gray"],
                 label="Mechanical Axis (Tweezer Beam)", zorder=1, fmt="+")

circle1 = plt.Circle((tweezer_position[0], tweezer_position[1]), 100e-6, color=uc_colors["gray"], fill=False, lw=4,
                     linestyle="--",
                     label="Spec. FOV around Mechanical Axis")

main_ax.add_patch(circle1)

cb = fig.colorbar(im, location="left", anchor=(-0.1, 0), label="Resolution [µm]", shrink=0.72)
cb.set_label("Resolution [m]", labelpad=10)
main_ax.legend(fontsize=22)
# main_ax.grid(alpha=0.3)
fig.canvas.mpl_connect('motion_notify_event', on_move)
plt.show()

# %% ---------------------------------------- PLOT DATAPOINTS -------------------------------------------------------
plt.rcParams.update({'font.size': 30})
from matplotlib.colors import Normalize

norm = Normalize(vmin=data["all_resolutions"].min(), vmax=3.5e-6)

fig, main_ax = plt.subplots(figsize=(23, 15))
sc = main_ax.scatter(data["all_coordinates"][0], data["all_coordinates"][1],
                     s=450,
                     # s=data["all_resolutions_error"] * 5 / (100e-6 / 16000),  # Just to display the error nicely...
                     c=data["all_resolutions"],
                     cmap="RdBu_r",
                     # vmin=0.775e-6,
                     # vmax=2.85e-6,
                     # vmax=0.5e-6,
                     # alpha=0.9,
                     label="PSF Fits",
                     norm=LogNorm(rsecondmin, 2.85e-6)
                     )

main_ax.set_xlabel('x Atomplane [µm]', labelpad=15)
main_ax.set_ylabel('y Atomplane [µm]', labelpad=15)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x * 1e6))
# ticks_cb = ticker.FuncFormatter(lambda x, pos: '{:.e}'.format(x * 1e6))
main_ax.xaxis.set_major_formatter(ticks)
main_ax.yaxis.set_major_formatter(ticks)

circle1 = plt.Circle((tweezer_position[0], tweezer_position[1]), 100e-6, color=uc_colors["gray"], fill=False, lw=6,
                     linestyle="--",
                     label="Specified FOV around Mechanical Axis (Tweezer Beam)")
# circle2 = plt.Circle((minimum.x[0], minimum.x[1]), 100e-6, color='#C61A27', fill=False, lw=6,
#                      label="FOV around Fitted Minimum")
#
main_ax.add_patch(circle1)
# main_ax.add_patch(circle2)

main_ax.errorbar(tweezer_position[0], tweezer_position[1], tweezer_error, tweezer_error, capsize=10, capthick=6, lw=6,
                 color=uc_colors["gray"],
                 label="Mechanical Axis (Tweezer Beam)", zorder=1, fmt="+")
main_ax.grid(alpha=0.3)
main_ax.set_ylim([0, camera_y_pixels * pixel_size / magnification])
main_ax.set_xlim([0, camera_x_pixels * pixel_size / magnification])
main_ax.invert_yaxis()
main_ax.set_aspect(1)
main_ax.legend(fontsize=22)

divider = make_axes_locatable(main_ax)
cax = divider.append_axes("left", size="3%", pad=2.8)

cb = fig.colorbar(sc, cax=cax, shrink=0.83)
cb.set_label("Resolution [m]", labelpad=15)
cax.yaxis.tick_left()
cax.yaxis.set_label_position('left')
plt.show()

# %% ---------------------------------------- QUIVER PLOT ----------------------------------------
angle_theta = []
for x, y, theta in zip(data["psf_x_stddev"], data["psf_y_stddev"], data["psf_angles"]):
    if x < y:
        angle_theta = np.append(angle_theta, theta - np.pi / 2) if theta > np.pi / 2 else np.append(angle_theta,
                                                                                                    theta + np.pi / 2)
    else:
        angle_theta = np.append(angle_theta, theta - np.pi)

aspect_ratio = []
major_axis = []
for h, value in enumerate(data["psf_y_stddev"]):
    if data["psf_y_stddev"][h] > data["psf_x_stddev"][h]:
        aspect_ratio = np.append(aspect_ratio, data["psf_y_stddev"][h] / data["psf_x_stddev"][h])
        major_axis = np.append(major_axis, data["psf_y_stddev"][h])
    else:
        aspect_ratio = np.append(aspect_ratio, data["psf_x_stddev"][h] / data["psf_y_stddev"][h])
        major_axis = np.append(major_axis, data["psf_x_stddev"][h])

fig, main_ax = plt.subplots(figsize=(23, 15))
sc = main_ax.quiver(data["all_coordinates"][0], data["all_coordinates"][1],
                    np.cos(angle_theta) * -major_axis, np.sin(angle_theta) * major_axis,
                    aspect_ratio, headaxislength=0, headlength=0, pivot='middle', width=0.008, scale=270,
                    cmap="Reds", clim=(1, 2)
                    )

main_ax.set_xlabel('x Atomplane [µm]', labelpad=15)
main_ax.set_ylabel('y Atomplane [µm]', labelpad=15)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x * 1e6))
ticks_cb = ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x * 1e6))
main_ax.xaxis.set_major_formatter(ticks)
main_ax.yaxis.set_major_formatter(ticks)

circle1 = plt.Circle((tweezer_position[0], tweezer_position[1]), 100e-6, color=uc_colors["gray"], fill=False, lw=6,
                     linestyle="--",
                     label="Specified FOV around Mechanical Axis (Tweezer Beam)")
# circle2 = plt.Circle((minimum.x[0], minimum.x[1]), 100e-6, color='#C61A27', fill=False, lw=6,
#                      label="FOV around Fitted Minimum")

main_ax.errorbar(tweezer_position[0], tweezer_position[1], tweezer_error, tweezer_error, capsize=10, capthick=6, lw=6,
                 color=uc_colors["gray"],
                 label="Mechanical Axis (Tweezer Beam)", zorder=1, fmt="+")
main_ax.add_patch(circle1)
# main_ax.add_patch(circle2)
main_ax.grid(alpha=0.3)
main_ax.set_ylim([0, camera_y_pixels * pixel_size / magnification])
main_ax.set_xlim([0, camera_x_pixels * pixel_size / magnification])
main_ax.invert_yaxis()
main_ax.set_aspect(1)
main_ax.legend(fontsize=22)

divider = make_axes_locatable(main_ax)
cax = divider.append_axes("left", size="3%", pad=2.8)

cb = fig.colorbar(sc, cax=cax, label="Resolution [µm]", shrink=0.83)
cb.set_label(r"PSF Ellipse Aspect Ratio $\epsilon$", labelpad=15)
cax.yaxis.tick_left()
cax.yaxis.set_label_position('left')
plt.show()

# %% ------------------------------------- Calculate Mean PSF Aspect Ratio in FOV -------------------------------------
# print(aspect_ratio)
aspect_in_FOV = np.median(aspect_ratio[k])
aspect_in_FOV_error = np.std(aspect_ratio[k]) / np.sqrt(len(aspect_ratio[k]))

print("eps_mean", aspect_in_FOV)
print("eps_mean error", aspect_in_FOV_error)

# ------------------------------------------------ PLOTLY POLYNOM PLOT -------------------------------------------------
fit_plot = px.imshow(polynom2d_third_deg([x_i, y_i], c_values), x=x, y=y, range_color=[0.6e-6, 2.4e-6])

fit_plot.update_layout(
    # yaxis=dict(scaleanchor="x", dtick=5e-6, scaleratio=1, range=[0, max(data["all_coordinates"][1])]),
    #                    xaxis=dict(range=[0, max(data["all_coordinates"][0])]),
    coloraxis_colorbar=dict(title="Resolution (PSF radius) [µm]", ticksuffix="m"),
    scene_aspectratio=dict(x=1.5, y=1, z=1),
    # plot_bgcolor="rgba(0,0,0,0)",
    paper_bgcolor='rgba(0,0,0,0)',
    title_text="PSF FOV Polynomial Fit",
    title_x=0.5,
    yaxis_dtick=100e-6, yaxis_ticksuffix="m",
    yaxis_tickfont=dict(size=15),
    #                    # xaxis_range=[0, max(data["all_coordinates"][0])],
    #                    # yaxis_range=[0, max(data["all_coordinates"][1])],
    xaxis_dtick=100e-6, xaxis_ticksuffix="m",
    xaxis_tickfont=dict(size=15),
    xaxis_title="Atomplane x [µm]",
    # yaxis_autorange="reversed",
    yaxis_title="Atomplane y [µm]"
)

fit_plot.write_html("FOV_Fit.html", include_plotlyjs="cdn", full_html=False,
                    include_mathjax='cdn')

# %% ------------------------------------------ PLOTLY DATAPOINTS FOV PLOT ------------------------------------------

x = np.linspace(0, camera_x_pixels * pixel_size / magnification, 500)
y = np.linspace(0, camera_y_pixels * pixel_size / magnification, 500)
x_i, y_i = np.meshgrid(x, y)
fig4 = px.imshow(polynom2d_third_deg([x_i, y_i], c_values), x=x, y=y, range_color=[0.6e-6, 3e-6],
                 # labels={'x': "x (atomplane) [µm]", 'y': "y (atomplane) [µm]", 'color': 'Fitted Resolution'}
                 )
# fig4.update_traces(opacity=0.5)

# Data
fig1 = px.scatter(data["all_coordinates"], x=data["all_coordinates"][0], y=data["all_coordinates"][1],
                  color=data["all_resolutions"], size=[0.1 for i in data["all_resolutions"]], size_max=16,
                  labels={'x': "x (atomplane) [µm]", 'y': "y (atomplane) [µm]",
                          'color': 'Resolution (PSF Radius) [µm]'},
                  # range_color=[0.6e-6, 3e-6],
                  color_continuous_scale="viridis"
                  )
fig1.update_traces(marker=dict(line=dict(width=0.001), opacity=0.95))

# Tweezer Circle:
t = np.linspace(0, 100, 5000)
x_circ, y_circ = 100e-6 * np.cos(t) + tweezer_position[0], 100e-6 * np.sin(t) + tweezer_position[1]
fig2 = px.scatter(x=x_circ, y=y_circ, size=np.array([2 for x in x_circ]))
fig2.update_traces(marker_size=4, marker_line_width=0, marker_color='#FDE74C')

# Fitted Circle:
x_circ, y_circ = 100e-6 * np.cos(t) + minimum.x[0], 100e-6 * np.sin(t) + minimum.x[1]
fig3 = px.scatter(x=x_circ, y=y_circ, size=np.array([2 for x in x_circ]))
fig3.update_traces(marker_size=4, marker_line_width=0, marker_color='#C61A27', hoverinfo='skip')

fig5 = go.Figure(data=fig1.data + fig2.data + fig3.data, layout=fig1.layout)

fig5.update_layout(yaxis=dict(scaleanchor="x", dtick=5e-6, scaleratio=1, range=[0, max(data["all_coordinates"][1])]),
                   xaxis=dict(range=[0, max(data["all_coordinates"][0])]),
                   coloraxis_colorbar=dict(title="Resolution (PSF radius) [µm]", ticksuffix="m"),
                   scene_aspectratio=dict(x=1.5, y=1, z=1),
                   plot_bgcolor="rgba(0,0,0,0)",
                   paper_bgcolor='rgba(0,0,0,0)',
                   title_text="PSF FOV Fit <br>"
                              "Distance between FOV Minimum and Tweezer Position {:.2f} µm => {:.2f} mrad".format(
                       distance * 1e6, distance_to_angle * 1e3),
                   title_x=0.5,
                   # title_y=0.98,
                   yaxis_dtick=100e-6, yaxis_ticksuffix="m",
                   yaxis_tickfont=dict(size=15),
                   xaxis_range=[0, max(data["all_coordinates"][0])],
                   yaxis_range=[0, max(data["all_coordinates"][1])],
                   xaxis_dtick=100e-6, xaxis_ticksuffix="m",
                   xaxis_tickfont=dict(size=15),
                   xaxis_title="Atomplane x [µm]",
                   yaxis_autorange="reversed",
                   yaxis_title="Atomplane y [µm]"
                   )

fig5.update_xaxes(showline=True, linewidth=1, linecolor='grey', gridcolor='lightgrey')
fig5.update_yaxes(showline=True, linewidth=1, linecolor='grey', gridcolor='lightgrey')

fig5.write_html("PSF_FOV_Datapoints.html", include_plotlyjs="cdn", full_html=False,
                include_mathjax='cdn')

# %% ------------------------------------------ OLD Interactive PSF 3D Plot ------------------------------------------
# fig1 = px.scatter_3d(
#     data["all_coordinates"], x=data["all_coordinates"][0], y=data["all_coordinates"][1], z=data["all_resolutions"],
#     #labels={"x":"","y":"","z":""},
#     labels={'x': "x (atomplane) [µm]", 'y': "y (atomplane) [µm]", 'z': 'Resolution (PSF Radius) [µm]'},
#     color=data["all_resolutions"],
#     range_color=[0.6e-6, 2e-6],
#     custom_data=[data["psf_filenames"]]
#     # title='<b>Mitutoyo Imaging Setup Characterization: Nanoholes imaged across CCD FOV</b><br>'
#     #       # '  Theoretical Resolution: (0.61*{:.0f}nm/NA) = {:.2f} µm | '
#     #       '  Magnification = {:.1f}<br>'
#     #       '  Median Resolution = {:.2f} µm | Min. Resolution = {:.2f} µm <br>'
#     #       '  Objective FOV on CCD: 1070 x 800  |  Total CCD FOV: 1288 x 964'
#     #       .format(# wavelength*1e9, resolution_theo * 1e6,
#     #               magnification, np.median(all_resolutions) * 1e6, np.min(all_resolutions) * 1e6),
#     # opacity=0.7
# )
#
# fig1.update_traces(hovertemplate="<br>".join([
#     "Resolution (PSF radius): %{z:.2s}m",
#     "[x, y] = [%{x:.2s}m, %{y:.2s}m]",
#     """<a href="https://lithium6hd.github.io/HQA/ObjectiveCharacterization/PSFs/%{customdata[0]}.html">PSF Plot</a>"""
# ]),
#     marker_size=5,
#     marker_line_width=0
# )
#
# # FOV Circle:
# t = np.linspace(0, 100, 5000)
# x_circ, y_circ, z_circ = 100e-6 * np.cos(t) + tweezer_position[0], 100e-6 * np.sin(t) + tweezer_position[1],\
#                          [0.6239e-6 for x in t]  # 0.6e-6+0.000000002*t
# fig2 = px.scatter_3d(x=x_circ, y=y_circ, z=z_circ, size=np.array([2 for x in x_circ]), range_color=[0.6e-6, 2e-6])
# fig2.update_traces(marker_size=4, marker_line_width=0, marker_color="black")
#
# fig3 = go.Figure(data=fig1.data + fig2.data, layout=fig1.layout)
#
# fig3.add_annotation(text=r'$\textsf{ Rayleigh criterion: }0.61 \cdot \frac{671\text{ nm}}{\text{NA}}=0.62 \text{µm}$',
#                     align='left',
#                     showarrow=False,
#                     xref='paper',
#                     yref='paper',
#                     font=dict(size=18),
#                     x=0.01,
#                     y=0.95,
#                     borderwidth=0)
#
# fig3.update_layout(yaxis=dict(scaleanchor="x", dtick=5e-6, scaleratio=1, range=[0, max(data["all_coordinates"][1])]),
#                    xaxis=dict(range=[0, max(data["all_coordinates"][0])]),
#                    scene=dict(zaxis=dict(range=[0.59e-6, 2e-6])),
#                    coloraxis_colorbar=dict(title="Resolution (PSF radius) [µm]", ticksuffix="m"),
#                    scene_aspectratio=dict(x=1.5, y=1, z=1),
#                    plot_bgcolor="rgba(0,0,0,0)",
#                    paper_bgcolor='rgba(0,0,0,0)',
#                    title_text="PSF FOV measurement",
#                    title_x=0.5
#                    )
#
# fig3.update_scenes(xaxis_autorange="reversed", xaxis_dtick=100e-6, xaxis_ticksuffix="m",
#                    xaxis_tickfont=dict(size=15),
#                    yaxis_dtick=100e-6, yaxis_ticksuffix="m",
#                    yaxis_tickfont=dict(size=15),
#                    zaxis_ticksuffix="m",
#                    zaxis_tickfont=dict(size=15)
#                    )
#
# # # py.offline.plot(fig, filename="plotly version of an mpl figure.html")
# # fig3.show()
# fig3.write_html(working_directory + "PSF_Resolution_Interactive_2.html", include_plotlyjs="cdn", full_html=False,
#                 include_mathjax='cdn')
