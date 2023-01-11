import time
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import numpy as np
from simple_pyspin import Camera
from PIL import Image
from astropy.modeling import models, fitting
from pylablib.devices import Thorlabs

wavelength = 671e-9  # [m]
numerical_aperture = 0.656  # Special Optics: 0.656, Mitutoyo: 0.3
pixel_size = 2.4e-6  # [m], Blackfly: 2.4µm Chameleon: 3.7µm
camera_x_pixels, camera_y_pixels = 5472, 3648
resolution_theo = 0.61 * wavelength / numerical_aperture  # Rayleigh Resolution
focus_error = 1e-6  # Error of the nanohole position along the optical axis
focal_length_obj = 18.8e-3  # Focal length of objective
focal_length_lens = 300e-3  # The lens determining the magnification used to focus light on the camera
effective_focal_length = ((1 / focal_length_obj + 1 / focal_length_lens) ** (-1))
magnification = focal_length_lens / focal_length_obj
magnification_error = effective_focal_length / (effective_focal_length - focal_length_obj) ** 2 * focus_error

# Set the working directory
working_directory = "files/22072022/4/"

# Record only a single PSF without stage movement and without writing results in files?
single_psf = True if input("Record only a single PSF without writing results in files? (y/n)") == "y" else False

# Set False if only the analysis processing should be run
record_mode = False if input("Only run fitting algorithm without recording? (y/n)") == "y" else True

# If the data still has to be calculated or can be read from file:
psf_fitting_true = True if input("Run fitting after recording? (y/n)") == "y" else False


if single_psf:
    line_reset_steps = 0
    matrix = (list(range(0, 8, 4)), list(range(0, 50, 25)))  # Stage won't move, but a single value is still needed
else:
    line_reset_steps = -425  # <----------- SET STEPS THE STAGE TRAVELS BACK TO RESET A ROW ALONG X --------------------
    matrix = (list(range(0, 3765, 75)),  # X Axis -> Second Value: Max. steps in a row | Third value: no. of images
              list(range(0, 1385, 35)))  # Y Axis -> Second Value: Max. steps in a line | Third value: no. of images

steps_x_y = [-1 * np.diff(matrix[0]), np.diff(matrix[1])]  # Array that contains the steps for each loop
number_of_pictures = (len(steps_x_y[0]) * len(steps_x_y[1]))  # Total number of pictures for ets. time calculation

print("Recording matrix: {} x {} steps = {} images"
      .format(len(steps_x_y[0]), len(steps_x_y[1]), number_of_pictures))

imgs = []
darks = []

# # %% Testing the camera:
#
# with Camera() as cam:
#     cam.start()  # Start recording
#     test_images = []
#
#     # ------------------------------- Images recording --------------------------------
#     # input("Press Enter to start image recording...")
#
#     # Initialize plot:
#     plt.figure(figsize=(10, 2.5))
#
#     for i in range(2):
#         test_images.append(cam.get_array())
#         # time.sleep(2)
#
#         # Show the test images:
#         plt.subplot(1, len(range(2)), i+1)
#         plt.imshow(test_images[i])
#         plt.xlabel('x [pixel]')
#         plt.ylabel('y [pixel]')
#         plt.title("Test " + str(i+1))
#
#     plt.colorbar(label="Pixel intensity [a.u.]")
#     plt.tight_layout()
#     plt.show()
#     cam.stop()
#
# sys.exit(0)

# ------------------------------- RECORDING PART  --------------------------------
if record_mode:
    with Camera() as cam, Thorlabs.KinesisPiezoMotor("97101928") as stage:
        # ------------------------------- Dark images recording --------------------------------

        dark_record = input("Dark image recording? (y/n)")

        if dark_record == "y":
            cam.start()  # Start recording
            for i in range(4):
                darks.append(cam.get_array())
                time.sleep(0.2)

            dark_med = sum(darks) / len(darks)
            dark_med_img = Image.fromarray(dark_med.astype(np.uint8))
            dark_med_img.save("Images/Dark_med.bmp")
            cam.stop()  # Stop recording

        # ------------------------------- Nanohole Images recording --------------------------------

        input("Press Enter to continue with PSF recording...")

        stage.enable_channels(1)
        stage.enable_channels(2)

        # Set starting position to zero (no movement):
        stage.set_position_reference(channel=1)
        stage.set_position_reference(channel=2)
        cam.start()  # Ready camera

        first = True
        l_it = 0  # image number iteration
        row_it = 0  # row iteration

        for j in steps_x_y[0]:
            if first:
                start_time = time.time()
            else:
                # Move X coordinate (j) in every loop but the first
                stage.move_by(j, channel=1)  # Positive step -> Down
                stage.wait_move()

            for k in steps_x_y[1]:

                # ------------------------------- Record image  --------------------------------
                time.sleep(0.4)  # Wait for stage to stabilize
                newimage = cam.get_array()
                img = Image.fromarray(newimage)
                img.save(working_directory+"Images/PSF_" + str(l_it) + ".bmp")
                # time.sleep(0.1)  # Wait for image recording

                l_it += 1  # Increase number of recorded images by 1

                print("Process: {:.2f}%".format(100 * l_it / number_of_pictures))

                # ------------------------------- Move Y coordinate (k) one step --------------------------------
                if not single_psf:
                    time.sleep(0.2)
                    stage.move_by(k, channel=2)
                    stage.wait_move()
                    time.sleep(0.1)

            # --------------------------- Move back in X direction  ------------------------
            time.sleep(3)

            # At the end, move by [line_reset_steps] along X to the start of the row:
            if not single_psf:
                stage.enable_channels(2)
                stage.move_by(line_reset_steps, channel=2)
                stage.wait_move()

            print("Finished row {} of {}.".format(row_it+1, len(steps_x_y[0])+1))

            row_it += 1

            time.sleep(3)

            if first:
                end_time = time.time()  # Measure time for single row to approximate total time
                first = False

            print("Approx. finish in ~{:.1f}min + {:.0f}min Analysis"
                  .format((end_time-start_time)*(len(steps_x_y[0])-row_it)/60, 0.008 * number_of_pictures))

        cam.stop()  # stop recording
        stage.close()


# %%
# ------------------------------- PROCESSING PART  --------------------------------

def psf_fitting(directory, px_size, magnification, magnification_err):
    """
    Gives two resolutions for a single nanohole image and their coordinates in a list [[x], [y]].
    The last three return values are only useful if the fit model is gaussian
    :param directory: Directory of image files incl. PSF and Dark pictures [String]
    :param px_size: Pixel size [m]
    :param magnification: Imaging magnification
    :param magnification_err: Uncertainty of the magnification
    :return (Resolution, Resolution error, coordinates, psf_names, angles, x_stddev, y_stddev)
    """

    resolutions = []
    resolution_errors = []
    coordinates = [[], []]
    psf_names = []
    angles = []
    x_stddev = []
    y_stddev = []

    offset_list = [0]  # Upper one too dark , hole_offset]  # Upper hole is darker: -hole_offset

    fit_model = "Gaussian"

    if fit_model == "Airy":
        # Fit model:
        f_init = models.AiryDisk2D(radius=4.75, fixed={"radius": True}
                                   # bounds={"radius": (4.3, 4.8), "x_0": (10, 11), "y_0": (10, 11)}
                                   )  #
    else:  # fit_model == "Gaussian":
        f_init = models.Gaussian2D()

    # Declare what fitting function you want to use
    fit_f = fitting.LevMarLSQFitter()

    for offset in offset_list:

        for filename in sorted(os.listdir(directory)):

            if filename.startswith("Dark"):

                dark = Image.open(directory + '/' + filename)
                dark = np.array(dark)

            elif filename.startswith("PSF_"):

                print("Evaluating ", filename)
                image = Image.open(directory + '/' + filename)
                image = np.array(image)
                med = np.median(image)
                image = image - med

                image_new = image  # - dark

                image_new[1794, 2717] = 0  # Hot pixel
                image_new[2311, 4109] = 0  # Hot pixel
                image_new[1494, 1252] = 0  # Hot pixel
                image_new[512, 1825] = 0  # Hot pixel

                # Find maximum of image, i.e. the brightest nanohole, i.e. the center of the PSF
                cents_i = np.where(image_new == np.max(image_new))

                xc_i = int(cents_i[1][0])
                yc_i = int(cents_i[0][0]) - offset
                # print("x: " + str(xc_i) + ", y: " + str(yc_i))

                # Cut out smaller box around PSF
                bb = 12
                box = image_new[yc_i - bb:yc_i + bb, xc_i - bb:xc_i + bb]
                # box = np.flip(box, 1)  # Display PSF picture according to real picture
                yp_i, xp_i = box.shape

                # Show PSF:
                # plt.imshow(box)
                # plt.xlabel('x [pixel]')
                # plt.ylabel('y [pixel]')
                # plt.title("Test Image at " + str(xc_i) + "," + str(yc_i))
                # plt.colorbar(label="Pixel intensity [a.u.]")
                # plt.show()

                # Generate grid of same size like box to put the fit on
                y_i, x_i, = np.mgrid[:yp_i, :xp_i]

                # Check for border issues (Blackfly Camera pixel: 5472px x 3648px), Middle point: 2736,1824
                if yc_i > bb and camera_y_pixels - yc_i > bb and xc_i > bb and camera_x_pixels - xc_i > bb:

                    # Fit process
                    f = fit_f(f_init, x_i, y_i, box)

                    # Covariance Matrix
                    cov = fit_f.fit_info["param_cov"]
                    param_std = np.sqrt(np.diag(cov)) * px_size / magnification

                    if fit_model == "Airy":

                        resolution_i = f.radius * px_size / magnification
                        resolution_error = param_std[0]

                    else:  # fit_model == "Gaussian":
                        resolution_i = (f.y_stddev + f.x_stddev) / 2 * px_size / magnification * 2.66
                        # Factor 2.66 to convert Gauss to Airy (see master thesis Micha)

                        ellipse_aspect_error = np.mean([abs(f.x_stddev * px_size / magnification * 2.66 - resolution_i),
                                                        abs(f.y_stddev * px_size / magnification * 2.66 - resolution_i)])
                        resolution_error = (ellipse_aspect_error +
                                            ((param_std[3] + param_std[4]) / 2) +
                                            (resolution_i * magnification_err / magnification))

                        angles = np.append(angles, f.theta)
                        x_stddev = np.append(x_stddev, f.x_stddev)
                        y_stddev = np.append(y_stddev, f.y_stddev)

                    # Throw away bad fits:
                    if 5e-6 > resolution_i > 0.4e-6 and cov is not None:

                        # # Plot a single PSF fit for diagnosis
                        # plt.figure(figsize=(8, 2.5))
                        # plt.subplot(1, 3, 1)
                        # plt.imshow(box)
                        # plt.title("Data")
                        # plt.subplot(1, 3, 2)
                        # plt.imshow(f(x_i, y_i))
                        # plt.title("Gaussian 2D Model")
                        # plt.subplot(1, 3, 3)
                        # plt.imshow(box - f(x_i, y_i))
                        # plt.title("Residual")
                        # plt.suptitle(
                        #     "x,y = [" + str(xc_i) + "," + str(yc_i) + "] Fit to nanohole | Airy radius = {:.2f} µm"
                        #     .format(resolution_i * 1e6), fontsize=13, y=0.96)
                        # plt.tight_layout()
                        # # plt.savefig("Nanohole Fit.png", format="png", dpi=300)
                        # plt.show()

                        # Write results in a list
                        resolutions.append(resolution_i)
                        resolution_errors.append(resolution_error)

                        # Also add the coordinates for each resolution:
                        coordinates = np.append(coordinates, [[xc_i], [yc_i]], axis=1)
                        psf_names = np.append(psf_names, filename[:-4])
                        angles = np.append(angles, f.theta)
                        x_stddev = np.append(x_stddev, f.x_stddev)  # max(f.y_stddev, f.x_stddev))
                        y_stddev = np.append(y_stddev, f.y_stddev)

                        print("Resolution (Radius to first minimum): {:.2f} +- {:.2f}µm at [{:.0f},{:.0f}]"
                              .format(resolution_i * 1e6, resolution_error * 1e6, xc_i, yc_i))

    return resolutions, resolution_errors, coordinates, psf_names, angles, x_stddev, y_stddev


# %%


if psf_fitting_true:

    # Fitting:
    all_resolutions, all_resolutions_error, all_coordinates, psf_filenames, psf_angles, psf_x_stddev, psf_y_stddev = \
        psf_fitting(working_directory + "Images", pixel_size, magnification, magnification_error)

    if not single_psf:
        # Writing the results in text files:
        data = {"all_resolutions": all_resolutions,
                "all_resolutions_error": all_resolutions_error,
                "all_coordinates": all_coordinates,
                "psf_filenames": psf_filenames,
                "psf_angles": psf_angles,
                "psf_x_stddev": psf_x_stddev,
                "psf_y_stddev": psf_y_stddev,
                }
        for col_name in data:
            np.savetxt(working_directory + "Results/{}.txt".format(col_name), data[col_name], fmt='%s')

# else:
#     # For reading the results from already produced files
#     data = {}
#     for filename in sorted(os.listdir(working_directory + "Results")):
#         if filename.endswith("txt"):
#             if filename.startswith("psf_filenames"):
#                 data[filename[:-4]] = np.loadtxt(working_directory + "Results/" + filename, dtype='U')
#             else:
#                 data[filename[:-4]] = np.loadtxt(working_directory + "Results/" + filename)

# # Convert Coordinates: image plane [px] to object plane [µm] conversion:
# data["all_coordinates"] = data["all_coordinates"] * pixel_size / magnification
