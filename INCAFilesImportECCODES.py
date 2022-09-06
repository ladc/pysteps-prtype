# ---------------------------------------------------------------------------
# READ inca and pysteps files.
# Format, Reproject, interpolate and plot
# ---------------------------------------------------------------------------

# Import
import os
import datetime
import numpy as np

from ptype_functions import *
from IncaGribImport import IncaGribImporter

from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

# --------------------------------------------------------------------------
# Parameters

# Date
startdate = datetime.datetime.strptime('202201090200', "%Y%m%d%H%M")

# File paths
gribPath = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\09\basic'
fn_zs = '_ZS_FC_INCA.grb'
fn_tt = '_TT_FC_INCA.grb'
fn_tg = '_TG_FC_INCA.grb'
netCDF4Filename = r'C:\Users\talen\Desktop\Student JOB\Data\pysteps_blended_nowcast_20220108-09\blended_nowcast_202201090205.nc'
topoFilename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-TOPO\inca_topo.asc'
dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs\membersTest'

# GRIB import selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']

# Projection strings
inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
nwc_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '

# Time period settings (these are default values in all function parameters)
timeBase = 60
timeStep = 5

# ---------------------------------------------------------------------------
# Import ZS Snow level
filename_ZS = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_zs)
importer = IncaGribImporter()
incaDictionary_ZS = importer.retrieve_grib_data(filename=filename_ZS, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_ZS = inca_dictionary_to_3Dmatrix(incaDictionary_ZS)

del (importer)

# ---------------------------------------------------------------------------
# Import TT Temperature
filename_TT = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_tt)
importer = IncaGribImporter()
incaDictionary_TT = importer.retrieve_grib_data(filename=filename_TT, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TT = inca_dictionary_to_3Dmatrix(incaDictionary_TT)

# Transform Kelvin to Celsius
R_inca_TT[:, :, :] = R_inca_TT[:, :, :] - 273.15

del (importer)

# ---------------------------------------------------------------------------
# Import TG Ground temperature
filename_TG = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_tg)
importer = IncaGribImporter()
incaDictionary_TG = importer.retrieve_grib_data(filename=filename_TG, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TG = inca_dictionary_to_3Dmatrix(incaDictionary_TG)

# Transform Kelvin to Celsius
R_inca_TG[:, :, :] = R_inca_TG[:, :, :] - 273.15

# ---------------------------------------------------------------------------
# Build inca metadata
metadata_inca = {}
metadata_inca['projection'] = inca_projectionString
metadata_inca['xpixelsize'] = incaDictionary_ZS['Metadata']['DxInMetres']
metadata_inca['ypixelsize'] = incaDictionary_ZS['Metadata']['DyInMetres']
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
metadata_inca['yorigin'] = 'upper'

# --------------------------------------------------------------------------
# LOAD INCA TOPOGRAPHY (this might be a different file format in the future)
inca_topo_grid = np.loadtxt(topoFilename)
inca_topo_grid = inca_topo_grid[::-1, :]  # Reorientation

# --------------------------------------------------------------------------
# Clean
del importer, incaDictionary_ZS, incaDictionary_TT, incaDictionary_TG

# ---------------------------------------------------------------------------
# LOAD netCDF file
# Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc['projection'] = nwc_projectionString
metadata_nwc['cartesian_unit'] = 'km'
print('netCDF4 load done!')

# --------------------------------------------------------------------------
# Reproject
# r_nwc has an extra dimension for members r_nwc[index,i,:,:] let's use the first one for now
inca_reprojected_ZS, inca_metadata_reprojected = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca,
                                                                 metadata_nwc)
inca_reprojected_TT, _ = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TG, _ = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
topo_reprojected, _ = reproject_grids(np.array([inca_topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
print('Reprojection done!')

# --------------------------------------------------------------------------
# Diagnosis (Single selection)
# --------------------------------------------------------------------------

# Interpolations (timeStep=5 and timeBase=60 are default in this function)
inca_interpolation_ZS = grid_interpolation(inca_reprojected_ZS[0], inca_reprojected_ZS[1], timeStep=timeStep,
                                           timeBase=timeBase)
inca_interpolation_TT = grid_interpolation(inca_reprojected_TT[0], inca_reprojected_TT[1], timeStep=timeStep,
                                           timeBase=timeBase)
inca_interpolation_TG = grid_interpolation(inca_reprojected_TG[0], inca_reprojected_TG[1], timeStep=timeStep,
                                           timeBase=timeBase)
print('Interpolation done!')

# Create a new function to create full interpolation list with a selected period
# AFTER interpolation we don't need the reprojected data anymore
del inca_reprojected_ZS, inca_reprojected_TT, inca_reprojected_TG

# interp_idx = 0
members_ts_index = 0

# tests #
if metadata_nwc['timestamps'][members_ts_index].minute == 0:
    interp_idx = -1
else:
    interp_idx = \
    np.where(np.arange(0, (timeBase + timeStep), timeStep) == metadata_nwc['timestamps'][members_ts_index].minute)[0][
        0]
print("Precipitation member index ", members_ts_index, " timestamp ", metadata_nwc['timestamps'][members_ts_index])
print("Match to interpolation GRIB index: ", interp_idx)

# Members Mean matrix.
r_nwc_mean = calculate_members_mean(r_nwc[:, members_ts_index, :, :])  # Input: an Array of member grids (3D matrix)

# calculate precipitation type result
ptype = calculate_precip_type(incaZnow=inca_interpolation_ZS[interp_idx, :, :],
                              incaTemp=inca_interpolation_TT[interp_idx, :, :],
                              incaGroundTemp=inca_interpolation_TG[interp_idx, :, :],
                              precipGrid=r_nwc_mean,
                              topographyGrid=topo_reprojected[0, :, :])

# PLOT
plot_ptype(np.array([ptype]), inca_metadata_reprojected, 0, metadata_nwc['timestamps'][members_ts_index], dir_gif)

# filenames = []
# filenames.append(plot_ptype(np.array([ptype]), metadata_nwc, 0, startdate, dir_gif))


# --------------------------------------------------------------------------
# Diagnosis (Loop over available PYSTEPS timestamps)
# --------------------------------------------------------------------------

# Create a timestamp index array for INCA interpolation matrix
inca_timestamps = create_timestamp_indexing(inca_reprojected_ZS.shape[0], startdate, timeStep=timeStep,
                                            timeBase=timeBase)

# Convert metadata_nwc['timestamps'] to datetime
nwc_timestamps = [datetime.datetime.strptime(ts.strftime('%Y%m%d%H%M'), '%Y%m%d%H%M') for ts in
                  metadata_nwc['timestamps']]

# Inca message index
message_idx = 1

# Calculate interpolation matrices
inca_interpolations_ZS, timestamps_idxs = generate_inca_interpolations(inca_reprojected_data=inca_reprojected_ZS,
                                                                       inca_timestamps=inca_timestamps,
                                                                       nwc_timestamps=nwc_timestamps,
                                                                       timeStep=timeStep,
                                                                       timeBase=timeBase)
inca_interpolations_TT, _ = generate_inca_interpolations(inca_reprojected_data=inca_reprojected_TT,
                                                         inca_timestamps=inca_timestamps,
                                                         nwc_timestamps=nwc_timestamps,
                                                         timeStep=timeStep,
                                                         timeBase=timeBase)
inca_interpolations_TG, _ = generate_inca_interpolations(inca_reprojected_data=inca_reprojected_TG,
                                                         inca_timestamps=inca_timestamps,
                                                         nwc_timestamps=nwc_timestamps,
                                                         timeStep=timeStep,
                                                         timeBase=timeBase)

# Clean (After interpolation, we don't need the reprojected data anymore)
del inca_reprojected_ZS, inca_reprojected_TT, inca_reprojected_TG

# Members mean precipitation type over time
print("Calculate precipitation type over members mean...")
filenames = []

for ts in range(len(timestamps_idxs)):
    print("Calculating precipitation type at: ", str(timestamps_idxs[ts]))

    # Members Mean matrix
    r_nwc_mean = calculate_members_mean(r_nwc[:, ts, :, :])  # Input: Array of member grids

    # calculate precipitation type result for the members mean
    ptype_mean = calculate_precip_type(incaZnow=inca_interpolations_ZS[ts, :, :],
                                       incaTemp=inca_interpolations_TT[ts, :, :],
                                       incaGroundTemp=inca_interpolations_TG[ts, :, :],
                                       precipGrid=r_nwc_mean,
                                       topographyGrid=topo_reprojected[0, :, :])
    # Plot
    # print('Plotting members mean precipitation type...')
    filenames.append(plot_ptype(ptype_mean, inca_metadata_reprojected, ts, startdate, dir_gif))
    # here could be a list to store the results

# Calculate precipitation type per member

dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs\member0'
print("Calculate precipitation type per members...")
filenames = []

for mem in range(r_nwc.shape[0]):
    # TEST
    mem = 0
    for ts in range(len(timestamps_idxs)):
        print("Calculating precipitation type at: ", str(timestamps_idxs[ts]), " for member : ", mem)
        # Calculate precipitation type for each member
        ptype_mem = calculate_precip_type(incaZnow=inca_interpolations_ZS[ts, :, :],
                                          incaTemp=inca_interpolations_TT[ts, :, :],
                                          incaGroundTemp=inca_interpolations_TG[ts, :, :],
                                          precipGrid=r_nwc[mem, ts, :, :],
                                          topographyGrid=topo_reprojected[0, :, :])
        # Plot
        filenames.append(plot_ptype(ptype_mem, inca_metadata_reprojected, ts, startdate, dir_gif))


# ------- SAVE

# --------------------------------------------------------------------------
# Members MEAN precipitation type over time
#
# print("Calculate precipitation type over members mean...")
# filenames = []
# for ts in range(len(timestamps_idxs)):
#     print("Calculating precipitation type at: ", str(timestamps_idxs[ts]))
#
#     # Members Mean matrix
#     r_nwc_mean = calculate_members_mean(r_nwc[:, ts, :, :])  # Input: Array of member grids
#
#     # calculate precipitation type result for the members mean
#     ptype_mean = calculate_precip_type(incaZnow=inca_interpolations_ZS[ts, :, :],
#                                        incaTemp=inca_interpolations_TT[ts, :, :],
#                                        incaGroundTemp=inca_interpolations_TG[ts, :, :],
#                                        precipGrid=r_nwc_mean,
#                                        topographyGrid=topo_reprojected[0, :, :])
#     # Plot
#     # print('Plotting members mean precipitation type...')
#     filenames.append(plot_ptype(ptype_mean, inca_metadata_reprojected, ts, timestamps_idxs[ts], dir_gif))
#     # here could be a list to store the results
#
# # Build gif
# kargs = {'duration': 0.4}
# with imageio.get_writer(os.path.join(dir_gif, r'INCA_mem_mean_%s.gif' % (startdate.strftime('%Y%m%d%H%M'),)), mode='I',
#                             **kargs) as writer:
#         for filename in filenames:
#             image = imageio.imread_v2(os.path.join(dir_gif, filename))
#             writer.append_data(image)
#
# # Close gif writer
# writer.close()
#
# # Remove temporary files
# for filename in set(filenames):
#     os.remove(os.path.join(dir_gif, filename))


# --------------------------------------------------------------------------
# Calculate precipitation type per member
#
# print("Calculate precipitation type per member...")
# filenames = []
# for mem in range(members):
#     for ts in range(len(timestamps_idxs)):
#         print("Calculating precipitation type at: ", str(timestamps_idxs[ts]), " for member : ", mem)
#         # Calculate precipitation type for each member
#         ptype_mem = calculate_precip_type(incaZnow=inca_interpolations_ZS[ts, :, :],
#                                           incaTemp=inca_interpolations_TT[ts, :, :],
#                                           incaGroundTemp=inca_interpolations_TG[ts, :, :],
#                                           precipGrid=r_nwc[mem, ts, :, :],
#                                           topographyGrid=topo_reprojected[0, :, :])
#         # Plot
#         filenames.append(plot_ptype(ptype_mem, inca_metadata_reprojected, ts, timestamps_idxs[ts], dir_gif))
#
#     # Build gif
#     kargs = {'duration': 0.4}
#     with imageio.get_writer(os.path.join(dir_gif, (r'INCA_mem_' + str(mem) + '_' + startdate.strftime('%Y%m%d%H%M') + '.gif')), mode='I', **kargs) as writer:
#         for filename in filenames:
#             image = imageio.imread_v2(os.path.join(dir_gif, filename))
#             writer.append_data(image)
#
#     # Close gif writer
#     writer.close()
#
#     # Remove temporary files
#     for filename in set(filenames):
#         os.remove(os.path.join(dir_gif, filename))

