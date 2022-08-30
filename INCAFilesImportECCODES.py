# ---------------------------------------------------------------------------
# READ inca files Format, Plot, Reproject, interpolate
# ---------------------------------------------------------------------------

# CLEAR ALL
for v in dir():
    exec('del '+ v)
    del v
###########

# Import
import os
import datetime
import numpy as np

from ptype_functions import *
from IncaGribImport import IncaGribImporter

from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

# -----------------------------------------------------------------------------
# ------------------------------ # START # ------------------------------------
# -----------------------------------------------------------------------------
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
dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs'

# GRIB selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']
inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
nwc_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '

# ---------------------------------------------------------------------------
# Import ZS Snow level
filename_ZS = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_zs)
importer = IncaGribImporter()
incaDictionary_ZS = importer.retrieve_grib_data(filename=filename_ZS, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_ZS = inca_dictionary_to_3Dmatrix(incaDictionary_ZS)

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

# Clean
del (importer)
del (incaDictionary_ZS)
del (incaDictionary_TT)
del (incaDictionary_TG)

# --------------------------------------------------------------------------
# LOAD INCA TOPOGRAPHY
inca_topo_grid = np.loadtxt(topoFilename)
inca_topo_grid = inca_topo_grid[::-1, :]  # Reorientation

# ---------------------------------------------------------------------------
# LOAD netCDF FILE

# Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc['projection'] = nwc_projectionString
metadata_nwc['cartesian_unit'] = 'km'
print('done!')

# ---------------------------------------------------------------------------
# Reproject
# r_nwc has an extra dimension for members r_nwc[index,i,:,:] let's use the first one for now
inca_reprojected_ZS, inca_metadata_reprojected = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca,
                                                                 metadata_nwc)
inca_reprojected_TT, _ = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TG, _ = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
topo_reprojected, topo_metadata_reprojected = reproject_grids(np.array([inca_topo_grid]), r_nwc[0, 0, :, :],
                                                              metadata_inca, metadata_nwc)



# PLOT TEST
# filenames = []
# filenames.append(worker(inca_reprojected_ZS, inca_metadata_reprojected, 0))
# filenames.append(worker(inca_reprojected_TT, inca_metadata_reprojected, 0))
# filenames.append(worker(inca_reprojected_TG, inca_metadata_reprojected, 0))


# ---------------------------------------------------------------------------
# PLOT# PLOT# PLOT# PLOT# PLOT# PLOT
# normalize between 0-1 and plot probability
# minVal = np.min(R_inca)
# maxVal = np.max(R_inca)
# diff = np.max(R_inca) - np.min(R_inca)
# a = np.copy(R_inca[0, :, :])
#
# for i in range(a.shape[0]):
#     a[i, :] = (a[i, :] - minVal) / diff
# a_r_reprojected, metadata_reprojected = reproject_grids(np.array([a]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
# filenames.append(worker(a_r_reprojected, metadata_reprojected, 0))

# ---------------------------------------------------------------------------
# Interpolations
inca_interpolation_ZS = grid_interpolation(inca_reprojected_ZS[0], inca_reprojected_ZS[1], timeStep=5)
inca_interpolation_TT = grid_interpolation(inca_reprojected_TT[0], inca_reprojected_TT[1], timeStep=5)
inca_interpolation_TG = grid_interpolation(inca_reprojected_TG[0], inca_reprojected_TG[1], timeStep=5)



####
# from scipy import interpolate
#
# a = np.zeros((3,3))
# b = np.zeros((3,3))
# b[0,:]=5
# b[1,:]=10
# b[2,:]=15
#
# interPoints = np.arange(0, (60 + 10), 10)
# interpolationGrid = np.zeros((a.shape[0], a.shape[1], len(interPoints)))
# interpolationGrid[:, :, :] = np.nan
# for i in range(a.shape[1]):
#     for j in range(b.shape[1]):
#         interpolationFunction = interpolate.interp1d([0, 60], [a[i, j], b[i, j]])
#         interpolationGrid[i, j, :] = interpolationFunction(interPoints)
# print('Done')

##############
# timeStep = 7
# td = 50
# interPoints[0]
# t1 = ((timeStep-td)/(timeStep-1)) * a + ((0-1)/(timeStep-1)) * b
#
# y = y1 + ((x – x1) / (x2 – x1)) * (y2 – y1)
# (0,0) , (60,15)
# 0 + ((td - 0) / (60-0)) * (15 - 0)
#
# a + ((interPoints[-1] - a) / (interPoints[-1] - a)) * (b - a)
#
# interpolationGrid[:,:,6]

# inca_interpolation_ZS[210, 210, 0] == inca_reprojected_ZS[0, 210, 210]
# inca_interpolation_ZS[210, 210, 12] == inca_reprojected_ZS[1, 210, 210]
# ~np.isnan(inca_interpolation_ZS[650, 660, 0])
# np.isnan(inca_interpolation_ZS[650, 661, 0])
# np.isnan(inca_interpolation_ZS[651, 660, 0])

# --------------------------------------------------------------------------
# calculate the interpolation points (index selection between inca and pysteps) using the timestamps as index from the pysteps
# How to determine the starting point , we should only use a index seq and determine the right grid by timestamp in metadata_nwc
# timeBase = 60
# timeStep = 5
#
# # Loop over
# for i in range(len(metadata_nwc['timestamps'])):
#     incaGrid = inca_interpolation[0][:, :, (i + common_index)]  # interpolation 0
#     precipGrid = r_nwc[0, i - common_index, :, :]  # CONSIDERING ONLY THE FIRST MEMBER! = 0

# -------------------------------------------------------
# Diagnosis

## HERE find the index of 15:00 in the inca file and netCDF
# metadata_nwc['timestamps'][10]  # 15:00

# calculate
interp_idx = 0
precip_idx = 0

# tests
if metadata_nwc['timestamps'][precip_idx].minute == 0:
    interp_idx = -1
else:
    interp_idx = np.where(np.arange(0, 65, 5) == metadata_nwc['timestamps'][precip_idx].minute)[0][0]
print("Precipitation index ", precip_idx, " timestamp ", metadata_nwc['timestamps'][precip_idx])
print("Match to interpolation GRIB index: ", interp_idx)


# calculate precipitation type result
ptype = calculate_precip_type(incaZnow=inca_interpolation_ZS[:, :, interp_idx],
                              incaTemp=inca_interpolation_TT[:, :, interp_idx],
                              incaGroundTemp=inca_interpolation_TG[:, :, interp_idx],
                              precipGrid=r_nwc[0, precip_idx, :, :],
                              topographyGrid=topo_reprojected[0, :, :])

# PLOT
# filenames = []
# filenames.append(worker(np.array([ptype]), metadata_nwc, 0, startdate, dir_gif))


