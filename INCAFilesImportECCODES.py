# ---------------------------------------------------------------------------
# READ inca files Format, Plot, Reproject, interpolate
# ---------------------------------------------------------------------------

# CLEAR ALL
for v in dir():
    exec('del ' + v)
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


# ------------------------------------------------------------------------------------------------------------
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

# GRIB import selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']

# Projection strings
inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
nwc_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '

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
# Clean
del (importer)
del (incaDictionary_ZS)
del (incaDictionary_TT)
del (incaDictionary_TG)

# --------------------------------------------------------------------------
# LOAD INCA TOPOGRAPHY (this might be a different file format in the future)
inca_topo_grid = np.loadtxt(topoFilename)
inca_topo_grid = inca_topo_grid[::-1, :]  # Reorientation

# ---------------------------------------------------------------------------
# LOAD netCDF file

# Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc['projection'] = nwc_projectionString
metadata_nwc['cartesian_unit'] = 'km'
print('netCDF4 load done!')

# ---------------------------------------------------------------------------
# Reproject
# r_nwc has an extra dimension for members r_nwc[index,i,:,:] let's use the first one for now
inca_reprojected_ZS, inca_metadata_reprojected = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TT, _ = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TG, _ = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
topo_reprojected, _ = reproject_grids(np.array([inca_topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
print('Reprojection done!')

# ---------------------------------------------------------------------------
# Interpolations
inca_interpolation_ZS = grid_interpolation(inca_reprojected_ZS[0], inca_reprojected_ZS[1], timeStep=5)
inca_interpolation_TT = grid_interpolation(inca_reprojected_TT[0], inca_reprojected_TT[1], timeStep=5)
inca_interpolation_TG = grid_interpolation(inca_reprojected_TG[0], inca_reprojected_TG[1], timeStep=5)
print('Interpolation done!')

# -------------------------------------------------------
# Diagnosis
#interp_idx = 0
members_time_index = 0

# tests #
if metadata_nwc['timestamps'][members_time_index].minute == 0:
    interp_idx = -1
else:
    interp_idx = np.where(np.arange(0, 65, 5) == metadata_nwc['timestamps'][members_time_index].minute)[0][0]
print("Precipitation member index ", members_time_index, " timestamp ", metadata_nwc['timestamps'][members_time_index])
print("Match to interpolation GRIB index: ", interp_idx)

# Members Mean matrix.
r_nwc_mean = members_mean_matrix_at(r_nwc[:, members_time_index, :, :])     # Input: Member list of grids (3D matrix)

# calculate precipitation type result
ptype = calculate_precip_type(incaZnow=inca_interpolation_ZS[interp_idx, :, :],
                              incaTemp=inca_interpolation_TT[interp_idx, :, :],
                              incaGroundTemp=inca_interpolation_TG[interp_idx, :, :],
                              precipGrid=r_nwc_mean,
                              topographyGrid=topo_reprojected[0, :, :])


# PLOT
plot_ptype(np.array([ptype]), inca_metadata_reprojected, 0, metadata_nwc['timestamps'][members_time_index], dir_gif)

# filenames = []
# filenames.append(plot_ptype(np.array([ptype]), metadata_nwc, 0, startdate, dir_gif))

# Think how to put the loop logic in a index function ....loop until the last available COMMON time
# len(metadata_nwc['timestamps'])
# startdate
# inca_interpolation_ZS.shape[0]

# message_idx = 1
# calculate interpolations [message_idx-1] [message_idx]
#for t in timestamps[]:
#    if t.minute == 0 and t != timestamps[0]: # to avoid recalculation if starts at timestamp = 00
        # message_idx = message_idx + 1
        # calculate interpolation [message_idx-1] [message_idx]
    # calculate mean value for nwc[:,t,:,:]
    # interp_idx = get minute value using np.arange(0,65,5) <> [t] (value between 0 to 12)
    # precipitation type [interp_idx, meanMatrix, topo ]






