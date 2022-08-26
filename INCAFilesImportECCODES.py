# ---------------------------------------------------------------------------
# READ inca files Format, Plot, Reproject, interpolate
# ---------------------------------------------------------------------------

# Import
from IncaGribImport import IncaGribImporter

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

import pysteps
from pysteps.visualization import plot_precip_field, quiver
from pysteps.utils.reprojection import reproject_grids


# Dictionary to matrix function
def inca_dictionary_to_3Dmatrix(incaDict):
    resultMatrix = np.empty(shape=(
        len(incaDict['Messages'].keys()), incaDict['Messages'][1]['Grid'].shape[0],
        incaDict['Messages'][1]['Grid'].shape[1]))
    for i in range(len(resultMatrix)):
        resultMatrix[i, :, :] = incaDict['Messages'][i + 1]['Grid']
    return resultMatrix


# Plot function
def worker(R_seq, metadata, i):
    title = 'INCA ' + startdate.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i)
    fig = plt.figure(figsize=(15, 15))
    fig.add_subplot(1, 1, 1)
    plot_precip_field(R_seq[i, :, :], geodata=metadata, title=str(i))  # , units='mm', ptype='depth')
    plt.suptitle(title)
    plt.tight_layout()
    filename = f'{i}.png'
    #  filenames.append(filename)
    plt.savefig(os.path.join(dir_gif, filename), dpi=72)
    plt.close()
    return filename


# PLOT!
# for i in range(R_inca.shape[0]):
# filenames = []
# filenames.append(worker(R_inca, metadata_inca, 0))


# test = [[float("nan") for i in range(10)] for j in range(10)]
# test = np.array(test)
# test[3:8,2:7 ] = 1
# np.where(~np.isnan(test))[0][0]
# np.where(~np.isnan(test))[-1][0]
# test[3,1]


# Find the start and end indexes of the inca grid in the reprojected grid (to loop over a smaller grid)
def get_inca_indexes(incaGrid, reprojectedGrid):
    # Find start index
    x_start = np.where(~np.isnan(reprojectedGrid))[0][0]
    y_start = np.where(~np.isnan(reprojectedGrid))[-1][0]
    # Calculate end index
    r_x, r_y = reprojectedGrid.shape
    i_x, i_y = incaGrid.shape
    x_end = (r_x - (r_x - (i_x + x_start))) - 1
    y_end = (r_y - (r_y - (i_y + y_start))) - 1

    return x_start, x_end, y_start, y_end


# Linear interpolation function
# IMPROVE
#   -Building the result matrix during interpolation calculation loop!
#   -Dummy interpolation
#   - use function vectorization to apply over each row (fewer loops)

def grid_interpolation(numpyGridStart, numpyGridEnd, timeStep, timeBase=60, applyOver=None):
    """ Time interpolation between 2 2D grids
    numpyGridStart:
        Numpy 2-D grid of start values (INCA)
    numpyGridEnd:
        Numpy 2-D grid of end values (INCA)
    timeStep:
        time steps considered for interpolation
    timeBase:
        Every hour (60 min)
    applyOver:
        Array with sub-indexes to calculate interpolation (inner grid)
    ----

    Return:
        Returns a list of 3D numpy interpolation matrix
    """
    if numpyGridStart.shape != numpyGridEnd.shape:
        raise ValueError("ERROR: Grids have different dimensions")
    if applyOver is None:
        _x1, _x2, _y1, _y2 = [0, numpyGridStart.shape[0], 0, numpyGridStart.shape[1]]
    else:
        _x1, _x2, _y1, _y2 = applyOver

    interPoints = np.arange(0, (timeBase + timeStep), timeStep)
    inter_x1 = 0
    inter_x2 = timeBase
    interpolationGrid = np.zeros((numpyGridStart.shape[0], numpyGridStart.shape[1], len(interPoints)))
    interpolationGrid[:, :, :] = np.nan

    print('Calculating linear interpolation..')
    for i in range(_x1, _x2):
        for j in range(_y1, _y2):
            interpolationFunction = interpolate.interp1d([inter_x1, inter_x2],
                                                         [numpyGridStart[i, j], numpyGridEnd[i, j]])
            interpolationGrid[i, j, :] = interpolationFunction(interPoints)
    print('Done')

    return interpolationGrid


def calculate_precip_type(incaZnow, incaTemp, incaGroundTemp, precipGrid, topographyGrid, DZML=100., TT0=2., TG0=0.,
                          RRMIN=0):
    """Precipitation type algorithm, returns a 2D matrix with categorical values:
    # PT=0  no precip
    # PT=1  rain
    # PT=2  rain/snow mix
    # PT=3  snow
    # PT=4  freezing rain

    incaZnow:
        INCA snow level 2D grid
    incaTemp:
        INCA temperature 2D grid
    incaGroundTemp:
        INCA ground temperature 2D grid
    precipGrid:
        Precipitation (netCDF PYSTEPS) 2D grid
    topographyGrid:
        Topography grid 2D

    returns:
        2D matrix with categorical data
    """

    # Result grid
    result = np.zeros((precipGrid.shape[0], precipGrid.shape[1]))
    topoZSDiffGrid = (incaInterpGrid_ZS - topographyGrid)  # dzs
    precipMask = (precipGrid > RRMIN)

    # SNOW ((dzs<-1.5*DZML) || ( (ZH[i][j] <= 1.5*DZML) && (dzs<=0)))
    snowMask = (topoZSDiffGrid < (-1.5 * DZML)) | ((topographyGrid <= (1.5 * DZML)) & (topoZSDiffGrid <= 0))
    result[snowMask & precipMask] = 3

    # RAIN+SNOW DIAGNOSIS (dzs < 0.5 * DZML) = 2
    rainSnowMask = ~snowMask & (topoZSDiffGrid < (0.5 * DZML))
    result[rainSnowMask & precipMask] = 2

    # RAIN
    rainMask = ~snowMask & ~rainSnowMask
    result[rainMask & precipMask] = 1

    # FREEZING RAIN DIAGNOSIS 4
    # if ((PT[i][j]==1) && ( (tg_<TG0 && TT[i][j]<TT0) || TT[i][j]<TG0))
    freezingMask = rainMask & ((incaGroundTemp < TG0) & (incaTemp < TT0) | (incaTemp < TG0))
    result[freezingMask] = 4

    return result


# -----------------------------------------------------------------------------
# ------------------------------ # START # ------------------------------------
# -----------------------------------------------------------------------------
# Parameters

# Date
startdate = datetime.datetime.strptime('202201081400', "%Y%m%d%H%M")

# File paths
gribPath = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\08\basic'
fn_zs = '_ZS_FC_INCA.grb'
fn_tt = '_TT_FC_INCA.grb'
fn_tg = '_TG_FC_INCA.grb'
netCDF4Filename = r'C:\Users\talen\Desktop\Student JOB\Data\pysteps_blended_nowcast_20220108-09\blended_nowcast_202201081405.nc'
topoFilename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-TOPO\inca_topo.asc'
dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs'

# GRIB selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']
inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

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
from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

# 2. Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc[
    'projection'] = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
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
# Find the origin and end index of the inca grid over the reprojected grid
x1, x2, y1, y2 = get_inca_indexes(R_inca_ZS[0, :, :], inca_reprojected_ZS[0, :, :])

# PLOT TEST
# filenames = []
# filenames.append(worker(inca_reprojected_ZS, inca_metadata_reprojected, 0))
# filenames.append(worker(inca_reprojected_TT, inca_metadata_reprojected, 0))
# filenames.append(worker(inca_reprojected_TG, inca_metadata_reprojected, 0))

# Fill the values with 0
# for i in range(r_reprojected.shape[0]):
#     r_reprojected[i, :, :] = np.nan_to_num(r_reprojected[i, :, :])

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
subIndexes = [x1, (x2 + 1), y1, (y2 + 1)]  # +1 for the range loop index
inca_interpolation_ZS = grid_interpolation(inca_reprojected_ZS[0], inca_reprojected_ZS[1], timeStep=5,
                                           applyOver=subIndexes)
inca_interpolation_TT = grid_interpolation(inca_reprojected_TT[0], inca_reprojected_TT[1], timeStep=5,
                                           applyOver=subIndexes)
inca_interpolation_TG = grid_interpolation(inca_reprojected_TG[0], inca_reprojected_TG[1], timeStep=5,
                                           applyOver=subIndexes)

# Test
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
# time_step_index = [metadata_nwc['timestamps'][i].minute for i in range(len(metadata_nwc['timestamps']))]
# inca_interValue_indexing = np.arange(0, (timeBase + timeStep), 5)  # np.arange(0, (timeBase + timeStep), timeStep)
# inca_indexing = np.arange(0, len(inca_interValue_indexing), 1)
# common_index = np.where(inca_interValue_indexing == time_step_index[0])[0][0]
#
# # Loop over
# for i in range(len(metadata_nwc['timestamps'])):
#     incaGrid = inca_interpolation[0][:, :, (i + common_index)]  # interpolation 0
#     precipGrid = r_nwc[0, i - common_index, :, :]  # CONSIDERING ONLY THE FIRST MEMBER! = 0

# -------------------------------------------------------
# Diagnosis

## HERE find the index of 15:00 in the inca file and netCDF
# metadata_nwc['timestamps'][10]  # 15:00

# common_index = 10 # index of pysteps timestamp that match the interpolation value in INCA
incaInterpGrid_ZS = np.copy(inca_interpolation_ZS[:, :, 0])
incaInterpGrid_TT = np.copy(inca_interpolation_TT[:, :, 0])
incaInterpGrid_TG = np.copy(inca_interpolation_TG[:, :, 0])
precipGrid = np.copy(r_nwc[0, 0, :, :])
topoGrid = np.copy(topo_reprojected[0, :, :])

# plot test
# filenames = []
# filenames.append(worker(np.array([incaInterpGrid_ZS]), metadata_nwc, 0))
# filenames.append(worker(np.array([incaInterpGrid_TT]), metadata_nwc, 0))
# filenames.append(worker(np.array([incaInterpGrid_TG]), metadata_nwc, 0))
# filenames.append(worker(r_nwc[0,:,:,:], metadata_nwc, 0))

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
#filenames = []
#filenames.append(worker(np.array([ptype]), metadata_nwc, 0))

# Interpolation tests
# inca_interpolation[0][200,200,:]
# inca_reprojected[0,200,200]
# inca_reprojected[1,200,200]
# [inca_reprojected[0,200,200] + i/10*(inca_reprojected[1, 200, 200] - inca_reprojected[0, 200, 200]) for i in range(12)]

# incaInterpGrid_ZS[65, 65]
# inca_interpolation[0][60, 60, 0] == R_inca[0, 0, 0]
# topoGrid[60, 60] == inca_topo_grid[0, 0]
# precipGrid[60, 60] > 0
# topoZSDiffGrid[60, 60]

# PLOT 650 X 660
# result[650, 660] = 1000
# precipGrid[651,661] = 1000
