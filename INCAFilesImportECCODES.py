
# ---------------------------------------------------------------------------
# READ inca files Format, Plot, Reproject, interpolate

from IncaGribImport import IncaGribImporter

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

import pysteps
from pysteps.visualization import plot_precip_field, quiver
from pysteps.utils.reprojection import reproject_grids

startdate = datetime.datetime.strptime('202201081400', "%Y%m%d%H%M")

path = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\08\basic'
filename_ZS = os.path.join(path, startdate.strftime("%Y%m%d%H%M") + '_ZS_FC_INCA.grb')

dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs'
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees', 'Latin2InDegrees']
projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

importer = IncaGribImporter()
incaDictionary_ZS = importer.retrieve_grib_data(filename=filename_ZS, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_ZS = np.empty(shape=(len(incaDictionary_ZS['Messages'].keys()), incaDictionary_ZS['Messages'][1]['Grid'].shape[0], incaDictionary_ZS['Messages'][1]['Grid'].shape[1]))
for i in range(len(R_inca_ZS)):
    R_inca_ZS[i, :, :] = incaDictionary_ZS['Messages'][i+1]['Grid']

del(importer)

metadata_inca = {}
metadata_inca['projection'] = projectionString
metadata_inca['xpixelsize'] = incaDictionary_ZS['Metadata']['DxInMetres']
metadata_inca['ypixelsize'] = incaDictionary_ZS['Metadata']['DyInMetres']
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
metadata_inca['yorigin'] = 'upper'


#---------------------------------------------------------------------------
# Import TT
filename_TT = os.path.join(path, startdate.strftime("%Y%m%d%H%M") + '_TT_FC_INCA.grb')

importer = IncaGribImporter()
incaDictionary_TT = importer.retrieve_grib_data(filename=filename_TT, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TT = np.empty(shape=(len(incaDictionary_TT['Messages'].keys()), incaDictionary_TT['Messages'][1]['Grid'].shape[0], incaDictionary_TT['Messages'][1]['Grid'].shape[1]))
for i in range(len(R_inca_TT)):
    R_inca_TT[i, :, :] = incaDictionary_TT['Messages'][i + 1]['Grid'] - 273.15  # to celsius

del(importer)

#---------------------------------------------------------------------------
# Import TG
filename_TG = os.path.join(path, startdate.strftime("%Y%m%d%H%M") + '_TG_FC_INCA.grb')
importer = IncaGribImporter()
incaDictionary_TG = importer.retrieve_grib_data(filename=filename_TG, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TG = np.empty(shape=(len(incaDictionary_TG['Messages'].keys()), incaDictionary_TG['Messages'][1]['Grid'].shape[0], incaDictionary_TG['Messages'][1]['Grid'].shape[1]))
for i in range(len(R_inca_TG)):
    R_inca_TG[i, :, :] = incaDictionary_TG['Messages'][i + 1]['Grid'] - 273.15  # to celsius

# Clean
del(importer)
del(incaDictionary_ZS)
del(incaDictionary_TT)
del(incaDictionary_TG)



# ---------------------------------------------------------------------------
# LOAD netCDF FILE
from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

netCDF4Filename = r'C:\Users\talen\Desktop\Student JOB\Data\pysteps_blended_nowcast_20220108-09\blended_nowcast_202201081405.nc'

# 2. Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc['projection'] = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
metadata_nwc['cartesian_unit'] = 'km'
print('done!')



# ---------------------------------------------------------------------------
# Reproject
# r_nwc has an extra dimension for members r_nwc[index,i,:,:]
# let's use the first one for now
inca_reprojected_ZS, inca_metadata_reprojected = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TT, inca_metadata_reprojected = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
inca_reprojected_TG, inca_metadata_reprojected = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

# Fill the values with 0
# perhaps do this after the interpolation to make it faster (smaller grid)
# for i in range(r_reprojected.shape[0]):
#     r_reprojected[i, :, :] = np.nan_to_num(r_reprojected[i, :, :])

# PLOT the reprojected thing! next!
# for i in range(r_reprojected.shape[0]):
#filenames.append(worker(r_reprojected, metadata_reprojected, 0))



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

# PLOT!
filenames = []

# 2 Plots
def worker(R_seq, metadata, i):
  title = 'INCA ' + startdate.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i)
  fig = plt.figure(figsize=(15, 15))
  fig.add_subplot(1, 1, 1 )
  plot_precip_field(R_seq[i, :, :], geodata=metadata, title=str(i))     # , units='mm', ptype='depth')
  plt.suptitle(title)
  plt.tight_layout()
  filename = f'{i}.png'
  #  filenames.append(filename)
  plt.savefig(os.path.join(dir_gif, filename), dpi=72)
  plt.close()
  return filename

# for i in range(R_inca.shape[0]):
#filenames.append(worker(R_inca, metadata_inca, 0))





# ---------------------------------------------------------------------------
# Linear interpolation function

from scipy import interpolate

# MAYBE consider to return only 1 grid

# IMPROVE by building the result matrix during interpolation calculation loop!
# SEPARATE in two different functions (utils):
#   -Dummy interpolation
#   -change the parameters to receive the x_start, y_start.. as parameters (default = 0), end = just len()
#       (def grid_interpolation(numpyGridList, timeStep, timeBase=60, row_extent=[60, 651], col_extent=[60, 661] (0,0) (0,0)):
#       maybe the indexes can be a parameter also like = [0, 5, 10, .... 60].. or maybe not, depends on the index selection (if starts at 10?)

def grid_interpolation(numpyGridList, timeStep, timeBase=60, findReprojectedData=True, incaDimensions=[591, 601]):
    """
    numpyGridList:
        Numpy list of 2-D grids for each hour (INCA)
    timeStep:
        time steps considered for interpolation
    timeBase:
        Every hour (60 min)
    findReprojectedData:
        Find the index of the 0,0 corner of the original INCA file in the reprojected grid (smaller grid)
    incaDimensions:
        Dimensions of the INCA grid
    ----

    Return:
        Returns a list of 3D numpy interpolation matrix
    """

    # Find start index of the non nan value over the reprojection grid (to loop over the smaller GRIB grid only)
    # get indexes FUNCTION
    if findReprojectedData:
        for i in range(numpyGridList.shape[1]):
            if not all(np.isnan(numpyGridList[0, i, :])):
                x_start = i
                break
        for i in range(numpyGridList.shape[2]):
            if not all(np.isnan(numpyGridList[0, :, i])):
                y_start = i
                break

    # Calculate end index
    x_end = numpyGridList.shape[1] - (numpyGridList.shape[1] - (incaDimensions[0] + x_start))
    y_end = numpyGridList.shape[2] - (numpyGridList.shape[2] - (incaDimensions[1] + y_start))
    #

    interpolationPoints = int(timeBase / timeStep)
    time_point_range = np.arange(0, (timeBase+timeStep), timeStep)
    inter_x1 = 0
    inter_x2 = timeBase
    gridTemplate = [[[float("nan") for k in range(interpolationPoints + 1)] for i in range(numpyGridList.shape[1])] for j
                       in range(numpyGridList.shape[2])]
    gridTemplate = np.array(gridTemplate)
    interpolationList = []

    for g in range(1, numpyGridList.shape[0]):
        print('Calculating linear interpolation for index', g-1, ' to ', g)
        interpolationGrid = np.copy(gridTemplate)
        for i in range(x_start, x_end):
            for j in range(y_start, y_end):
                interpolationFunction = interpolate.interp1d([inter_x1, inter_x2], [numpyGridList[g-1][i, j], numpyGridList[g][i, j]])
                interpolationGrid[i, j] = interpolationFunction(time_point_range)
        interpolationList.append(interpolationGrid)
    print('Done')

    return interpolationList


inca_interpolation_ZS = grid_interpolation(inca_reprojected_ZS[0:2], timeStep=5)
inca_interpolation_TT = grid_interpolation(inca_reprojected_TT[0:2], timeStep=5)
inca_interpolation_TG = grid_interpolation(inca_reprojected_TG[0:2], timeStep=5)

# Test
# inca_interpolation[0][210, 210, 0] == r_reprojected[0, 210, 210]
# inca_interpolation[0][210, 210, 12] == r_reprojected[1, 210, 210]


# --------------------------------------------------------------------------
# LOAD INCA TOPO
topoFilename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-TOPO\inca_topo.asc'

inca_topo_grid = np.loadtxt(topoFilename)

# Reorientate
inca_topo_grid = np.flip(inca_topo_grid)
for i in range(inca_topo_grid.shape[0]):
  inca_topo_grid[i, :] = np.flip(inca_topo_grid[i, :])

# Reproject
topo_reprojected, topo_metadata_reprojected = reproject_grids(np.array([inca_topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

# PLOT
#filenames.append(worker(topo_reprojected, topo_metadata_reprojected, 0))



# --------------------------------------------------------------------------
# calculate the interpolation points (index selection between inca and pysteps) using the timestamps as index from the pysteps

# How to determine the starting point , we should only use a index seq and determine the right grid by timestamp in metadata_nwc

timeBase = 60
timeStep = 5

time_step_index = [metadata_nwc['timestamps'][i].minute for i in range(len(metadata_nwc['timestamps']))]
inca_interValue_indexing = np.arange(0, (timeBase + timeStep), 5) # np.arange(0, (timeBase + timeStep), timeStep)
inca_indexing = np.arange(0, len(inca_interValue_indexing), 1)
common_index = np.where(inca_interValue_indexing==time_step_index[0])[0][0]

# Loop over
for i in range(len(metadata_nwc['timestamps'])):
    incaGrid = inca_interpolation[0][:,:, (i + common_index)]     # interpolation 0
    precipGrid = r_nwc[0, i-common_index, :, :]                 # CONSIDERING ONLY THE FIRST MEMBER! = 0



# -------------------------------------------------------
# Diagnosis

# Let's just try the filtering process #
# PT=0  no precip
# PT=1  rain
# PT=2  rain/snow mix
# PT=3  snow
# PT=4  freezing rain


## HERE find the index of 15:00 in the inca file and netCDF
metadata_nwc['timestamps'][10]  # 15:00


#common_index = 10 # index of pysteps timestamp that match the interpolation value in INCA
incaInterpGrid_ZS = np.copy(inca_interpolation_ZS[0][:, :, 0])
incaInterpGrid_TT = np.copy(inca_interpolation_TT[0][:, :, 0])
incaInterpGrid_TG = np.copy(inca_interpolation_TG[0][:, :, 0])
precipGrid = np.copy(r_nwc[0, 0, :, :])
topoGrid = np.copy(topo_reprojected[0, :, :])

# Result grid
result = np.zeros((precipGrid.shape[0], precipGrid.shape[1]))

DZML = 100.
TT0 = 2.
TG0 = 0.
RRMIN = 0

topoZSDiffGrid = (incaInterpGrid_ZS - topoGrid)     # dzs

# incaInterpGrid_ZS[65, 65]
# inca_interpolation[0][60, 60, 0] == R_inca[0, 0, 0]
# topoGrid[60, 60] == inca_topo_grid[0, 0]
# precipGrid[60, 60] > 0
# topoZSDiffGrid[60, 60]

# RAIN
result[precipGrid > RRMIN] = 1


# RAIN+SNOW DIAGNOSIS
# (dzs < 0.5 * DZML) = 2
result[(topoZSDiffGrid < (0.5 * DZML)) & (precipGrid > RRMIN)] = 2


# SNOW
# # ( (dzs<-1.5*DZML) || ( (ZH[i][j] <= 1.5*DZML) && (dzs<=0)) )
result[(topoZSDiffGrid < (-1.5 * DZML)) & (precipGrid > RRMIN)] = 3
result[(topoGrid <= (1.5 * DZML)) & (topoZSDiffGrid <= 0) & (precipGrid > RRMIN)] = 3

# PLOT 650 X 660
# result[650, 660] = 1000
# precipGrid[651,661] = 1000

# FREEZING RAIN DIAGNOSIS 4
# if ((PT[i][j]==1) && ( (tg_<TG0 && TT[i][j]<TT0) || TT[i][j]<TG0) )
result[(result == 1) & ((incaInterpGrid_TG < TG0) & (incaInterpGrid_TT < TT0) | (incaInterpGrid_TT < TG0))] = 4

# PLOT
filenames.append(worker(np.array([result]), inca_metadata_reprojected, 0))

# Interpolation tests
# inca_interpolation[0][200,200,:]
# inca_reprojected[0,200,200]
# inca_reprojected[1,200,200]
#
# [inca_reprojected[0,200,200] + i/10*(inca_reprojected[1, 200, 200] - inca_reprojected[0, 200, 200]) for i in range(12)]









