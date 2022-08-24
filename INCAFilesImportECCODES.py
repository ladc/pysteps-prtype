
# ---------------------------------------------------------------------------
# READ inca files Format, Plot, Reproject, interpolate

from IncaGribImport import IncaGribImporter

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

import pysteps
from pysteps.visualization import plot_precip_field, quiver

filename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\08\basic\202201080000_ZS_FC_INCA.grb'

startdate = datetime.datetime.strptime('202201080000', "%Y%m%d%H%M")
dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs'
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees', 'Latin2InDegrees']
projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

importer = IncaGribImporter()
incaDictionary = importer.retrieve_grib_data(filename=filename, metadata_keys=keys)

# Transform to a cube
R_inca = np.empty(shape=(len(incaDictionary['Messages'].keys()), incaDictionary['Messages'][1]['Grid'].shape[0], incaDictionary['Messages'][1]['Grid'].shape[1]))

for i in range(len(R_inca)):
    R_inca[i,:,:] = incaDictionary['Messages'][i+1]['Grid']

metadata_inca = {}
metadata_inca['projection'] = projectionString
metadata_inca['xpixelsize'] = incaDictionary['Metadata']['DxInMetres']
metadata_inca['ypixelsize'] = incaDictionary['Metadata']['DyInMetres']
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
#metadata_inca['yorigin'] = 'lower'
metadata_inca['yorigin'] = 'upper'

#del(incaDictionary)


# ---------------------------------------------------------------------------
# PLOT
from pysteps.utils.reprojection import reproject_grids

filenames = []

# normalize between 0-1 and plot probability
minVal = np.min(R_inca)
maxVal = np.max(R_inca)
diff = np.max(R_inca) - np.min(R_inca)
a = np.copy(R_inca[0, :, :])

for i in range(a.shape[0]):
    a[i, :] = (a[i, :] - minVal) / diff

a_r_reprojected, metadata_reprojected = reproject_grids(np.array([a]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
filenames.append(worker(a_r_reprojected, metadata_reprojected, 0))


# PLOT!
filenames = []

# TEST TEST
def worker(R_seq, metadata, i):
  title = 'INCA Z0 ' + startdate.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i)
  fig = plt.figure(figsize=(15, 15))
  fig.add_subplot(1, 1, 1 )
  plot_precip_field(R_seq[i, :, :], geodata=metadata, ptype="prob",  probthr=0.0, title=str(i)) # , units='mm', ptype='depth')
  plt.suptitle(title)
  plt.tight_layout()
  filename = f'{i}.png'
  #  filenames.append(filename)
  plt.savefig(os.path.join(dir_gif, filename), dpi=72)
  plt.close()
  return filename


###
# 2 Plots
def worker(R_seq, metadata, i):
  title = 'INCA Z0 ' + startdate.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i)
  fig = plt.figure(figsize=(15, 15))
  fig.add_subplot(1, 1, 1 )
  plot_precip_field(R_seq[i, :, :], geodata=metadata, title=str(i)) # , units='mm', ptype='depth')
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
# LOAD netCDF FILE!!
from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

netCDF4Filename = r'C:\Users\talen\Desktop\Student JOB\Data\pysteps_blended_nowcast_20220108-09\blended_nowcast_202201081405.nc'

# 2. Set the filename and load file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)
metadata_nwc['projection'] = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
metadata_nwc['cartesian_unit'] = 'km'
print('done!')


# --------------------------------------------------------------------------- maybe do this after the interpolation!
# Reproject
# r_nwc has a extra dimension for members r_nwc[index,i,:,:]
# let's use the first one for now
r_reprojected, metadata_reprojected = reproject_grids(R_inca, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

# Fill the values with 0
# perhaps do this after the interpolation to make it faster (smaller grid)
# for i in range(r_reprojected.shape[0]):
#     r_reprojected[i, :, :] = np.nan_to_num(r_reprojected[i, :, :])

# PLOT the reprojected thing! next!
# for i in range(r_reprojected.shape[0]):
#filenames.append(worker(r_reprojected, metadata_reprojected, 0))




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


# Test
# interpolationResult = grid_interpolation(r_reprojected[0:2], timeStep=5)
# interpolationResult[0][210, 210, 0] == r_reprojected[0, 210, 210]
# interpolationResult[0][210, 210, 12] == r_reprojected[1, 210, 210]


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
inca_interValue_indexing = np.arange(0, 65, 5) # np.arange(0, (timeBase + timeStep), timeStep)
inca_indexing = np.arange(0, len(inca_interValue_indexing), 1)
common_index = np.where(inca_interValue_indexing==time_step_index[0])[0][0]

# Loop over
for i in range(len(metadata_nwc['timestamps'])):
    incaGrid = interpolationResult[0][:,:,i + common_index]     # interpolation 0
    precipGrid = r_nwc[0, i-common_index, :, :]                 # CONSIDERING ONLY THE FIRST MEMBER! = 0






if time_step_index[10] == 0:
    # re-calculate the interpolations from inca to the next hour
else:
    interpolationResult[0].[:,:, time_step_index[0]]

time_point = [int(i) for i in metadata_nwc['leadtimes'] % 60]
#r_interpolation[200][200](time_point[0:12])














