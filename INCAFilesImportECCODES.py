
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


# ---------------------------------------------------------------------------
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
filenames.append(worker(r_reprojected, metadata_reprojected, 0))



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


# ---------------------------------------------------------------------------
# Linear interpolation

# example!
from scipy import interpolate

# x = r_reprojected[0, 200, 200]
# y = r_reprojected[1, 200, 200]
# f = interpolate.interp1d([0, 60], [r_reprojected[0, 200, 200], r_reprojected[1, 200, 200]])
# xnew = np.arange(5, 65, 5)
# for i in xnew:
#     print(f(i))


### This should be a function that receives as parameter the reprojection matrix ###

# Find start index of the non nan value over the reprojection
for i in range(r_reprojected.shape[1]):
    if not all(np.isnan(r_reprojected[0, i, :])):
        x_start = i
        break

for i in range(r_reprojected.shape[2]):
    if not all(np.isnan(r_reprojected[0, :, i])):
        y_start = i
        break

# Create interpolations point by point in r_interpolation matrix
# Interpolation grid list
#r_interpolation = [[[float("nan") for i in range(r_reprojected.shape[1])] for j in range(r_reprojected.shape[2])] for t in range(r_reprojected.shape[0]-1)]
# Interpolation grid 3-D
#r_interpolation = [[[float("nan")] for i in range(r_reprojected.shape[1])] for j in range(r_reprojected.shape[2])]
r_interpolation = [[[float("nan") for k in range(12)] for i in range(r_reprojected.shape[1])] for j in range(r_reprojected.shape[2])]
r_interpolation = np.array(r_interpolation)

# Time profile
# from datetime import datetime
# now1 = datetime.now().time()

# Calculate interpolation cube for R_inca
time_point_range = [int(i) for i in metadata_nwc['leadtimes'] % 60][0:12]
inter_x1 = 0
inter_x2 = 60

t = 1  # slice index
for i in range(x_start, (R_inca.shape[1] + x_start)):
    for j in range(y_start, (R_inca.shape[2] + y_start)):
        interpolationFunction = interpolate.interp1d([inter_x1, inter_x2], [r_reprojected[t - 1, i, j], r_reprojected[t, i, j]])  # how to assign this function
        r_interpolation[i][j] = interpolationFunction(time_point_range)
print('DONE', sep='')

# TEST
# r_reprojected[0, 200, 200]
# r_reprojected[1, 200, 200]
# r_interpolation[0][200][200](np.arange(0, 65, 5))
# r_interpolation[200,200]

# m = np.array([[1,2,3],[4,5,6],[7,8,9]])
# f = np.array([[1,0,0],[0,1,0],[0,0,1]])
# m[f == 1] = 100

# Let's use the time (x) as interpolation parameter
time_point = [int(i) for i in metadata_nwc['leadtimes'] % 60]
#r_interpolation[200][200](time_point[0:12])


for i in range(x_start, (R_inca.shape[1] + x_start)):
    for j in range(y_start, (R_inca.shape[2] + y_start)):
        r_interpolation[i][j]()


if time_point[i] == 0:  # update t = t + 1
    # Recalculate the interpolation function
else:
    #












