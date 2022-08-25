
# ----------------------------------# ----------------------------------# ----------------------------------
# plot PT INCA FILE

from IncaGribImport import IncaGribImporter

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

import pysteps
from pysteps.visualization import plot_precip_field, quiver

filename  = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\08\precip\PT_FC_INCA.grb'
startdate = datetime.datetime.strptime('202201081500', "%Y%m%d%H%M")


dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs'

# Keys to extract from the GRIB messages
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees', 'Latin2InDegrees']
projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

# Import GRIB
importer = IncaGribImporter()
incaDictionary = importer.retrieve_grib_data(filename=filename, metadata_keys=keys)

# Transform to a list of grids
R_inca = np.empty(shape=(len(incaDictionary['Messages'].keys()), incaDictionary['Messages'][1]['Grid'].shape[0], incaDictionary['Messages'][1]['Grid'].shape[1]))

for i in range(len(R_inca)):
    R_inca[i,:,:] = incaDictionary['Messages'][i+1]['Grid']

# Fill GRIB metadata
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


# ---------------------------------------------------------------------------
# PLOT!
filenames = []

# 2 Plots
def worker(R_seq, metadata, i):
  title = 'INCA Z0 ' + startdate.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i)
  fig = plt.figure(figsize=(15, 15))
  fig.add_subplot(1, 1, 1 )
  plot_precip_field(R_seq[i, :, :], geodata=metadata, title=str(i)) # , units='mm', ptype='depth')
  plt.suptitle(title)
  plt.tight_layout()
  filename = f'{i} - ' + startdate.strftime("%Y-%m-%d %H%M") + '.png'
  #  filenames.append(filename)
  plt.savefig(os.path.join(dir_gif, filename), dpi=72)
  plt.close()
  return filename



# ---------------------------------------------------------------------------
# LOAD netCDF FILE to reproject
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
r_reprojected, metadata_reprojected = reproject_grids(R_inca, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

#for i in range(R_inca.shape[0]):
filenames.append(worker(r_reprojected, metadata_inca, 0))





# PLOT TOPO # ---------------------------------------------------------------------------

topoFilename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-TOPO\inca_topo.asc'

inca_topo_grid = np.loadtxt(topoFilename)

# Reorientate
inca_topo_grid = np.flip(inca_topo_grid)
for i in range(inca_topo_grid.shape[0]):
  inca_topo_grid[i, :] = np.flip(inca_topo_grid[i, :])

# reproject
topo_reprojected, topo_metadata_reprojected = reproject_grids(np.array([inca_topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

# plot
filenames.append(worker(topo_reprojected, topo_metadata_reprojected, 0))


