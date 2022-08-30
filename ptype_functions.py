
# Import
from IncaGribImport import IncaGribImporter

import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

from pysteps.visualization import plot_precip_field, quiver
from visualization.precipitationTypeFields import plot_precipType_field

# Dictionary to matrix function
def inca_dictionary_to_3Dmatrix(incaDict):
    resultMatrix = np.empty(shape=(
        len(incaDict['Messages'].keys()), incaDict['Messages'][1]['Grid'].shape[0],
        incaDict['Messages'][1]['Grid'].shape[1]))
    for i in range(len(resultMatrix)):
        resultMatrix[i, :, :] = incaDict['Messages'][i + 1]['Grid']
    return resultMatrix


# Plot function
def worker(R_seq, metadata, i, startdate, dir_gif):
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

# Linear interpolation function
def grid_interpolation(numpyGridStart, numpyGridEnd, timeStep=5, timeBase=60):
    """ Time interpolation between 2 2D grids

    numpyGridStart:
        Numpy 2-D grid of start values (INCA)
    numpyGridEnd:
        Numpy 2-D grid of end values (INCA)
    timeStep:
        Size of the time step for interpolation (every 5, 10, 15.. min)
    timeBase:
        Time period considered in minutes (e.g. over one hour = 60, 2 hours = 120)
    applyOver:
        Array with sub-indexes to calculate interpolation (inner grid)
    ----

    Return:
        Returns a list of 3D numpy interpolation matrix
    """
    if numpyGridStart.shape != numpyGridEnd.shape:
        raise ValueError("ERROR: Grids have different dimensions")

    interPoints = np.arange(0, (timeBase + timeStep), timeStep)
    interpolationGrid = np.zeros((numpyGridStart.shape[0], numpyGridStart.shape[1], len(interPoints)))
    interpolationGrid[:, :, :] = np.nan

    print('Calculating linear interpolation..')
    for i in range(len(interPoints)):
        interpolationGrid[:, :, i] = numpyGridStart + ((numpyGridEnd - numpyGridStart) / interPoints[-1]) * interPoints[i]
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
        2D matrix with categorical data for each type
    """

    # Result grid
    result = np.zeros((precipGrid.shape[0], precipGrid.shape[1]))
    topoZSDiffGrid = (incaZnow - topographyGrid)  # dzs
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
    freezingMask = (result == 1) & (((incaGroundTemp < TG0) & (incaTemp < TT0)) | (incaTemp < TG0))
    result[freezingMask] = 4

    return result

