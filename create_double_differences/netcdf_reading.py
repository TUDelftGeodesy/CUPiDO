## netcdf_reading.py
#
#   This function reads a NetCDF file with geodetic data (single differences).
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    19 September 2016
#   Modified:   
#

import scipy.io as _sp_io;
import numpy as _np;

def read_netCdf(f):
    with _sp_io.netcdf_file(f, 'r') as a:

        # Point data
        PointData = _np.zeros(a.variables['x'].shape[0], dtype=[
                                ('station_name', 'a10'), ('x', 'f4'), ('y', 'f4'), ('station_class', 'a10')]);
        PointData['station_name']  = _np.asfortranarray(a.variables['station_name'].data).view('a10');
        PointData['x']             = a.variables['x'].data;
        PointData['y']             = a.variables['y'].data;
        PointData['station_class'] = _np.asfortranarray(a.variables['station_class'].data).view('a10');
        PointData = PointData.copy(); # Copy data to memory (from netCDF mmap).


        # Project data
        ProjectData = _np.zeros(a.variables['project_epoch'].shape[0], dtype=[
                                ('project_name', 'a10'), ('project_epoch', 'f8'), ('project_class', 'a10')]);
        ProjectData['project_name']  = _np.asfortranarray(a.variables['project_name'].data).view('a10');
        ProjectData['project_epoch'] = a.variables['project_epoch'].data + 719529; # Plus 1970-01-01 00:00:00
        ProjectData['project_class'] = _np.asfortranarray(a.variables['project_class'].data).view('a10');
        ProjectData = ProjectData.copy();  # Copy data to memory (from netCDF mmap).

        # Observations
        Observations = _np.zeros(a.variables['epoch'].shape[0], dtype=[
                                 ('from_index', 'uint32'), ('to_index', 'uint32'), ('project_index', 'uint32'),
                                 ('epoch', 'f8'), ('sdObs_flag', 'b'), ('sensitivity', 'f4', 3)]);
        Observations['from_index']    = a.variables['stationFromIndex'].data -1; # -1 for Matlab to Python index
        Observations['to_index']      = a.variables['stationToIndex'].data -1; # -1 for Matlab to Python index
        Observations['project_index'] = a.variables['projectIndex'].data -1; # -1 for Matlab to Python index
        Observations['epoch']         = a.variables['epoch'].data + 719529;
        Observations['sdObs_flag']    = a.variables['sdObsFlag'].data;
        Observations['sensitivity']   = a.variables['sensitivity'].data.T;
        Observations = Observations.copy(); # Copy data to memory (from netCDF mmap).

        # Matrices
        SDObs = a.variables['sdObs'].data.copy();
        SDCov = a.variables['sdCov'].data.copy();

    return {'Stations': PointData,
            'Projects': ProjectData,
            'Observations': Observations,
            'SDObs': SDObs,
            'SDCov': SDCov};

if __name__ == '__main__':
    d = read_netCdf('leveling.nc');

    for k, v in d.iteritems():
        print '{:=^50s}'.format(' '+ str(k) +' ' +str(v.shape));
        if v.dtype.names is not None:
            print ', '.join(v.dtype.names);

        print v;
