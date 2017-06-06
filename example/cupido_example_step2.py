## cupido_example_step2.py
#
#   This is an example script to demonstrate the cupido.py function.
#
#   The user can easily add more configurations to evaluate the functionality
#   of the tool.
#
#
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.1 
#   Created:    19 September 2016
#   Modified:   
#

import sys
sys.path.insert(0, '../create_double_differences');
from cupido import cupido as _cupido;
from scipy.io import savemat as _savemat;

GPS1D = 1;
GPS3D = 3; # Not 2, as 3D enables 1D too.
LEV   = 4;
INSAR = 8; # For future use.


## Example 1
example_configuration = {
    'netcdf_levelling': './example_netcdf_leveling.nc',
    'netcdf_gps':       './example_netcdf_gps.nc',
    'csv_idealization': '../create_double_differences/cupido_example_idealization_precision.csv',
    'region_of_interest':       [-1000, 3000, 17000, 18000],
    'period_of_interest':       ['19700101', '20500101'], # from <= t <= to
    'techniques_of_interest':   GPS3D|LEV,
    'include_flagged_obs': False,
    };

output = _cupido(**example_configuration);
output['Configuration'] = example_configuration;
_savemat('./cupido_example_result1.mat', output);


## Example 2
example_configuration = {
    'netcdf_levelling': './example_netcdf_leveling.nc',
    'netcdf_gps':       './example_netcdf_gps.nc',
    'csv_idealization': '../create_double_differences/cupido_example_idealization_precision.csv',
    'region_of_interest':       [-1000, 3000, 17000, 18000],
    'period_of_interest':       ['19700101', '20000101'], # from <= t <= to
    'techniques_of_interest':   LEV,
    'include_flagged_obs': False,
    };

output = _cupido(**example_configuration);
output['Configuration'] = example_configuration;
_savemat('./cupido_example_result2.mat', output);



## Example 3
example_configuration = {
    'netcdf_levelling': './example_netcdf_leveling.nc',
    'netcdf_gps':       './example_netcdf_gps.nc',
    'csv_idealization': '../create_double_differences/cupido_example_idealization_precision.csv',
    'region_of_interest':        './cupido_example_polygon.csv',
    'period_of_interest':       ['19700101', '20000101'], # from <= t <= to
    'techniques_of_interest':   GPS3D|LEV,
    'include_flagged_obs': False,
    };

output = _cupido(**example_configuration);
output['Configuration'] = example_configuration;
_savemat('./cupido_example_result3.mat', output);


