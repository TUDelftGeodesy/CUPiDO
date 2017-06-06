## cupido.py
#
#   This is the main function for the selection of geodetic measurement from 
#   NetCDF files.
#
#   The script makes a selection of single difference levelling and/or GPS
#   and returns a dataset of double differences.
#
#   Input:
#     - netcdf_levelling:        NetCDF file with levelling data,
#     - netcdf_gps:              NetCDF file with GPS data,
#     - csv_idealization:        CSV file with idealization precision model parameters,
#     - period_of_interest:      the period of interest,
#     - region_of_interest:      the region of interest,
#     - techniques_of_interest:  one or more technique indicators, e.g.,
#                                LEV = levelling,
#                                GPS1D = GPS in Up-direction,
#                                GPS3D = GPS in North-, East-, Up-direction.
#                                EXAMPLE usage: LEV|GPS3D
#
#   Output:
#     - BENCHMARKS:              table of used benchmarks,
#     - DD_TABLE:                table of double differences,
#     - DATES:                   table of used projects/epochsl,
#     - DD_OBS:                  vector of double difference observations,
#     - DD_COV_MX:               covariance matrix of double difference observations.
#
#
#   The function calls various other functions, see
#
#   netcdf_reading.py, csv_idealization_reading.py, construct_sd2dd_transformation.py,
#   construct_sd2dd_transformation_gps.py, construct_st_idealization_covmx.py,
#   construct_t_idealization_covmx.py and cupido_logging.py.
#
#   Furthermore, a number of standard Python2 modules are used.
#
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.1.4
#   Created:    19 September 2016
#   Modified:   v1.1, 29 September 2016, enabling different benchmark classes, selection
#                                        with or without outliers, polygon selection
#               v1.1.1, 7 October 2016, better handling of nan values in region of interest
#                                       check (avoiding warning)
#               v1.1.2, 9 October 2016, bug fix in case SD observations of a single epoch
#                                       are selected (and therefore no DD can be formed).
#               v1.1.3, 4 November 2016, bug fix in indexing.
#               v1.1.4, 14 February 2017, implementation of unique benchmark list
#                                         (for combined levelling/gps data)
#

import numpy as _np;
import time as _time;

from netcdf_reading                     import read_netCdf                        as _read_netCdf;
from csv_idealization_reading           import read_idealization                  as _read_idealization;
from csv_polygon_reading                import read_polygon                       as _read_polygon;
from csv_polygon_reading                import inside_polygon                     as _inside_polygon;
from construct_sd2dd_transformation     import construct_sd2dd_transformation     as _construct_sd2dd_transformation;
from construct_sd2dd_transformation_gps import construct_sd2dd_transformation_gps as _construct_sd2dd_transformation_gps;
from construct_st_idealization_covmx    import construct_st_idealization_covmx    as _construct_st_idealization_covmx;
from construct_t_idealization_covmx     import construct_t_idealization_covmx     as _construct_t_idealization_covmx;
from scipy.linalg                       import block_diag                         as _blkdiag;

from cupido_logging                     import logging                            as _logger;
from datetime                           import date                               as _date;
from datetime                           import datetime                           as _datetime;

# Initiate logging.
log_file = _logger('cupido', 'txt', True);

# Define the constants used for the configuration.
# Options are defined using binary OR, thus using levelling and 1D GPS: GPS1D | LEV.
# Resulting integer value (5) is interpreted using binary logic.
GPS1D = 1;
GPS3D = 3; # Not 2, as 3D enables 1D too.
LEV   = 4;
INSAR = 8; # For future use.

# date to datenum
def period_to_epoch(period):
    return _date.toordinal(_datetime.strptime(period, '%Y%m%d')) + 366 #366 is Matlab - Python offset

# datenum to date
datefnc = _np.vectorize(lambda d: _date.fromordinal(int(d)-366).strftime('%Y%m%d'));

# datenum to decimal year
def _dec_year(date_num):
    date =  _datetime.strptime(date_num, '%Y%m%d');
    def sinceEpoch(date): # returns seconds since epoch
        return _time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = _datetime(year=year, month=1, day=1)
    startOfNextYear = _datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

dec_year = _np.vectorize(_dec_year);

def cupido(netcdf_levelling=None,
             netcdf_gps=None,
             csv_idealization=None,
             period_of_interest=['19640101', '21000101'],
             region_of_interest=[-7e3, 289e3, 300e3, 629e3],
             techniques_of_interest=(GPS3D | LEV),
             include_flagged_obs=False):

    # At least one of the two datasets should be defined.
    if netcdf_levelling is None and netcdf_gps is None:
        raise ValueError('Neither levelling nor GPS dataset provided. No data to process.');

    # Convert epoch of interest to ordinal dates (Matlab's datenum).
    period_of_interest[0] = period_to_epoch(period_of_interest[0]) if isinstance(period_of_interest[0], str) else period_of_interest[0];
    period_of_interest[1] = period_to_epoch(period_of_interest[1]) if isinstance(period_of_interest[1], str) else period_of_interest[1];

    assert period_of_interest[0] <= period_of_interest[1];

    # Check the supplied region of interest.
    if isinstance(region_of_interest, (list, tuple)):
        assert region_of_interest[0] <= region_of_interest[2];
        assert region_of_interest[1] <= region_of_interest[3];

    # Check the 'include_flagged_obs' option.
    assert isinstance(include_flagged_obs, bool), 'Option \'include_flagged_obs\' should be either True or False (boolean).';

    # Write basic logging.
    log_file.write('NetCDF file Leveling: \'{:s}\'.'.format(netcdf_levelling));
    log_file.write('NetCDF file GPS: \'{:s}\'.'.format(netcdf_gps));
    log_file.write('');
    log_file.write('CSV file idealization precision: \'{:s}\'.'.format(csv_idealization));
    log_file.write('');
    log_file.write('Period of interest: {:.1f} - {:.1f} (ordinal date).'.format(*period_of_interest));
    log_file.write('');
    if isinstance(region_of_interest, (list, tuple)):
        log_file.write('Region of interest: SWx: {:.1f}; SWy: {:.1f}; NEx: {:.1f}; NEy: {:.1f}.'.format(*region_of_interest));
    else:
        log_file.write('Region of interest defined in: {:s}.'.format(region_of_interest));
    log_file.write('');
    log_file.write('Techniques of interest: {:s}.'.format(', '.join(
        [k for i, k in enumerate(('GPS1D', 'GPS3D', 'LEV', 'INSAR')) if techniques_of_interest & 2**i])));
    log_file.write('');
    log_file.write('Use flagged observations: {:s}.'.format(str(include_flagged_obs)));

    # Read idealization precision parameters
    idealization_parameters = _read_idealization(csv_idealization);

    # Read the Region of Interest polygon.
    if isinstance(region_of_interest, str):
        region_of_interest = _read_polygon(region_of_interest);

    # The initial steps are done for both GPS and levelling data.
    def read_data(f, sd2dd_fnc=_construct_sd2dd_transformation, direction_flag=[0, 0, 0], axes=(0, 1, 2)):
        # Simple function to ease memory management for Python.
        # Read a netCdf file and select the stations/observations of interest. Return the data.

        # Get data from netCdf.
        d = _read_netCdf(f);

        # BENCHMARKS = [bm_id, x, y, tech];
        temp_BENCHMARKS = _np.asarray([[r['station_name'].strip(), r['x'], r['y'], r['station_class']] for r in d['Stations']], dtype=_np.object_);

        # Select stations within the area of interest.
        if isinstance(region_of_interest, (list, tuple)): # As [SWx, SWy, NEx, NEy]
            # Ignore erros due to the NaN.
            with _np.errstate(invalid='ignore'):
                points_of_interest = _np.where(
                    (d['Stations']['x'] >= region_of_interest[0]) & (d['Stations']['x'] <= region_of_interest[2]) & \
                    (d['Stations']['y'] >= region_of_interest[1]) & (d['Stations']['y'] <= region_of_interest[3]))[0];
        elif isinstance(region_of_interest, _np.ndarray) and region_of_interest.shape[1] == 2: # As polygon of interest.
            points_of_interest = _np.where(_inside_polygon(_np.column_stack([d['Stations']['x'], d['Stations']['y']]),
                                                           region_of_interest))[0];
        else:
            raise ValueError('Could not detect the type of region of interest.');


        # Select observations within the selected period.
        observations_period_of_interest = (period_of_interest[0] <= d['Observations']['epoch']) & \
                                          (d['Observations']['epoch'] < period_of_interest[1]+1);

        # Select only observations that are not 'flagged', i.e. unfit for use.
        if include_flagged_obs is True:
            observations_not_flagged = _np.ones(d['Observations']['sdObs_flag'].shape, dtype=_np.bool);
        else:
            observations_not_flagged = d['Observations']['sdObs_flag'] <= 0;

        # Select single differences within the period of interest,
        #   taking into account the points of interest.
        observations_of_interest = (_np.in1d(d['Observations']['from_index'], points_of_interest) | \
                                    _np.in1d(d['Observations']['to_index'], points_of_interest));

        # For each of the sensitivity groups, execute the sd2dd algorithm.
        # Initiate variables as list.
        Sfinal = list(); STfromto_final = list(); SD_ind_final = list();
        DD_table = list(); SDObs = list();  SDCov = list();
        for s in axes:
            # Select observations based on the points of interest and the sensitivity of interest (looped).
            observations_sensitivity_of_interest = d['Observations']['sensitivity'][:, s] == 1;

            # Select general list of observations (not constrained to the area of interst).
            temp_obs_all_oi = _np.where(observations_not_flagged &
                                        observations_period_of_interest &
                                        observations_sensitivity_of_interest)[0];

            # Select the points within the area of interest from this list.
            temp_obs_oi = _np.where(observations_of_interest[temp_obs_all_oi])[0];

            if len(temp_obs_oi) <= 0:
                # No points found in this sensitivity class.
                continue;

            assert len(temp_obs_all_oi) >= len(temp_obs_oi), \
                    'There should be as least as many observations as observations of interest.';

            # Construct the DD transformation
            temp_SD_all = _np.asarray([[r['from_index'], r['to_index'], r['project_index']] for r in d['Observations']], dtype=_np.int32);

            temp_Sfinal, temp_STfromto_final, temp_SD_ind_final = sd2dd_fnc(BENCHMARKS_all=temp_BENCHMARKS,
                                                                            dates=d['Projects']['project_epoch'],
                                                                            SD_all=temp_SD_all[temp_obs_all_oi],
                                                                            BM_ind=points_of_interest,
                                                                            SD_list_ind=temp_obs_oi);

            if len(temp_Sfinal) > 0:
                # Add to the list of results.
                Sfinal.append(temp_Sfinal); STfromto_final.append(temp_STfromto_final); SD_ind_final.append(temp_SD_ind_final);

                # Select the used observations.
                if direction_flag[s]==0: #levelling
                    usefullSD_ind = temp_obs_all_oi[temp_SD_ind_final];
                else: #gps
                    usefullSD_ind = temp_obs_all_oi[temp_obs_oi[temp_SD_ind_final]];

                # Select the observations and corresponding covariance matrix.
                SDObs.append(d['SDObs'][usefullSD_ind]);
                SDCov.append(d['SDCov'][_np.meshgrid(usefullSD_ind, usefullSD_ind)]);

                # Create DD_table
                temp_DD_table = temp_STfromto_final.copy();
                temp_DD_table[:, 0:4] += 1;
                temp_DD_table = _np.column_stack((temp_DD_table,
                                                  _np.tile(_np.asarray([s==0, s==1, s==2]), (temp_DD_table.shape[0], 1)), # Sensitivity
                                                  _np.tile(direction_flag[s], (temp_DD_table.shape[0], 1)))); #                      

                DD_table.append(temp_DD_table);

        #create output
        Sout = _blkdiag(*Sfinal);
            
        if len(STfromto_final) > 1:
            STfromto_out = _np.vstack(STfromto_final);
        elif len(STfromto_final)==0:
            STfromto_out = STfromto_final;
        else:
            STfromto_out = STfromto_final[0];

        if len(DD_table) > 1:
            DD_table_out = _np.vstack(DD_table);
        elif len(DD_table)==0:
            DD_table_out = DD_table;
        else:
            DD_table_out = DD_table[0];

        if len(SDObs) > 1:
            SDObs_out = _np.hstack(SDObs);
        elif len(SDObs)==0:
            SDObs_out = SDObs;
        else:
            SDObs_out = SDObs[0];
                
        SDcov_out = _blkdiag(*SDCov);
                
        return {'S':              Sout,
                'STfromto':       STfromto_out,
                'DD_table':       DD_table_out,
                'y':              SDObs_out,
                'Q':              SDcov_out,
                'BENCHMARKS':     temp_BENCHMARKS,
                'Dates':          d['Projects']};

    # Read the levelling and/or GPS data, with previously defined function.
    if (netcdf_levelling is not None) and (techniques_of_interest & 4):
        data_levelling = read_data(f=netcdf_levelling,
                                   sd2dd_fnc=_construct_sd2dd_transformation,
                                   direction_flag = [0, 0, 0], #direction indicator
                                   axes=(2,));
    else:
        data_levelling = None;

    if (netcdf_gps is not None) and (techniques_of_interest & 1):
        # Use separte (non-standard) GPS function for the construction of the sd2dd transformation.
        data_gps = read_data(f=netcdf_gps,
                             sd2dd_fnc=_construct_sd2dd_transformation_gps,
                             direction_flag = [1, 2, 3],
                             axes=((0, 1, 2) if techniques_of_interest & 2 else (2,)));
    else:
        data_gps = None;

    # Combine both datasets, or take one as the base.
    if (((data_levelling is not None) + (data_gps is not None)) == 2) and len(data_levelling['y'])>0 and len(data_gps['y'])>0: # Both used, combine both datasets.
        # Combine BENCHMARKS.
        # Although the points exist multiple times, they come from different techniques (different class).
        BENCHMARKS     = _np.vstack((data_gps['BENCHMARKS'],
                                     data_levelling['BENCHMARKS']));
        bm_list = _np.arange(BENCHMARKS.shape[0]);
        
        # Get cross-indices of levelling and GPS (to remove double benchmarks)
        bm_lev_idx = _np.arange(data_levelling['BENCHMARKS'].shape[0])[_np.in1d(data_levelling['BENCHMARKS'][:,0], 
                           data_gps['BENCHMARKS'][:,0])];
        bm_gps_idx = _np.arange(data_gps['BENCHMARKS'].shape[0])[_np.in1d(data_gps['BENCHMARKS'][:,0], 
                           data_levelling['BENCHMARKS'][bm_lev_idx,0])];
        bm_lev_idx += data_gps['BENCHMARKS'].shape[0];
                
        # Update benchmark list
        bm_list[bm_lev_idx] = bm_gps_idx;
        bm_list += 1; # Benchmark indices start at 1              
                
        # Combine S and Q.
        # [ gps 0         ]
        # [ 0   levelling ]
        S_temp = (data_levelling['S'].shape[0] + data_gps['S'].shape[0],
                  data_levelling['S'].shape[1] + data_gps['S'].shape[1]);
        S = _np.zeros(S_temp);
        Q = _np.zeros((S_temp[1], S_temp[1]));
        S_temp = data_gps['S'].shape;
        S[:S_temp[0], :S_temp[1]] = data_gps['S'];
        Q[:S_temp[1], :S_temp[1]] = data_gps['Q'];
        S[S_temp[0]:, S_temp[1]:] = data_levelling['S'];
        Q[S_temp[1]:, S_temp[1]:] = data_levelling['Q'];
        del S_temp;

        # Combine DD_table.
        # [ gps; levelling ]
        # Map BENCHMARK and Project reference to combined table.
        DD_table_temp = data_levelling['DD_table'].copy();
        DD_table_temp[:, 0:2] += data_gps['BENCHMARKS'].shape[0]; # Reindex to combined data.

        # Reindex levelling benchmarks to GPS benchmarks
        DD_table_temp[:, 0] = bm_list[DD_table_temp[:, 0]-1]; # -1 to account for correct indexing
        DD_table_temp[:, 1] = bm_list[DD_table_temp[:, 1]-1]; # -1 to account for correct indexing
        
        DD_table_temp[:, 2:4] += data_gps['Dates'].shape[0];
        DD_table = _np.vstack((data_gps['DD_table'],
                               DD_table_temp));
        del DD_table_temp;

        # Combine Dates(_table).
        # [ gps; levelling ]
        Dates = _np.append(data_gps['Dates'], data_levelling['Dates']);

        # Combine observations (y).
        # [ gps; levelling ]
        y = _np.append(data_gps['y'], data_levelling['y']);

    else:
        # Only a single dataset is used, copy the reference.
        if (data_levelling is not None) and len(data_levelling['y'])>0:
            var_temp = data_levelling;
        elif (data_gps is not None) and len(data_gps['y'])>0:
            var_temp = data_gps;
        else:
            raise ValueError('There is no data selected in the region of interest or period of interest.')


        S           = var_temp['S'];
        Q           = var_temp['Q'];
        DD_table    = var_temp['DD_table'];
        BENCHMARKS  = var_temp['BENCHMARKS'];
        Dates       = var_temp['Dates'];
        y           = var_temp['y'];


    # Convert from SD to DD.
    Ydd   = S.dot(y);
    Qm_dd = S.dot(Q).dot(S.T);
    Ndd = Ydd.shape[0]; # Number of DD
    

    # Search for the synthetic benchmark. This benchmark should be removed before calculating the idealisation covmx.
    syn_bm_id = _np.where(_np.isnan(BENCHMARKS[:, 1].astype(_np.float32)))[0];
    assert len(syn_bm_id) <= 1, 'Only a single SYN_BM should be present in the dataset.';
    temp_DD_table = DD_table[:, :4].astype(_np.int32);
    # Remove SYN_BM only if one is present, the assert will catch situations where more than one SYN_BM is present.
    if (len(syn_bm_id) == 1):
        temp_DD_table[temp_DD_table[:, 0] == syn_bm_id[0] +1, 0] = 0;
        temp_DD_table[temp_DD_table[:, 1] == syn_bm_id[0] +1, 1] = 0;

    # Select DD in Up direction (no idealization precision for horizontal components).
    up_id = _np.min(_np.where(DD_table[:, 6] == 1)[0], axis=0);
    temp_DD_table = temp_DD_table[up_id:,:];    
    
    # Calculate the idealization matrix.
    idealization_parameters_temp = [_np.float64(idealization_parameters['ST']['variance']),
                                    _np.float64(idealization_parameters['ST']['range']),
                                    _np.float64(idealization_parameters['ST']['power'])];
    BENCHMARKS_temp = BENCHMARKS[:, 1:3].astype(_np.float64);
    Dates_dec_year = dec_year(datefnc(Dates['project_epoch']));

    temp_Qid_ST = _construct_st_idealization_covmx(
                                            IdealNoiseParam=idealization_parameters_temp,
                                            DDtable=temp_DD_table,
                                            BMxy=BENCHMARKS_temp,
                                            Dates=Dates_dec_year);
    Qid_ST = _np.zeros((Ndd, Ndd)); # Create full-size matrix (including horizontal components).
    Qid_ST[up_id:, up_id:] = temp_Qid_ST;

    # Strip trailing whitespace from names of benchmark types.
    # As long as they are not forced to be 10 characters of ASCII, this is the most stable method.
    strip_type = _np.vectorize(lambda r: r.rstrip());

    # As a simple check, show the missing classes (except SYN_BM, as it has no parameters).
    benchmark_types = strip_type(BENCHMARKS[:, 3]);
    param_missing = set(benchmark_types) - set(idealization_parameters['T'].keys()) - {'SYN_BM'};
    if len(param_missing) > 0:
        log_file.write('Missing idealisation parameters for the following classes: {:s}.'.format(param_missing));
    del param_missing;

    temp_Qid_T  = _construct_t_idealization_covmx(
                                            IdealNoiseParam=idealization_parameters['T'],
                                            DDtable=temp_DD_table,
                                            BMxy=BENCHMARKS_temp,
                                            Types=benchmark_types,
                                            Dates=Dates_dec_year);
    Qid_T = _np.zeros((Ndd, Ndd));  # Create full-size matrix (including horizontal components).
    Qid_T[up_id:, up_id:] = temp_Qid_T;
    del benchmark_types, strip_type;


    # Combine both idealization matrices.
    Qid = (Qid_ST +Qid_T)*1e-6;
    del Qid_ST, Qid_T, idealization_parameters_temp, BENCHMARKS_temp;

    # Combine all covariance matrices.
    Q_final = Qid +Qm_dd;

    #TECH_FLAG = {0:'LEV',1:'GPSN',2:'GPSE',3:'GPSU'};
    TECH_FLAG = _np.asarray(['LEV','GPSN','GPSE','GPSU']);
    
    # Convert DD_table to string references.
    DD_table_num = DD_table[:, :4].astype(_np.int32)
    DD_table_final = _np.column_stack((BENCHMARKS[DD_table_num[:, 0] -1, 0], # from
                                       BENCHMARKS[DD_table_num[:, 1] -1, 0], # to
                                       Dates['project_name'][DD_table_num[:, 2:4] -1], # epoch_from, epoch_to
                                       DD_table[:,4:7],
                                       TECH_FLAG[DD_table[:, 7].astype(_np.int32)]));

    # Select used stations within the area of interest. [SWx, SWy, NEx, NEy]
    used_benchmarks = DD_table_final[:, 0:2].ravel();
    used_benchmarks_ind = _np.in1d(BENCHMARKS[:, 0], used_benchmarks);
    unique_bm, unique_bm_ind = _np.unique(BENCHMARKS[used_benchmarks_ind, 0], return_index=True);
    BENCHMARKS_final = BENCHMARKS[used_benchmarks_ind][unique_bm_ind, :3];

    # Select used projects/epochs
    used_epochs = DD_table_final[:, 2:4].ravel();
    used_epochs_ind = _np.in1d(Dates['project_name'], used_epochs);
    unique_epochs, unique_epochs_ind = _np.unique(Dates['project_name'][used_epochs_ind], return_index=True);
    DATES_final = _np.column_stack((Dates['project_name'][used_epochs_ind][unique_epochs_ind], datefnc(Dates['project_epoch'][used_epochs_ind][unique_epochs_ind])));
    DATES_final = _np.asarray([[r[0], r[1]] for r in DATES_final], dtype=_np.object_);
    
    # Produce the combined output
    return {'BENCHMARKS':   BENCHMARKS_final,
            'DD_TABLE':     DD_table_final,
            'DATES':        DATES_final,
            'DD_OBS':       Ydd,
            'DD_COV_MX':    Q_final};


if __name__ == '__main__':
    example_configuration = {
        # Locations of the netCDF databases of the GPS and levelling information.
        # If it is either is undefined or None, the dataset will be ignored. If a dateset is provided,
        #   techniques_of_interest determines whether it will be used or (partially) ignored.
        # Setting the data sets the methods of interest.
        'netcdf_levelling':         './leveling.nc',
        'netcdf_gps':               './gps.nc',

        # Location of the CSV file with the idealization precision parameters.
        'csv_idealization':         'cupido_example_idealization_precision.csv',

        # Region of Interest (ROI)
        # Either defined by a rectangular grid or by a file (string) with a polygon definition.
        # A rectangular grid is defined as a list or tuple of the form: [SWx, SWy, NEx, NEy].
        # SWx <= x <= NEx; SWy <= y <= NEy
        # A polygon definition is given in a separate file (reference the file), in the form:
        #    x,   y
        #    123, 456
        #    789, 584
        #    ..., ...
        # The polygon will be closed by the algorithm (ie. last point does not necessarily have
        # to be the same as the starting point).
        # An additional buffer of 1 meter is applied to the polygon, to avoid numerical effects.
        'region_of_interest':       [-1000, 2000, 28000, 24000],
        #'region_of_interest':       'cupido_polygon.csv',

        # Define the time interval of interest.
        # Format: yyyymmdd (as string), or datenum (offset in days, relative to 1 Jan 0 AD).
        'period_of_interest':       ['19640101', '21000101'], # from <= t <= to
    
        # Define the techniques of interest.
        # Available: GPS1D, GPS3D and LEV.
        # Define as binary mask, e.g.: GPS1D | LEV.
        'techniques_of_interest':   LEV|GPS3D,

        # Include observations flagged as outliers/unfit/etc. in the netCdf.
        'include_flagged_obs':      False
        };
        

    output = cupido(**example_configuration);

    # Display output of this run.
    log_file.write('');
    log_file.write('{:=^45s}'.format(' BENCHMARKS '));
    log_file.write('{:^10s} + {:^10s} + {:^10s}'.format('Code', 'X', 'Y'));
    for r in output['BENCHMARKS']:
        log_file.write('{:^10s} | {:+10.1f} | {:+10.1f}'.format(*r));

    log_file.write('');
    log_file.write('{:=^45s}'.format(' DATES '));
    log_file.write('{:^10s} + {:^10s}'.format('Project', 'Epoch'));
    for r in output['DATES']:
        log_file.write('{:^10s} | {:10s}'.format(*r));

    log_file.write('');
    log_file.write('{:=^45s}'.format(' DD_TABLE '));
    log_file.write('{:^10s} + {:^10s} + {:^10s} + {:10s} + {:11s} + {:10s}'.format('From', 'To', 'Prjct Frm', 'Prjct To', 'Sensitivity', 'Class'));
    for r in output['DD_TABLE']:
        log_file.write('{:^10s} | {:^10s} | {:^10s} | {:^10s} | {:.1f} {:.1f} {:.1f} | {:10s}'.format(*r));

    log_file.write('');
    log_file.write('{:=^45s}'.format(' DD_OBS '));
    log_file.write(_np.array_str(output['DD_OBS'], precision=3));

    log_file.write('');
    log_file.write('{:=^45s}'.format(' DD_COV_MX '));
    log_file.write(_np.array_str(output['DD_COV_MX'], precision=2));

