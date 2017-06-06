function data = cupido_read_netcdf2struct(netcdf_file)
%% Read LTS2 Netcdf file to struct variable
%
% *Sami Samiei Esfahany, Delft University of Technology, 5 September 2016 *
%
% This Matlab function reads an LTS2 Netcdf file and output the content
% in a struct variable.
%

% (c) Sami Samiei Esfahany, Delft University of Technology, 2016.

% Created:    12 October 2016 by Sami Samiei Esfahany
% Modified:   
%

%CUPIDO_READ_NETCDF2STRUCT  Read CUPiDO Netcdf file to struct variable
%   CUPIDO_READ_NETCDF2STRUCT(NETCDF_FILE) reads the CUPiDO Netcdf file NETCDF_FILE
%   and provides a structure array as output.
%
%   data = cupido_read_netcdf2struct(netcdf_file)
%   reads the CUPiDO Netcdf file NETCDF_FILE and provides a structure array as output.
%
%   Example:
%
%      data = cupido_read_netcdf2struct('cupido_gps.nc');
%
%   See also cupido_write_netcdf and cupido_merge_netcdf.
%

% (c) Sami Samiei Esfahany, Delft University of Technology, 2016.

% Created:    12 October 2016 by Sami Samiei Esfahany
% Modified:   
%

data.PointData.station_name = ncread(netcdf_file,'station_name');
data.PointData.x=ncread(netcdf_file,'x');
data.PointData.y=ncread(netcdf_file,'y');
data.PointData.station_class=ncread(netcdf_file,'station_class');
  
data.ProjectData.project_name=ncread(netcdf_file,'project_name');
data.ProjectData.project_epoch=ncread(netcdf_file,'project_epoch')+datenum('1-Jan-1970 00:00:00');
data.ProjectData.project_class=ncread(netcdf_file,'project_class');

data.Observations.from_index=ncread(netcdf_file,'stationFromIndex');
data.Observations.to_index=ncread(netcdf_file,'stationToIndex');
data.Observations.project_index=ncread(netcdf_file,'projectIndex');

data.Observations.epoch=ncread(netcdf_file,'epoch')+datenum('1-Jan-1970 00:00:00'); 
data.Observations.sdObs_flag=ncread(netcdf_file,'sdObsFlag');
data.Observations.sensitivity=ncread(netcdf_file,'sensitivity');

data.SDObs=ncread(netcdf_file,'sdObs');
data.SDCov=ncread(netcdf_file,'sdCov');


