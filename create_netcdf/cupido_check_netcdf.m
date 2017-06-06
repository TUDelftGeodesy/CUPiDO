function cupido_check_netcdf(netcdf_file,globalattributes, ...
               pntname,pntcrd,pntclass, ...
               prjname,prjepoch,prjclass, ...
               obstable,sdobs,sdcov)  
%CUPIDO_CHECK_NETCDF  Check netcdf file with CUPiDO schema, attributes and data.
%   This function checks a netcdf file with the CUPiDO schema
%   The calling syntax is
%
%      cupido_write_netcdf(netcdf_file,globalattributes, ...
%               pntname,pntcrd,pntclass, ...
%               prjname,prjepoch,prjclass, ...
%               obstable,sdobs,sdcov,sdobsflag)  
%
%   The input variables are
%
%      netcdf_file      netcdf file name with extention .nc
%      globalattributes cell array with the global attributes (see below)
%      pntname          cell array with the point name 
%      pntcrd           matrix with coordinates (x,y or lon,lat)
%      pntclass         cell array or string with the point class 
%      prjname          cell array with the project name 
%      prjepoch         mean project epoch [Matlab date number]
%      prjclass         cell array or string with the project class
%      obstable         observation table with index to from_point, to_point and project 
%      sdobs            vector with the observed height difference [m]
%      sdcov            covariance matrix [m]

%   (c) Hans van der Marel, Delft University of Technology, 2016.

%   Created:    12 Apr 2016 by Hans van der Marel


%% Get information about NetCDF file into structure finfo

finfo=ncinfo(netcdf_file);


%% Read netcdf file for verification

% Read point data

pntname2=cellstr(ncread(netcdf_file,'station_name'));
pntcrd2(:,1)=ncread(netcdf_file,'x');
pntcrd2(:,2)=ncread(netcdf_file,'y');
pntclass2=cellstr(ncread(netcdf_file,'station_class'));

% Project data

prjname2=cellstr(ncread(netcdf_file,'project_name'));
prjepoch2=ncread(netcdf_file,'project_epoch')+datenum('1-Jan-1970 00:00:00');
prjclass2=cellstr(ncread(netcdf_file,'project_class'));

% Observations

%  obstable       observation table with index to from_point, to_point and project 
%  epoch          array with epoch (Matlab date number)
%  sdobs          array with the observed height difference [m]
%  sdcov          covariance matrix [m]
%  sdobsflag      integer observation flag (default 0)
%  sensitivity    sensitivity matrix [0-1]
 
stationFromIndex=ncread(netcdf_file,'stationFromIndex');
stationToIndex=ncread(netcdf_file,'stationToIndex');
projectIndex=ncread(netcdf_file,'projectIndex');

sdobs2=ncread(netcdf_file,'sdObs');
sdcov2=ncread(netcdf_file,'sdCov');

epoch2=ncread(netcdf_file,'epoch')+datenum('1-Jan-1970 00:00:00');
sdobsflag2=ncread(netcdf_file,'sdObsFlag');
sensitivity2=ncread(netcdf_file,'sensitivity');

obstable2=[stationFromIndex stationToIndex projectIndex ];

%% Test the data read from netcdf
if all(strncmp(pntname2,pntname,8)) && ...
                all(all(abs(pntcrd2(~isnan(pntcrd2))-pntcrd(~isnan(pntcrd))) < 1e-1)) && ... 
                all(strncmp(pntclass2,pntclass,3))
   fprintf('Station data is ok.\n')
else
   fprintf('Station data is NOT ok!!!!\n')
   pntname2
   pntcrd2
   pntclass2
end

if all(strncmp(prjname2,prjname,8)) && ...
                all(abs(prjepoch2-prjepoch) < 1e-1) && ... 
                all(strncmp(prjclass2,prjclass,3))
   fprintf('Project data is ok.\n')
else
   fprintf('Project data is NOT ok!!!!\n')
   prjname2
   datestr(prjepoch2)
   prjclass2
end

if  all(abs(sdobs2 - sdobs) < 1e-6 ) && ...
                all(all(abs(sdcov2 - sdcov) < 1e-6 )) && ...
                all(all(obstable2 == obstable)) 
   fprintf('Observation data is ok.\n')
else
   fprintf('Observation data is NOT ok!!!!\n')
   obstable2
   sdobs2
   sdcov2
end

end
