function varargout = cupido_read_netcdf(netcdf_file)
%CUPIDO_READ_NETCDF  Read CUPiDO Netcdf file 
%   CUPIDO_READ_NETCDF(NETCDF_FILE) reads the CUPiDO Netcdf file NETCDF_FILE
%   and prints a summary of the content.
%
%   [pntname,pntcrd,pntclass, prjname,prjepoch,prjclass,obstable,sdobs,...
%   sdcov,sdobsflag,sensitivity,epoch, finfo] = cupido_read_netcdf(netcdf_file)
%   reads the CUPiDO Netcdf file NETCDF_FILE and outputs the content
%   as variables.
%
%   Example:
%
%      cupido_read_netcdf('cupido_gps.nc');
%
%      [pntname,pntcrd,pntclass, prjname,prjepoch,prjclass, ...
%              obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch, ...
%              finfo] = cupido_read_netcdf('cupido_gps.nc');  
%
%   See also cupido_write_netcdf and cupido_merge_netcdf.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016. 

% Created:  12 Apr 2016 by Hans van der Marel
% Modified: 12 Apr 2016 by Hans van der Marel
%              - Initial version
%           24 Aug 2016 by Hans van der Marel
%              - New structure
%           14 Sep 2016 by Hans van der Marel
%              - Additional plots 
%           10 Oct 2016 by Hans van der Marel
%              - converted script to function

%% Define netcdf file name 

if nargin ~=1, error('This function requires one input argument.');, end


%% Display NetCDF file schema in command window

% ncdisp(netcdf_file)

%% Get information about NetCDF file into structure finfo and print attributes

fprintf('Netcdf file: %s\n\n',netcdf_file)

finfo=ncinfo(netcdf_file);

fprintf('Global Attributes:\n\n')
for k=1:numel(finfo.Attributes)
  fprintf('%s: %s\n',finfo.Attributes(k).Name,finfo.Attributes(k).Value)
end
fprintf('\n\n')

%% Read netcdf file


% Read point data

pntname=ncread(netcdf_file,'station_name');
pntcrd(:,1)=ncread(netcdf_file,'x');
pntcrd(:,2)=ncread(netcdf_file,'y');
pntclass=ncread(netcdf_file,'station_class');

% Project data

prjname=ncread(netcdf_file,'project_name');
prjepoch=ncread(netcdf_file,'project_epoch')+datenum('1-Jan-1970 00:00:00');
prjclass=ncread(netcdf_file,'project_class');

% Observations

%  sdobstable     observation table with index to from_point, to_point and project 
%  epoch          array with epoch (Matlab date number)
%  sdobs          array with the observed height difference [m]
%  sdcov          covariance matrix [m]
%  sdobsflag      integer observation flag (default 0)
%  sensitivity    sensitivity matrix [0-1]
 
stationFromIndex=ncread(netcdf_file,'stationFromIndex');
stationToIndex=ncread(netcdf_file,'stationToIndex');
projectIndex=ncread(netcdf_file,'projectIndex');

sdobs=ncread(netcdf_file,'sdObs');
sdcov=ncread(netcdf_file,'sdCov');

epoch=ncread(netcdf_file,'epoch')+datenum('1-Jan-1970 00:00:00');
sdobsflag=ncread(netcdf_file,'sdObsFlag');
sensitivity=ncread(netcdf_file,'sensitivity');

sdobstable=[stationFromIndex stationToIndex projectIndex ];

%% 

if nargout >=1 
  varargout={  pntname,pntcrd,pntclass, ...
               prjname,prjepoch,prjclass, ...
               sdobstable,sdobs,sdcov,sdobsflag,sensitivity,epoch,finfo }; 
  return;
end

%% Print point data

numpnt=size(pntname,1);

fprintf('\nBenchmarks (%d points):\n',numpnt)

fprintf('\nPNTNAME             X_RD         Y_RD   CLASS\n\n')
for k=1:numpnt
   fprintf('%-10s  %12.3f %12.3f   %s\n',pntname(k,:),pntcrd(k,:),pntclass(k,:))
end
fprintf('\n');

%% Print project data

numprj=size(prjname,1);

fprintf('\nProjects (%d projects):\n',numprj)

fprintf('\nPRJNAME      MEAN_EPOCH   CLASS\n\n')
dfmt='yyyy-mm-dd';
for k=1:numprj
   fprintf('%-10s   %s   %s\n',prjname(k,:), ...
       datestr(prjepoch(k),dfmt),prjclass(k,:));
end
fprintf('\n');

%% Print observation data

numobs=size(sdobstable,1);

fprintf('\nObservations (%d observations):\n',numobs)

fprintf('\nFROM       TO         PROJECT      OBS [m] STDEV [mm]  FLAG   SENSITIVITY   EPOCH\n\n')
dfmt='yyyy-mm-dd HH:MM';
for k=1:numobs
   fprintf('%-10s %-10s %-10s%10.3f %10.3f  %4d   %3.1f %3.1f %3.1f   %s\n', ...
       pntname(sdobstable(k,1),:),pntname(sdobstable(k,2),:),prjname(sdobstable(k,3),:), ...
       sdobs(k),sqrt(sdcov(k,k))*1000,sdobsflag(k),sensitivity(k,:),datestr(epoch(k),dfmt));
end
fprintf('\n');


%% Done

end
