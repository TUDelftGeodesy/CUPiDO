function cupido_write_netcdf(netcdf_file,globalattributes, ...
               pntname,pntcrd,pntclass, ...
               prjname,prjepoch,prjclass, ...
               obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch)  
%CUPIDO_WRITE_NETCDF  Create netcdf file with CUPiDO schema, attributes and data.
%   This function creates a netcdf file with the CUPiDO schema, writes 
%   global attributes to the netcdf file, and fills the file with data.   
%   The calling syntax is
%
%      cupido_write_netcdf(netcdf_file,globalattributes, ...
%               pntname,pntcrd,pntclass, ...
%               prjname,prjepoch,prjclass, ...
%               obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch)  
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
%      sdobsflag        integer observation flag (default 0)
%      sensitivity      sensitivity matrix [0-1] (default [ 0 0 1 ])
%      epoch            array with epoch [Matlab date number] (default prjepoch)
%
%   The Global attitributes (NETCDF CF-1.6) must be defined outside this 
%   function according to the following rules:
%
%   globalattributes = { ...
%     'title'        , 'A concise description of what is in the dataset.'  ; ...
%     'institution'  , 'Specifies where the original data was produced.' ; ...
%     'source'       , 'The method of production of the original data. If it was model-generated, source should name the model and its version, as specifically as could be useful. If it is observational, source should characterize it, e.g., surface observation or radiosonde.' ; ...
%     'history'      , 'Provides an audit trail for modifications to the original data. Well-behaved generic netCDF filters will automatically append their name and the parameters with which they were invoked to the global history attribute of an input netCDF file. We recommend that each line begin with a timestamp indicating the date and time of day that the program was executed.' ; ...
%     'references'   , 'Published or web-based references that describe the data or methods used to produce it.' ; ...
%     'comment'      , 'Miscellaneous information about the data or methods used to produce it.' ; ...
%     'conventions'  , 'CF-1.6' ; ...
%     'email'        , 'Email address of the contact person.' ; ...
%     'version'      , ' ' ; ...
%     'termsForUse'  , 'Specifies the terms of use, if any.' ; ...
%     'featureType'  , 'Specifies the type of discrete sampling geometry to which the data in the file belongs, and implies that all data variables in the file contain collections of features of that type.' ; ... 
%   };
%
%   Example:
%     
%     netcdf_file='gps.nc';
%     globalattributes = { ...
%       'title'        , 'GPS observations in the Netherlands.'  ; ...
%       'institution'  , 'Delft University of Technology, Netherlands.' ; ...
%       'source'       , 'GPS height database.' ; ...
%       'technique'    , 'GPS' ; ...
%       'history'      , ' ' ; ...
%       'references'   , 'TU Delft, 2016.' ; ...
%       'comment'      , ' ' ; ...
%       'conventions'  , 'CF-1.6' ; ...
%       'featureType'  , 'timeSeries' ; ...
%       'email'        , 'h.vandermarel@tudelft.nl' ; ...
%       'version'      , ' ' ; ...
%       'termsForUse'  , 'These data can be used freely for research purposes.' ; ...
%       'disclaimer'   , 'This data is made available in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' ; ...
%     };
%     cupido_write_netcdf(netcdf_file,globalattributes, ...
%               pntname,pntcrd,pntclass, ...
%               prjname,prjepoch,prjclass, ...
%               obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch)  
%
%   The last two arguments are optional and pntclass and prjclass can
%   be string instead of string arrays, and for sdobsflag and sensitivity
%   default values can be given for all observations. For example 
%
%     cupido_write_netcdf(netcdf_file,globalattributes, ...
%               pntname,pntcrd,'LEV', ...
%               prjname,prjepoch,'LEV', ...
%               obstable,sdobs,sdcov,0,[ 0 0 1])  
%
%   You can use ncdisp(netcdf_file) at all time to inspect the netcdf
%   file.
%
%   See also NCDISP.

%   (c) Hans van der Marel, Delft University of Technology, 2016.

%   Created:    12 Apr 2016 by Hans van der Marel
%   Modified:   12 Apr 2016 by Hans van der Marel
%                - Initial dummy version
%               23 August 2016 by Hans van der Marel
%                - Changes to the netcdf structure
%                - Converted script into functions
%                - Use real variables
%               25 August 2016 by Hans van der Marel
%                - corrected few typo's in documentation

%%  Default length of character strings in netcdf

NAME_STRLEN=10;         

%% Check the input variables and provide defaults

if nargin < 11 || nargin > 14
    error('Incorrect number of input arguments')
end
if nargin < 12
   sdobsflag=0;
end
if nargin < 13
   sensitivity=[0 0 1];
end
if nargin < 14
   epoch=[];
end

%% Set the dimensions

numpnt=size(pntname,1);
numprj=size(prjname,1);
numobs=size(obstable,1);

%% Convert input variables to proper format

if ischar(pntclass) && size(pntclass,1) == 1
    tmp=pntclass;
    pntclass={};
   [pntclass{1:numpnt,1}]=deal(tmp);
end
if ischar(prjclass) && size(prjclass,1) == 1
    tmp=prjclass;
    prjclass={};
   [prjclass{1:numprj,1}]=deal(tmp);
end
if isscalar(sdobsflag)
   sdobsflag=repmat(sdobsflag,[numobs,1]);
end
if size(sensitivity,1) == 1 
   sensitivity=repmat(sensitivity,[numobs,1]);
end
if isempty(epoch)
   epoch=prjepoch(obstable(:,3));
end

%% Create netcdf schema file
%
% NETCDF Schema: Indexed ragged array representation of time series 
%
% This part needs only some elementary dimensions, 
%
% numpnt         number of points (scalar)
% numprj         number of projects (scalar)
% numobs         number of observations (scalar) 
%
% and global attributes defined above.

fprintf('Create CUPiDO netcdf schema ...\n')

% NETCDF schema

mySchema.Name   = '/';
mySchema.Format = 'classic';

% Global attributes

for k=1:size(globalattributes,1)
   mySchema.Attributes(k).Name = globalattributes{k,1};
   mySchema.Attributes(k).Value = globalattributes{k,2};
end

% Dimensions

mySchema.Dimensions(1).Name   = 'station';
mySchema.Dimensions(1).Length = numpnt;

mySchema.Dimensions(2).Name   = 'project';
mySchema.Dimensions(2).Length = numprj;

mySchema.Dimensions(3).Name   = 'obs';
mySchema.Dimensions(3).Length = numobs;  % Inf

mySchema.Dimensions(4).Name   = 'name_strlen';
mySchema.Dimensions(4).Length = NAME_STRLEN;  

mySchema.Dimensions(5).Name   = 'components';
mySchema.Dimensions(5).Length = 3;  

% Variables

kvar=0;

% Station variables

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'station_name';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([1 4]);
mySchema.Variables(kvar).Datatype    = 'char';
variableattributes = { ...
  'long_name'     , 'benchmark name' ; ...
  'standard_name' , 'station name' ; ...  
  'cf_role'       , 'timeseries_id' ; ...      
  'comment'       , 'There is one synthetic benchmark with id SYN_BM used for GPS, with XY-coordinate of inf. This benchmark is used as an imaginary reference benchmark for GPS data, and its location is assumed to be far from the subsidence  AOI (or in the edge of the simulation area) with zero deformation.' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'x';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(1);
mySchema.Variables(kvar).Datatype    = 'single';
variableattributes = { ...
  'long_name'     , 'x-coordinate'
  'standard_name' , 'projection_x_coordinate' ; ...  
  'units'         , 'm' ; ...      
  'axis'          , 'X' ; ...
  'comment'       , 'e.g., specify reference system' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'y';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(1);
mySchema.Variables(kvar).Datatype    = 'single';
variableattributes = { ...
  'long_name'     , 'y-coordinate'
  'standard_name' , 'projection_y_coordinate' ; ...  
  'units'         , 'm' ; ...      
  'axis'          , 'Y' ; ...
  'comment'       , 'e.g., specify reference system' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'station_class';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([1 4]);
mySchema.Variables(kvar).Datatype    = 'char';
variableattributes = { ...
  'long_name'     , 'Benchmark class' ; ...
  'standard_name' , 'station class' ; ...     
  'comment'       , 'Possible values are [SYN_BM|GPS|CORS|LEV|...]' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

% Project/epoch variables

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'project_name';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([2 4]);
mySchema.Variables(kvar).Datatype    = 'char';
variableattributes = { ...
  'long_name'     , 'Project name' ; ...
  'standard_name' , 'project name' ; ...  
  'comment'       , 'Name of NAM project as used in the NAM database.' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'project_epoch';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(2);
mySchema.Variables(kvar).Datatype    = 'double';
variableattributes = { ...
  'long_name'     , 'mean epoch of the project' ; ...
  'standard_name' , 'mean epoch' ;  ...
  'units'         , 'days since 1970-01-01 00:00:00' ; ...
  'comment'       , 'Mean epoch of measurement for a project. Actual measurement times are often not accurately known, most of the time this is the estimated middle time of the project. For the purpose of forming single differences in time all measurements in a project/epoch are considered to be taken at the same time.' ; ...  
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'project_class';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([2 4]);
mySchema.Variables(kvar).Datatype    = 'char';
variableattributes = { ...
  'long_name'     , 'Project class' ; ...
  'standard_name' , 'project class' ; ...  
  'comment'       , 'Possible values are [GPS|LEV[#]|...] with # the levelling order (optional).' ; ...
};
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

% Observation variables

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'stationFromIndex';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'int32';
variableattributes = { ...
  'long_name'     , 'index of the station this observation is from (reference)' ; ...
  'instance_dimension' , 'station' ; ...
  'comment'       , 'Note: GPS observations refer to an imaginary reference benchmark (SYN_BM) which location is assumed to be far from the subsidence AOI (or in the edge of the simulation area) with zero deformation. ' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'stationToIndex';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'int32';
variableattributes = { ...
  'long_name'     , 'index of the station this observation is to' ; ...
  'instance_dimension' , 'station' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end


kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'projectIndex';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'int32';
variableattributes = { ...
  'long_name'     , 'index of the project/epoch this observation is for' ; ...
  'instance_dimension' , 'project' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'epoch';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'double';
variableattributes = { ...
  'long_name'     , 'measurement epoch' ; ...
  'standard_name' , 'time' ;  ...
  'units'         , 'days since 1970-01-01 00:00:00' ; ...
  'comment'       , 'Time of measurement. The time of measurement may not be accurately known, most of the time this is the estimated middle time of the project, which is stored in the variable project_epoch. This field is used to store the exact measurement time (if available), or else it is a copy of the project_epoch.' ; ...  
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'sdObs';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'double';
variableattributes = { ...
  'long_name'     , 'single (spatial) difference observation' ; ...
  'units'         , 'm' ; ...      
  'positive'      , 'up' ; ...
  'coordinates'   , 'epoch x y' ; ...
  'comment'       , 'Observed height difference between two benchmarks. For GPS actually the observed ellipsoidal height is used as GPS uses a stable virtual reference benchmark outside the area of interest.' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'sdCov';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([3 3]);
mySchema.Variables(kvar).Datatype    = 'double';
variableattributes = { ...
  'long_name'     , 'single difference covariance matrix' ; ...
  'units'         , 'm^2' ; ...      
  'comment'       , 'Covariance matrix for the single difference observations, including idealization noise for benchmarks, but without spatial/temporal idealization noise for shallow compaction.' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'sdObsFlag';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions(3);
mySchema.Variables(kvar).Datatype    = 'int32';
variableattributes = { ...
  'long_name'     , 'Observation flag' ; ...
  'comment'       , 'Observation flag: 0 is good to use, 1 or higher indicates outliers (do not use).' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end

kvar=kvar+1;
mySchema.Variables(kvar).Name        = 'sensitivity';
mySchema.Variables(kvar).Dimensions  = mySchema.Dimensions([3 5]);
mySchema.Variables(kvar).Datatype    = 'single';
variableattributes = { ...
  'long_name'     , 'sensitivity matrix' ; ...
  'units'         , 'ratio' ; ...      
  'comment'       , 'Matrix with sensitivity on the measurements to North, East and Up components of deformation. For example [0 0 1] for levelling and GPS_up, [1 0 0] for GPS_north and [0 1 0] for GPS_east.' ; ...
  };
for k=1:size(variableattributes,1)
   mySchema.Variables(kvar).Attributes(k).Name = variableattributes{k,1};
   mySchema.Variables(kvar).Attributes(k).Value = variableattributes{k,2};
end


%% Delete old files

if exist(netcdf_file,'file')
   fprintf('Netcdf file %s already exists, will be deleted first to recreate it.\n',netcdf_file)
   % If you want to automatically update the history uncomment the
   % following lines.
   % history = ncreadatt(netcdf_file,'/','history');   
   % for k=1:size(globalattributes,1)
   %    if strcmp(mySchema.Attributes(k).Name,'history')
   %      mySchema.Attributes(k).Value = [ history ' ' datestr(now) ' ' mySchema.Attributes(k).Value ';'];
   %   end
   % end
   delete(netcdf_file);
end

%% Write schema to file

fprintf('Write CUPiDO netcdf schema to file...\n')
ncwriteschema(netcdf_file, mySchema);

%ncdisp(netcdf_file);

%% Write data to netcdf file

fprintf('Write data to CUPiDO netcdf...\n')

% Station data

%  pntname        cell array with the point name 
%  pntcrd         matrix with x- and y-coordinate
%  pntclass       cell array with the point class 

ncwrite(netcdf_file,'station_name',char(pntname))
ncwrite(netcdf_file,'x',single(pntcrd(:,1)))
ncwrite(netcdf_file,'y',single(pntcrd(:,2)))
%ncwrite(netcdf_file,'lat',single(lat))
%ncwrite(netcdf_file,'lon',single(lon))
%ncwrite(netcdf_file,'height',single(pntcrd(:,3)))
ncwrite(netcdf_file,'station_class',char(pntclass))

% Project data

%  prjname        cell array with the project name 
%  prjepoch       mean project epoch (Matlab date number)
%  prjclass       cell array with the project class

ncwrite(netcdf_file,'project_name',char(prjname))
ncwrite(netcdf_file,'project_epoch',double(prjepoch-datenum('1-Jan-1970 00:00:00')))
ncwrite(netcdf_file,'project_class',char(prjclass))

% Observations

%  obstable       observation table with index to from_point, to_point and project 
%  epoch          array with epoch (Matlab date number)
%  sdobs          array with the observed height difference [m]
%  sdcov          covariance matrix [m]
%  sdobsflag      integer observation flag (default 0)
%  sensitivity    sensitivity matrix [0-1]
 
ncwrite(netcdf_file,'stationFromIndex',int32(obstable(:,1)));
ncwrite(netcdf_file,'stationToIndex',int32(obstable(:,2)));
ncwrite(netcdf_file,'projectIndex',int32(obstable(:,3)));
ncwrite(netcdf_file,'epoch',double(epoch-datenum('1-Jan-1970 00:00:00')));
ncwrite(netcdf_file,'sdObs',double(sdobs));
ncwrite(netcdf_file,'sdCov',double(sdcov));
ncwrite(netcdf_file,'sdObsFlag',int32(sdobsflag));
ncwrite(netcdf_file,'sensitivity',single(sensitivity));

end

