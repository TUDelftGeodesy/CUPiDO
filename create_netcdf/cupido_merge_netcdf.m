function cupido_merge_netcdf(netcdf_file_in1,netcdf_file_in2,netcdf_file_out)
%CUPIDO_MERGE_NETCDF  Merge two CUPiDO Netcdf files 
%   CUPIDO_MERGE_NETCDF(NETCDF_FILE_IN1,NETCDF_FILE_IN2,NETCDF_FILE_OUT) 
%   merges the CUPiDO Netcdf files NETCDF_FILE_IN1 and NETCDF_FILE_IN2
%   into a new output Netcdf file NETCDF_FILE_OUT.
%
%   This function takes the global attributes from the first Netcdf file
%   and copies these to the output file. Global attributes from the second
%   Netcdf are ignored, though a comment is added to the output Netcdf 
%   that it is a merged product.
%
%   Example:
%
%      cupido_merge_netcdf('lts2_gps.nc','lts2_gpsbaseline.nc','lts2_gpsmerged.nc');
%
%   See also cupido_write_netcdf.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016. 

% Created:  30 Sep 2016 by Hans van der Marel
% Modified: 10 Oct 2016 by Hans van der Marel
%            - converted script to function


%% Define netcdf file names 

if nargin ~= 3
    error('Incorrect number of input arguments.')
end


%% Read first netcdf file

finfo1=ncinfo(netcdf_file_in1);

pntname1=cellstr(ncread(netcdf_file_in1,'station_name'));
pntcrd1(:,1)=ncread(netcdf_file_in1,'x');
pntcrd1(:,2)=ncread(netcdf_file_in1,'y');
pntclass1=cellstr(ncread(netcdf_file_in1,'station_class'));

prjname1=cellstr(ncread(netcdf_file_in1,'project_name'));
prjepoch1=ncread(netcdf_file_in1,'project_epoch')+datenum('1-Jan-1970 00:00:00');
prjclass1=cellstr(ncread(netcdf_file_in1,'project_class'));

stationFromIndex=ncread(netcdf_file_in1,'stationFromIndex');
stationToIndex=ncread(netcdf_file_in1,'stationToIndex');
projectIndex=ncread(netcdf_file_in1,'projectIndex');

sdobstable1=[stationFromIndex stationToIndex projectIndex ];

sdobs1=ncread(netcdf_file_in1,'sdObs');
sdcov1=ncread(netcdf_file_in1,'sdCov');

epoch1=ncread(netcdf_file_in1,'epoch')+datenum('1-Jan-1970 00:00:00');
sdobsflag1=ncread(netcdf_file_in1,'sdObsFlag');
sensitivity1=ncread(netcdf_file_in1,'sensitivity');


%% Read second netcdf file

finfo2=ncinfo(netcdf_file_in2);

pntname2=cellstr(ncread(netcdf_file_in2,'station_name'));
pntcrd2(:,1)=ncread(netcdf_file_in2,'x');
pntcrd2(:,2)=ncread(netcdf_file_in2,'y');
pntclass2=cellstr(ncread(netcdf_file_in2,'station_class'));

prjname2=cellstr(ncread(netcdf_file_in2,'project_name'));
prjepoch2=ncread(netcdf_file_in2,'project_epoch')+datenum('1-Jan-1970 00:00:00');
prjclass2=cellstr(ncread(netcdf_file_in2,'project_class'));

stationFromIndex=ncread(netcdf_file_in2,'stationFromIndex');
stationToIndex=ncread(netcdf_file_in2,'stationToIndex');
projectIndex=ncread(netcdf_file_in2,'projectIndex');

sdobstable2=[stationFromIndex stationToIndex projectIndex ];

sdobs2=ncread(netcdf_file_in2,'sdObs');
sdcov2=ncread(netcdf_file_in2,'sdCov');

epoch2=ncread(netcdf_file_in2,'epoch')+datenum('1-Jan-1970 00:00:00');
sdobsflag2=ncread(netcdf_file_in2,'sdObsFlag');
sensitivity2=ncread(netcdf_file_in2,'sensitivity');

%% Merge attributes
%
% Take the attributes from the first file and add a comment

finfo=finfo1;

fprintf('Global Attributes:\n\n')
for k=1:numel(finfo.Attributes)
  if strcmpi(finfo.Attributes(k).Name,'comment')
      comment=sprintf('Merged from %s and %s.',netcdf_file_in1,netcdf_file_in2);
      finfo.Attributes(k).Value=[finfo.Attributes(k).Value ' ' comment];
  end
  fprintf('%s: %s\n',finfo.Attributes(k).Name,finfo.Attributes(k).Value)
end
fprintf('\n\n')
globalattributes=[ {finfo.Attributes.Name}' {finfo.Attributes.Value}' ];

%% Merge point data

[pntname,ia,ib]=unique([pntname1;pntname2],'first');

pntcrdtmp=[pntcrd1;pntcrd2];
pntcrd=pntcrdtmp(ia,:);
pntclasstmp=[pntclass1;pntclass2];
pntclass=pntclasstmp(ia);

% Check the coordiantes

if ~isempty(find(abs(pntcrd(ib,:)-pntcrdtmp)>0.5))
    fprintf('Some points have different coordinates\n')
end

% Print point data

numpnt=size(pntname,1);

fprintf('\nBenchmarks (%d points):\n',numpnt)

fprintf('\nPNTNAME             X_RD         Y_RD   CLASS\n\n')
for k=1:numpnt
   fprintf('%-10s  %12.3f %12.3f   %s\n',pntname{k},pntcrd(k,:),pntclass{k})
end
fprintf('\n');

%% Merge project data

if ~isempty(intersect(prjname1,prjname2))
   error('cannot merge two netcdf with identical project names.')
end
prjname=[prjname1;prjname2];
prjepoch=[prjepoch1;prjepoch2];
prjclass=[prjclass1;prjclass2];

%% Merge observation data

stationFromIndex=mkindex([pntname1(sdobstable1(:,1));pntname2(sdobstable2(:,1))],pntname);
stationToIndex=mkindex([pntname1(sdobstable1(:,2));pntname2(sdobstable2(:,2))],pntname);
projectIndex=mkindex([prjname1(sdobstable1(:,3));prjname2(sdobstable2(:,3))],prjname);

sdobstable=[stationFromIndex stationToIndex projectIndex ];

sdobs=[sdobs1;sdobs2];
sdcov=[sdcov1 zeros(size(sdcov1,1),size(sdcov2,2)) ; ...
        zeros(size(sdcov2,1),size(sdcov1,2)) sdcov2  ];
sdobsflag=[sdobsflag1;sdobsflag2];
sensitivity=[sensitivity1;sensitivity2];
epoch=[epoch1;epoch2];


%% write netcdf

cupido_write_netcdf(netcdf_file_out,globalattributes, ...
                pntname,pntcrd,pntclass, ...
                prjname,prjepoch,prjclass, ...
                sdobstable,sdobs,sdcov,sdobsflag,sensitivity,epoch)  

%% Text output

fprintf('\n\nFrom       To         Project       Obs [m] StDev [mm] Flg   N E U  Date\n\n')
for k=1:size(sdobstable,1)
   fprintf('%-10s %-10s %-10s %10.4f %10.2f  %2d  %2d%2d%2d  %s\n',...
       pntname{sdobstable(k,1)},pntname{sdobstable(k,2)},prjname{sdobstable(k,3)},...
       sdobs(k),sqrt(sdcov(k,k))*1000,sdobsflag(k),sensitivity(k,:),datestr(epoch(k),'yyyy-mm-dd HH:MM'))
end

%% Done

end
