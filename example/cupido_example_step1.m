%CUPIDO_EXAMPLE_STEP1  Script to create an example leveling and GPS CUPiDO netcdf
%   file.
%
%   This script creates two netcdf files with the CUPiDO schema based on
%   a simulated subsidence bowl. The user can change the simulation settings
%   in the settings section below.
%
%   (c) Freek van Leijen, Delft University of Technology, 2017.

%   Created:    26 Mar 2017 by Freek van Leijen



clear all
close all

addpath('../create_netcdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin settings section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
lev.netcdf_file = 'example_netcdf_leveling.nc';
gps.netcdf_file = 'example_netcdf_gps.nc';

lev.start = '19760527';
lev.end = '20161231';
lev.Nepoch = 12;
lev.Nbenchmarks = 200;
lev.std = 0.001; %m/sqrt(km)

gps.start = '19930821';
gps.end = '20161231';
gps.Nepoch = 100;
gps.Nstations = 5;
gps.std = 0.005; %m

aoi.x0 = -5000; %m
aoi.xN = 25000; %m
aoi.y0 = 0; %m
aoi.yN = 20000; %m
aoi.max_defo_rate = -0.003; %m/y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End settings section (You should not have to change anything below this line.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display('Simulating signal ...');

%% Simulate epochs
lev.start_datenum = datenum(lev.start,'yyyymmdd');
lev.end_datenum = datenum(lev.end,'yyyymmdd');
lev.datenums = lev.start_datenum + sort(rand(lev.Nepoch,1))*(lev.end_datenum-lev.start_datenum);
lev.epochs = datestr(lev.datenums,'yyyymmdd');

gps.start_datenum = datenum(gps.start,'yyyymmdd');
gps.end_datenum = datenum(gps.end,'yyyymmdd');
gps.datenums = gps.start_datenum + sort(rand(gps.Nepoch,1))*(gps.end_datenum-gps.start_datenum);
gps.epochs = datestr(gps.datenums,'yyyymmdd');


%% Simulate deformation signal
aoi.x_size = aoi.xN - aoi.x0;
aoi.y_size = aoi.yN - aoi.y0;

[x,y] = meshgrid(aoi.x0:aoi.x_size/500:aoi.xN,aoi.y0:aoi.y_size/500:aoi.yN);
x0 = aoi.x0 + 0.5*aoi.x_size; % center of bowl
y0 = aoi.x0 + 0.5*aoi.x_size; % center of bowl
r = min(aoi.x_size,aoi.y_size)/8; % radius
aoi.bowl = aoi.max_defo_rate*exp(-0.5*((x-x0).^2 + (y-y0).^2)/(r^2));


%% Simulate benchmark locations
lev.x = rand(lev.Nbenchmarks,1)*aoi.x_size+aoi.x0;
lev.y = rand(lev.Nbenchmarks,1)*aoi.y_size+aoi.y0;

unique_stations = 0;
while unique_stations == 0
  gps.station_idx = unique(randi(lev.Nbenchmarks,gps.Nstations,1));
  if length(gps.station_idx) == gps.Nstations
    unique_stations = 1;
  end
end
gps.x = lev.x(gps.station_idx);
gps.y = lev.y(gps.station_idx);


%% Simulate leveling data
display(' ');
display('Simulating leveling data set ...');

for v = 1:lev.Nepoch
  
  epoch.Nbenchmarks = 0;
  while epoch.Nbenchmarks < 3
    epoch.benchmark_idx = unique(ceil(rand(ceil(rand(1)*lev.Nbenchmarks),1)*lev.Nbenchmarks));
    epoch.Nbenchmarks = length(epoch.benchmark_idx);
  end   

  trian = delaunay(lev.x(epoch.benchmark_idx),lev.y(epoch.benchmark_idx));
  % to get a 'nice' triangulation
  trian = sort(trian,2); 
  % operations to extract the unique arcs from the Delaunay triangulation
  A = trian(:,1:2);
  B = trian(:,2:3);
  C = trian(:,1:2:3);
  a = union(A,B,'rows');
  epoch.arcs = union(a,C,'rows');
  epoch.Nobs = size(epoch.arcs,1);

  epoch.x1 = lev.x(epoch.benchmark_idx(epoch.arcs(:,1)));
  epoch.x2 = lev.x(epoch.benchmark_idx(epoch.arcs(:,2)));
  epoch.y1 = lev.y(epoch.benchmark_idx(epoch.arcs(:,1)));
  epoch.y2 = lev.y(epoch.benchmark_idx(epoch.arcs(:,2)));

  hdl = NaN(3,1);
  figure(v);hold on
  imagesc(0.001*x(1,:),0.001*y(:,1),1000*((lev.datenums(v)- lev.start_datenum)/365.25)*aoi.bowl);
  if aoi.max_defo_rate > 0
    caxis([0 1000*((lev.datenums(end)- lev.start_datenum)/365.25)*aoi.max_defo_rate]); 
  elseif aoi.max_defo_rate < 0
    caxis([1000*((lev.datenums(end)- lev.start_datenum)/365.25)*aoi.max_defo_rate 0]); 
  end
  k = colorbar;
  set(get(k,'title'),'string','mm');
  h = plot(0.001*lev.x,0.001*lev.y,'k.');
  hdl(2) = h(1);
  h = plot(0.001*gps.x,0.001*gps.y,'ko','linewidth',2);
  hdl(1) = h(1);
  h = plot(0.001*[epoch.x1 epoch.x2]',0.001*[epoch.y1 epoch.y2]','k');
  hdl(3) = h(1);
  set(gca,'Ydir','normal');
  axis equal
  axis tight
  xlabel('X [km]');
  ylabel('Y [km]');
  title(['Leveling campaign ' lev.epochs(v,:)])
  legend(hdl,'GPS station','leveling benchmark','leveling measurement',...
         'Location','SouthOutside', 'Orientation','horizontal');
  box on;

  epoch.signal = ((lev.datenums(v)- lev.start_datenum)/365.25)*aoi.max_defo_rate*...
             (-exp(-0.5*((epoch.x2-x0).^2 + (epoch.y2-y0).^2)/(r^2)) - ...
             -exp(-0.5*((epoch.x1-x0).^2 + (epoch.y1-y0).^2)/(r^2)));
  epoch.std = lev.std*sqrt(0.001*hypot(epoch.x2-epoch.x1,epoch.y2-epoch.y1));
  epoch.noise = randn(epoch.Nobs,1).*epoch.std;
  epoch.obs = epoch.signal + epoch.noise;

  %% Adjustment
  A = sparse(epoch.Nobs,epoch.Nbenchmarks);
  for w = 1:epoch.Nobs
    A(w,epoch.arcs(w,1)) = -1;
    A(w,epoch.arcs(w,2)) = 1;
  end
  A(:,1) = []; %set reference benchmark
  
  Qy = diag(epoch.std.^2);
  invQy = inv(Qy);
  Qxhat = inv(A'*invQy*A);
  xhat = Qxhat*A'*invQy*epoch.obs;

  lev_epochs(v).Nobs = epoch.Nbenchmarks - 1;
  lev_epochs(v).val = xhat;
  lev_epochs(v).cov = Qxhat;
  lev_epochs(v).table = [repmat(epoch.benchmark_idx(1),lev_epochs(v).Nobs,1) epoch.benchmark_idx(2:end)];
  clear epoch

end


%% Create leveling dataset
lev.Nobs = sum([lev_epochs.Nobs]);
lev.sdcov = zeros(lev.Nobs,lev.Nobs);
lev.sdobs = zeros(lev.Nobs,1);
lev.obstable = NaN(lev.Nobs,3);

count = 0;
for v = 1:lev.Nepoch
  lev.sdobs(count+1:count+lev_epochs(v).Nobs) = lev_epochs(v).val;
  lev.sdcov(count+1:count+lev_epochs(v).Nobs,count+1:count+lev_epochs(v).Nobs) = lev_epochs(v).cov;
  lev.obstable(count+1:count+lev_epochs(v).Nobs,:) = [lev_epochs(v).table repmat(v,lev_epochs(v).Nobs,1)];
  count = count+lev_epochs(v).Nobs;
end

figure;
imagesc(lev.sdcov);
k = colorbar;
set(get(k,'title'),'string','m^2');
title('Covariance matrix leveling single differences');


lev.pntcrd = [lev.x lev.y];
lev.pntclass = repmat({'LEV'},lev.Nbenchmarks,1);
lev.pntname = cellstr([repmat('PNT',lev.Nbenchmarks,1) num2str((1:lev.Nbenchmarks)','%05d')]);

lev.prjname = cellstr([repmat('lv',lev.Nepoch,1) lev.epochs]); %current version, max 10 characters
lev.prjepoch = lev.datenums;
lev.prjclass = repmat('Primary',lev.Nepoch,1);


%% Write leveling netcdf file
display(' ');
display('Writing leveling netcdf file ...');

lev.globalattributes = { ...
      'title'        , 'Example leveling dataset'  ; ...
      'institution'  , 'Delft University of Technology, Netherlands.' ; ...
      'source'       , 'Simulation' ; ...
      'technique'    , 'Leveling' ; ...
      'history'      , ' ' ; ...
      'references'   , 'TU Delft, 2016.' ; ...
      'comment'      , ' ' ; ...
      'Conventions'  , 'CF-1.6' ; ...
      'featureType'  , 'timeSeries' ; ...
      'email'        , 'f.j.vanleijen@tudelft.nl' ; ...
      'version'      , '1.0' ; ...
      'termsForUse'  , 'These data can be used freely for research purposes.' ; ...
      'disclaimer'   , 'This data is made available in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' ; ...
};

cupido_write_netcdf(lev.netcdf_file,lev.globalattributes, ...
                lev.pntname,lev.pntcrd,lev.pntclass, ...
                lev.prjname,lev.prjepoch,lev.prjclass, ...
                lev.obstable,lev.sdobs,lev.sdcov,0,[0 0 1]);

display(' ');
display('Checking leveling netcdf file ...');
cupido_check_netcdf(lev.netcdf_file,lev.globalattributes, ...
                lev.pntname,lev.pntcrd,lev.pntclass, ...
                lev.prjname,lev.prjepoch,lev.prjclass, ...
                lev.obstable,lev.sdobs,lev.sdcov);


%% Simulate GPS data
display(' ');
display('Simulating leveling data set ...');

for v = 1:gps.Nepoch

  epoch.signal = ((gps.datenums(v)- lev.start_datenum)/365.25)*aoi.max_defo_rate*...
             (-exp(-0.5*((gps.x-x0).^2 + (gps.y-y0).^2)/(r^2)));

  gps_epochs(v).Nobs = gps.Nstations;
  gps_epochs(v).up = epoch.signal + randn(gps.Nstations,1)*gps.std;
  gps_epochs(v).north = randn(gps.Nstations,1)*gps.std*0.5;
  gps_epochs(v).east = randn(gps.Nstations,1)*gps.std*0.5;
  
  gps_epochs(v).table = [repmat(gps.Nstations+1,gps.Nstations,1) (1:gps.Nstations)'];
 
end


%% Create GPS dataset
gps.Nobs = gps.Nstations*gps.Nepoch;
gps.sdcov = zeros(3*gps.Nobs,3*gps.Nobs);
gps.sdobs = zeros(3*gps.Nobs,1);
gps.obstable = NaN(3*gps.Nobs,3);

count = 0;
for v = 1:gps.Nepoch
  gps.sdobs(count+1:count+gps_epochs(v).Nobs) = gps_epochs(v).up;
  gps.sdobs(gps.Nobs+count+1:gps.Nobs+count+gps_epochs(v).Nobs) = gps_epochs(v).north;
  gps.sdobs(2*gps.Nobs+count+1:2*gps.Nobs+count+gps_epochs(v).Nobs) = gps_epochs(v).east;
  
  gps.obstable(count+1:count+gps_epochs(v).Nobs,:) = [gps_epochs(v).table repmat(v,gps_epochs(v).Nobs,1)];
  gps.obstable(gps.Nobs+count+1:gps.Nobs+count+gps_epochs(v).Nobs,:) = [gps_epochs(v).table repmat(v,gps_epochs(v).Nobs,1)];
  gps.obstable(2*gps.Nobs+count+1:count+2*gps.Nobs+gps_epochs(v).Nobs,:) = [gps_epochs(v).table repmat(v,gps_epochs(v).Nobs,1)];

  count = count+gps_epochs(v).Nobs;
end

gps.sdcov = diag([repmat(gps.std^2,1,gps.Nobs) repmat((0.5*gps.std)^2,1,2*gps.Nobs)]);
gps.sens = [repmat([0 0 1],gps.Nobs,1); repmat([1 0 0],gps.Nobs,1); repmat([0 1 0],gps.Nobs,1)];

figure;
imagesc(gps.sdcov);
k = colorbar;
set(get(k,'title'),'string','m^2');
title('Covariance matrix GPS single differences');


gps.pntcrd = [[gps.x gps.y];[NaN NaN]];
gps.pntclass = [repmat({'CORS'},gps.Nstations,1);{'SYN_BM'}];
gps.pntname = [cellstr([repmat('PNT',gps.Nstations,1) num2str(gps.station_idx,'%05d')]);{'SYN_BM'}];

gps.prjname = cellstr([repmat('cg',gps.Nepoch,1) gps.epochs]); %current version, max 10 characters
gps.prjepoch = gps.datenums;
gps.prjclass = repmat('Continuous',gps.Nepoch,1);


%% Write GPS netcdf file
display(' ');
display('Writing GPS netcdf file ...');

gps.globalattributes = { ...
      'title'        , 'Example GPS dataset'  ; ...
      'institution'  , 'Delft University of Technology, Netherlands.' ; ...
      'source'       , 'Simulation' ; ...
      'technique'    , 'Continuous GPS' ; ...
      'history'      , ' ' ; ...
      'references'   , 'TU Delft, 2016.' ; ...
      'comment'      , ' ' ; ...
      'Conventions'  , 'CF-1.6' ; ...
      'featureType'  , 'timeSeries' ; ...
      'email'        , 'f.j.vanleijen@tudelft.nl' ; ...
      'version'      , '1.0' ; ...
      'termsForUse'  , 'These data can be used freely for research purposes.' ; ...
      'disclaimer'   , 'This data is made available in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' ; ...
};

cupido_write_netcdf(gps.netcdf_file,gps.globalattributes, ...
                gps.pntname,gps.pntcrd,gps.pntclass, ...
                gps.prjname,gps.prjepoch,gps.prjclass, ...
                gps.obstable,gps.sdobs,gps.sdcov,0,gps.sens);

display(' ');
display('Checking GPS netcdf file ...');
cupido_check_netcdf(gps.netcdf_file,gps.globalattributes, ...
                gps.pntname,gps.pntcrd,gps.pntclass, ...
                gps.prjname,gps.prjepoch,gps.prjclass, ...
                gps.obstable,gps.sdobs,gps.sdcov);

display(' ');
display('Done.');
%% Done


