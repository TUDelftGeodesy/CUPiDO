%CUPIDO_EXAMPLE_STEP3_MODELING  Script to apply the modeling (inversion) of
%Mogi model parameters based on the results of steps 1 and 2
%   file.
%
%   This script apply the inversion of two parameters of the Mogi model 
%   (i.e. depth and rate of the volume change ), by direct search of a 
%   user-defined search space. The objective function of the inversion is 
%   the L2-norm of the misfit between DD observations and DD predicted by 
%   forward models. The user can change the simulation and search space 
%   settings in the settings section below.
%
%   (c) Sami Samiei Esfahany, Delft University of Technology, 2017.

%   Created:    18 July 2017 by Sami Samiei Esfahany

clc
clear all
close all


data=load('cupido_example_result1.mat');
%data=load('cupido_example_result2.mat')
%data=load('cupido_example_result3.mat')


%% DD observation vector (y)
y = data.DD_OBS'; % in [m]
NumDD=length(y);

%% extract the coordinates of from/to benchmarks for each DD observation

% give coordinate to the SYN_BM
[temp,ind_SYNBM]=ismember(cellstr('SYN_BM'),data.BENCHMARKS(:,1));
if ind_SYNBM
    data.BENCHMARKS(ind_SYNBM,2)=num2cell(single(10000000000));  % x coordinate of the SYN_BM (given outside of the AOI)
    data.BENCHMARKS(ind_SYNBM,3)=num2cell(single(10000000000));  % y coordinate of the SYN_BM (given outside of the AOI)
end

[temp,BMind_from] = ismember(data.DD_TABLE(:,1),data.BENCHMARKS(:,1));
DDx1=cell2mat(data.BENCHMARKS(BMind_from,2));
DDy1=cell2mat(data.BENCHMARKS(BMind_from,3));

[temp,BMind_to] = ismember(data.DD_TABLE(:,2),data.BENCHMARKS(:,1));
DDx2=cell2mat(data.BENCHMARKS(BMind_to,2));
DDy2=cell2mat(data.BENCHMARKS(BMind_to,3));


%% Extract covariance matrix of DD observations
Qdd=data.DD_COV_MX;
invQdd =inv(Qdd);
%% extract the date of from/to epochs for each DD observation
[temp,EPind_from] = ismember(data.DD_TABLE(:,3),data.DATES(:,1));
DDt1=datenum(cell2mat(data.DATES(EPind_from,2)),'yyyymmdd');

[temp,EPind_to] = ismember(data.DD_TABLE(:,4),data.DATES(:,1));
DDt2=datenum(cell2mat(data.DATES(EPind_to,2)),'yyyymmdd');

%% extract sensitivity vector for each DD
DDsens=double(cell2mat(data.DD_TABLE(:,5:7)));

%% extract technique of each DD
DDtech=data.DD_TABLE(:,8);


%% define the search space for the Mogi model
% only two parameters are estimated (depth and the rate of volume change)
% the coordinates of the center of the poit source are assumed to be known
% a-priori (deterministic)

x0=10000; % y coordinate of the point source [m] - fixed (deteministic)
y0=10000; % x coordinate of the point source [m] - fixed (deteministic)

depth_range=[2000 8000]; % range of the search space for the depth [m]
Vrate_range=[-800000 -300000]; % range of the search space for the rate of the volume chage [m^3/year]

depth_int=100;  % interval  of the search space for the depth [m]
Vrate_int=5000; % interval of the search space for the rate of the volume chage [m^3/year]


d=[depth_range(1):depth_int:depth_range(2)]';
v=[Vrate_range(1):Vrate_int:Vrate_range(2)]';
[dd,vv]=meshgrid(d,v);  % serach space

ddall=dd(:);
vvall=vv(:);
NumSearch=length(ddall);

%% loop on forward modeling
% i.e.: for each elemnt of the search space, we evaluate the Mogi model and
% construct the corresponding DD observations. Then we compute the misfit (residuals)
% between the model and DD observations and compute the L2-norm of the
% residuals. The model with minimum L2-norm is selected as the solution.
% Note that the L2-norm is computed in the metric of the inverse of the
% DD covariance matrix.

L2norm=NaN(NumSearch,1);

for i=1:NumSearch
    
    Depth=ddall(i);
    Volume_rate=vvall(i);
    
    y_model_UP= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi * Depth) .* ...
        (((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
        ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2));  % DD Mogi forward model (UP)
    
    y_model_E= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi) .* ...
        ( (DDx2-x0).* ((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
        (DDx1-x0).* ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2) );  %DD Mogi forward model (East)
    
    y_model_N= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi) .* ...
        ( (DDy2-y0).* ((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
        (DDy1-y0).* ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2) );  %DD Mogi forward model (North)
    
    y_model=sum(DDsens .* [y_model_N y_model_E y_model_UP],2);  % construct the final DD model
    e=(y_model-y);  %residuals or midfit
    L2norm(i)=e'*invQdd*e; % L2-norm of the residuals (in the metric of inverse of the DD covarianxce matrix)
    
    progress=100*i/NumSearch
end


%% find the mumimum L2-norm solution
L2norm_grid=reshape(L2norm,size(dd));
[minL2,ind_solution]=min(L2norm_grid(:));

%% Visualization

% the true values for visualization
Depth_True=6000;
Volume_rate_True=-550000;

%% search space
figure
imagesc(d/1000,v,log10(L2norm_grid))
hold on
plot(Depth_True/1000,Volume_rate_True,'ko','markersize',15,'MarkerFaceColor','w')
plot(dd(ind_solution)/1000,vv(ind_solution),'g.','markersize',30,'MarkerFaceColor','g')
temp = get(colorbar,'YTick');
colorbar('Ytick',temp,'YTickLabel',10.^temp);
xlabel('Depth [km]')
ylabel('volume change rate [m^3/year]')
title('L2-norm over the 2D search-space')
legend('True value','L2-norm solution')

%% Covariance matrix of DD observations
figure
imagesc(1e6*Qdd)
title('DD covariance matrix')
k = colorbar;
set(get(k,'title'),'string','mm^2');

%% DD model vs. DD observations

% compute DDs predicted by the solution model
Depth=dd(ind_solution);
Volume_rate=vv(ind_solution);

y_model_UP= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi * Depth) .* ...
    (((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
    ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2));  % DD Mogi forward model (UP)

y_model_E= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi) .* ...
    ( (DDx2-x0).* ((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
    (DDx1-x0).* ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2) );  %DD Mogi forward model (East)

y_model_N= (((DDt2-DDt1)/365.25)*Volume_rate*(1-0.27)/pi) .* ...
    ( (DDy2-y0).* ((DDx2-x0).^2 +(DDy2-y0).^2+Depth^2).^(-3/2) - ...
    (DDy1-y0).* ((DDx1-x0).^2 +(DDy1-y0).^2+Depth^2).^(-3/2) );  %DD Mogi forward model (North)

y_model_solution=sum(DDsens .* [y_model_N y_model_E y_model_UP],2);  % construct the final DD model
e=(y_model_solution-y);  %residuals or midfit
%%
figure
plot(y_model_solution*1000,y*1000,'r.')
xlabel('DD predicted by the L2-norm solution [mm]')
ylabel('DD observations [mm]')
hold on 
plot(abs(max([y_model_solution*1000;y*1000]))*[-1 1],abs(max([y_model_solution*1000;y*1000]))*[-1 1],'k--')
grid on
title('DD model vs. DD observations')
axis equal


