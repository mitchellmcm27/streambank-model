function eval_model_monthly(input_data,row,fitted_mdl)

addpath('utils')
addpath('bdv model')
addpath('gis')
shp_dir = fullfile('data','NHD','lines_Project_Join_new.shp');

% convert lat-long to utm
input_lat = input_data.lat;
input_lon = input_data.long;
[input_x, input_y] = deg2utm(input_lat,input_lon); % returns col vecs
pathrow=input_data.pathrow; % landsat path/row

A = input_data.A;
S = input_data.S;

% determines which regional curve to use
region = input_data.region;
[B, H, Q, U] = model_geometry(A,region);
A_bkf = B.*H;
Rh = A_bkf./(H + B + H);

%roughness and friction
manning_n = 0.087;
Cf0 = 9.81*manning_n.^2./(Rh.^(1/6));

cl_shp = read_nhd_shapefile(input_data.COMID(row),'COMID',shp_dir);
dpx = cl_shp.X(1:end-1);
dpy = cl_shp.Y(1:end-1);

% resample centerline, calculate curvature, normal vector (n), downstream
% coordinate (s), and bank points
[x,y,iR,n,s,bank_points] = model_curvature(dpx, dpy, B(row), input_x(row), input_y(row));

% clip to speed up modeling
% x=x(100:800);
% y=y(100:800);
% n=n(100:800);
% s=s(100:800);
% bank_points=bank_points(100:800,:);


%% 1 D model

% runoff/baseflow volume 
V_est = input_data.noah_runoff_2{row} + input_data.noah_baseflow_2{row}; % kg/m2;
V_est(V_est<0) = 0;
V_est = V_est.*A(row).*1000.*1000; % kg = L = dm3
V_est = V_est./1000; % m3

f=input_data.f_monthly_2{row};
Q_dates = input_data.dates_2{row}; % valid dates for observations at this site
dates = Q_dates;
ndays = eomday(year(Q_dates),month(Q_dates)); % number of days in each month
wet_days = f.*ndays; % number of wet days

% Streamflow
% Assuming events last 1 day and each event is the same magnitude
Q_est = V_est./wet_days./24./60./60; % estimated event Q m3/s
delta_t = wet_days(1:numel(dates))*24*60*60; % assuming 1 event = 1 day, this is the total length of events in seconds

% Normal depth
H_est = model_normal_depth(Q_est,H(row),B(row),S(row),manning_n);

% Hydraulic radius at this depth
Rh_est = H_est.*B(row)./(2.*H_est + B(row));

% Estimate channel roughness assuming constant Manning's n
Cf_est = 9.81*manning_n.^2./(Rh_est.^(1/6)); % Cf = gn^2/Rh^(1/6)

%% Flow model
disp('Monthly flow model...')
U_nb_full = nan(numel(dates),numel(x));
H_nb_full = nan(size(U_nb_full));
yrs = ceil(numel(dates)/12);

figure
plot(x,y,'b')
axis equal
drawnow
% for yr = 1:yrs
%     strt = (yr-1)*12 + 1 % 1, 13, 25...
%     ed = min(numel(dates),yr*12)           % 12, 24, 36,.. inclusive
%     sim_dates = dates(strt:ed)
sim_dates=dates;
for mnth = 1:numel(sim_dates) % only loop through dates where erosion data is available
    sim_dates(mnth)
    if mnth==1
        AR0 = 0;
        % simulate bankfull conditions to set bed topography
        [unl,~,~,~] = ...
            model_velocity_nonlinear_swe(x,y,iR,Q(row),B(row),H(row),s,Cf0(row),S(row),1,AR0,1);
        AR0 = unl.AR;
    end
    
    % Allow flow Q and depth to vary each month,
    % keep bed topography from the bankfull condition
    [unl,~,~,n] = ...
        model_velocity_nonlinear_swe(x,y,iR,Q_est(mnth),B(row),H_est(mnth),s,Cf_est(mnth),S(row),1,AR0,3);
    
    U_nb_full(mnth,:) = unl.q./unl.h.*unl.asR.*B(row)/2;
    H_nb_full(mnth,:) = H(row).*(unl.AR+unl.Fr2R).*B(row)/2;

end
sumU = nanmean(U_nb_full,1)';
bheight = nanmean(H_nb_full,1)' + H(row);

%% Sample data
disp('Sampling data...')
disp('size of vectors')
size(x)
size(n)
size(sumU)

% bank points where sumU is > 0
bank_points = [x - B(row)*sign(sumU').*cos(n); y - B(row)*sign(sumU').*sin(n)]';

pcolor([x+B(row)/2*cos(n);x-B(row)/2*cos(n)],...
       [y+B(row)/2*sin(n);y-B(row)/2*sin(n)],...
       0*mean(unl.q./unl.h) + [-sumU, sumU]');
   
shading interp
colorbar
box on
hold on
axis equal
plot(x,y,'b')
drawnow
scatter(bank_points(:,1),bank_points(:,2),'b.')

fc = sample_fc_at_points(bank_points,pathrow(row));
fc=fc(:);

tab = table(fc, sumU, bheight);
dn = .2835.*fc.^-1.0767.*exp(0.14095.*bheight).*sumU;
dn=dn';

pcolor([x+B(row)/2*cos(n); x-B(row)/2*cos(n)],...
       [y+B(row)/2*sin(n); y-B(row)/2*sin(n)],...
       [-dn; dn]);
   
shading interp
colorbar
box on
hold on
axis equal
plot(x,y,'b')
drawnow

%% plot
%scatter(bank_points(:,1),bank_points(:,2),20,abs(dn),'filled')
% l1=color_line(x+B(row)/2*cos(n)-B(row)/2*(dn).*cos(n),y+B(row)/2.*sin(n)-B(row)/2*(dn).*sin(n),-dn);
% l1.LineWidth = 2;
% l2=color_line(x-B(row)/2*cos(n)-B(row)/2*(dn).*cos(n),y-B(row)/2.*sin(n)-B(row)/2*(dn).*sin(n),dn);
% l2.LineWidth=2;

% l1=color_line(x+B(row)/2*cos(n)-B(row)/2.*cos(n),y+B(row)/2.*sin(n)-B(row)/2.*sin(n),-dn);
% l1.LineWidth = 2;
% l2=color_line(x-B(row)/2*cos(n)-B(row)/2.*cos(n),y-B(row)/2.*sin(n)-B(row)/2.*sin(n),dn);
% l2.LineWidth=2;

colormap hot
ax=gca;
ax.CLim=[0 max(dn)];


end