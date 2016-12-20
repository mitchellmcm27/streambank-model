function []=make_graphics(data,i,plots)


if any(ismember([1 2 3 4 5],plots))
nhd = shaperead('K:\GIS\v\sites\lines_Project_Join_new.shp');
A = data.drainage_A_sqkm(i);
S = data.channel_slope(i);

[W, H] = model_geometry(A);
A_bkf = W.*H;
[Q, U, Cf] = model_average_velocity(W, H, S);
[x_site, y_site,~] = deg2utm(data.lat(i),data.long(i));

if isempty(plots)
    plots = 1:100;
end

for j=1:numel(nhd)
    if data.COMID(i) == nhd(j).COMID
        nhd = nhd(j);
        break
    end
end

dpx = nhd.X;
dpy = nhd.Y;
bbox = [min(dpx) max(dpx);
    min(dpy) max(dpy)];
%bbox = [526187 526817;
%        3412317 3413347];

quad = data.quad(i);


lidar = char(strcat('K:\GIS\r\elevation\lidar01m',quad,'.tif'));
ned = char(strcat('K:\GIS\r\elevation\ned03m',quad,'.tif'));
try
    DEM = GRIDobj(lidar);
    DEM = crop(DEM,bbox(1,:),bbox(2,:));
catch
    DEM = GRIDobj(ned);
    DEM = crop(DEM,bbox(1,:),bbox(2,:));
end

[xi,yi,k,n,s,bank_points] = model_curvature(dpx, dpy, W, x_site, y_site);

index = dsearchn([xi; yi]',[x_site y_site]); % closest index to our site

[Ub_nl, H_nb, tau_nb] = model_velocity_nonlinear(xi,yi,k,Q,W,H,U,Cf,S,1,0);
end

%% digitized points
if ismember(1, plots)
figure
imageschs(DEM,[],'colormap','parula');
hold on
scatter(x_site,y_site,'r+');
plot(dpx,dpy,'w.');
xlim(bbox(1,:));
ylim(bbox(2,:));
colorbar off
end

%% spline
if ismember(2,plots)
figure
imageschs(DEM,[],'colormap','parula');
hold on
scatter(x_site,y_site,'r+');
plot(xi,yi,'-w');
xlim(bbox(1,:));
ylim(bbox(2,:));
colorbar off
end

%% bank points
if ismember(3,plots)
figure
imageschs(DEM,[],'colormap','parula');
hold on
scatter(x_site,y_site,'r+');
plot(xi,yi,'-w');
hold on
b=plot(bank_points(:,1),bank_points(:,2),'.');
b.Color = [0.5 0.5 0.5 0.5];
xlim(bbox(1,:));
ylim(bbox(2,:));
colorbar off
end
%% velocity
if ismember(4,plots)
figure
%
% m = imageschs(DEM,[],'colormap','parula');
% hold on
% scatter(x_site,y_site,'r+');

pcolor([xi+W/2*cos(n); xi; xi-W/2*cos(n)],...
    [yi+W/2*sin(n); yi; yi-W/2*sin(n)],...
    U + [-Ub_nl, zeros(size(xi))', +Ub_nl]');
shading interp
cmap = rgbmap('blue', 'white', 'red', 'yellow',256);
colormap(cmap)
ax=gca;
ax.CLim = [0.1 0.8];
axis equal
box on;
title('velocity (m/s)')
xlabel('x-coordinate (m)')
ylabel('y-coordinate (m)')

deltan = (Ub_nl')*5*0;
hold on;
plot(xi+W/2.*cos(n)-W/2*(deltan).*cos(n), yi+W/2.*sin(n)-W/2*(deltan).*sin(n),'k-',...
    xi-W/2.*cos(n)-W/2*(deltan).*cos(n),yi-W/2.*sin(n)-W/2*(deltan).*sin(n),'k-');
xlim(bbox(1,:));
ylim(bbox(2,:));
end

%% flow depth
if ismember(5,plots)
figure
imageschs(DEM,[],'colormap','viridis');
hold on
scatter(x_site,y_site,'r+');
S = data.channel_slope(i);
[~,zi,~,~] = demprofile(DEM,numel(xi),xi,yi);
zi = double(zi(:));
pcolor([xi+W/2*cos(n); xi; xi-W/2*cos(n)],...
    [yi+W/2*sin(n); yi; yi-W/2*sin(n)],...
    H + [H_nb + zi, zeros(size(xi))' + zi, -H_nb + zi]');
shading interp
axis equal
end
%% precip and discharge
if ismember(6,plots)
figure
dates=data.dates{i};
Q_norm = data.Q_norm{i};
P_actual = data.P_monthly{i};
P_premonth = data.P_premonth(i);
P_norm = data.P_norm{i};


F = [P_premonth P_actual]./[P_norm(end) P_norm];
F_hi = F;
F_hi(F<1) = 1;
F = [F(2:end) .* 0.5.*(F_hi(1:end-1))];

Q_est = F .* Q_norm;

Q_est2 = data.AR_monthly{i}; % mm
Q_est2 = Q_est2./1000.*data.drainage_A_sqkm(i).*1000.*1000; % m3
Q_est2 = Q_est2./2592000; % m3/s

Q_est3 = data.nldas_runoff{i} + data.nldas_baseflow{i}; % kg/m2;
Q_est3 = Q_est3.*data.drainage_A_sqkm(i).*1000.*1000; % kg = L = dm3
Q_est3 = Q_est3./1000; % m3
Q_est3 = Q_est3./2592000; % m3/s

h1 = axes;
hold on
d=datenum(dates);
p=bar([d(1)-30; d],[P_norm(end) P_norm],'hist');
b=bar([d(1)-30; d],[P_premonth P_actual],0.3);

b.EdgeColor='none';
b.FaceColor=rgb('cloudy blue');
p.FaceColor = b.FaceColor./[2 2 1.5];
p.EdgeColor = 'none';
set(h1, 'Ydir', 'reverse')
set(h1, 'YAxisLocation', 'Right')
h1.YLim = [0 max(P_actual)*3];
h1.YTick = 0:100:(max(P_actual));
h1.XTick=[];
h1.YColor = p.FaceColor;
h1.Box = 'off';
hold on
ylabel('Precipitation (mm)','Rotation',-90,...
       'VerticalAlignment','baseline','HorizontalAlignment','center')
legend('Normal monthly P','Actual monthly P','Location','east')

h2=axes;
area(datenum(dates),Q_norm,'EdgeColor','none','FaceColor',[0.8 0.8 0.8])
hold on
plot(dates,Q_est,'.k-',dates,Q_est2,'.g-',dates,Q_est3,'.b-')

set(h2, 'Color', 'None')
h2.YLim=[0 max(Q_est)*1.7];
h2.YTick = 0:2:12;
h2.Box = 'off';
ylabel('Discharge (m^3/s)',...
    'VerticalAlignment','baseline','HorizontalAlignment','center');


try
    Q_actual = data.Q_monthly_USGS{i};
    plot(dates,Q_actual,'-r.');
end
legend('Normal monthly Q','Recursive Q','Budyko Q','NLDAS Q','Actual Q (USGS)',...
    'Location','southwest');
    h2.XLim = h1.XLim;
end
end