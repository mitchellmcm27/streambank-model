
function [output] = train_model_monthly(input_table,varargin)

addpath('utils')
addpath('bdv model')
addpath('gis')

%% Process args
 plt=0;
 animate=0;
if ~isempty(varargin)
    for arg=1:length(varargin)
        switch(varargin{arg})
            case 'plot'
                plt = 1;
            case 'animate'
                animate = 1;
                plt=1;
        end
    end
end


%% Basic calculations all sites
input_lat = input_table.lat;
input_lon = input_table.long;
[input_x, input_y] = deg2utm(input_lat,input_lon); % returns column vectors

start_dates = input_table.date_y1; % initial measurement date = model simulation start date
dates_y2 = input_table.date_y2; % 1 year later
dates_y3 = input_table.date_y3; % final measurement date

end_dates = dates_y3; % assume model simulation end date is the final measurement date
end_dates(isnat(end_dates)) = dates_y2(isnat(end_dates)); % for sites that only have 1 year of data, take the 2nd year as end date (not applicable anymore)

% Round dates to the nearest first day of the month
% e.g. January 3rd becomes January 1st
% but January 28th becomes February 1st.
% The cutoff is the 15th (16th goes to first of next month).
start_dates = round_dates(start_dates); % inclusive
end_dates = round_dates(end_dates); % exclusive
input_dates = [start_dates end_dates]; % concatenate into a matrix with 2 cols

A = input_table.A; % drainage area sqkm
S = input_table.S; % downstream slope km/km
region = input_table.region; % code that determines which regional curve to use

[B, H, Q, U] = model_geometry(A,region); % regional curves

A_bkf = B.*H; % Bankfull cross-section area
Rh = A_bkf./(H + B + H); % Bankfull hydrualic radius
manning_n = 0.087; % assumed constant over study area
Cf0 = 9.81*manning_n.^2./(Rh.^(1/6)); % Friction factor (dimensionless)

pathrow=input_table.pathrow; % Landsat (WRS-2) path/row

bd = input_table.bd; % SSURGO bulk density
psand = input_table.psand; % SSURGO percent sand
K = input_table.K; % SSURGO erodibility factor
f=input_table.f_monthly_2; % Monthly storm frequency over 2 years
fc = input_table.fc; % Tree cover 

%% Allocate empty vectors for intermediate calculations

% Cell-vector for each site, will hold no. of days in each month fo
delta_t = cell(size(input_table,1),1);

% Cell-vector to hold near-bank velocity excess for each site
U_bank_monthly = cell(size(input_table,1),1);

% Near-bank depth excess
H_bank_monthly = cell(size(input_table,1),1);


%% Loop over each site (rows in input table)
for row = 1:size(input_table,1);
    try
        % only process rows that have gis == true, which means 
        % they are suitable for the model
        assert(input_table.gis(row))
    catch
        continue
    end
    
    disp('------------------------------------------')
    disp(strcat('Site: ',input_table.Properties.RowNames(row)))
    disp('------------------------------------------')
    disp('')
    disp('Cross-section Geometry:')
    disp(strcat('... Width:',32,num2str(B(row))))
    disp(strcat('... Depth:',32,num2str(H(row))))
    
    % Determine simulation dates based on length of observation (months)
    % For each row, this creates a vector of monthly dates between the 
    % start and end dates
    model_dates = input_dates(row,1):calmonths(1):(input_dates(row,2)-calmonths(1)) + caldays(14);
    
    %% Read shapefile and calculate curvature
    
    % Read shapefile into a matlab struct
    centerline_shp = read_nhd_shapefile(input_table.COMID(row),'COMID','data/NHD/lines_Project_Join_new.shp');
    
    % Get UTM X and Y coords from shapefile
    % DP = digitized points
    dpx = centerline_shp.X(1:end-1);
    dpy = centerline_shp.Y(1:end-1);
    
    % Curvature calculation
    % x,y: resampled channel centerline points, UTM
    % iR: curvature (1/R or inverse of R, which is radius of curvature)
    % n: unit normal vector at each (x,y) location
    % s: downstream coordinate, starts at 0
    % bank_points: points 1 channel-width outside of the centerline, used
    %   to sample raster data
    [x,y,iR,n,s,bank_points] = model_curvature(dpx, dpy, B(row), input_x(row), input_y(row));
    
    %% Select location for sampling raster data
    % Nearest centerline point to field GPS coordinate
    idx_cl = dsearchn([x; y]', [input_x(row) input_y(row)]);
    
    % Nearest bank_point to field GPS coordinate
    idx_bank = dsearchn(bank_points, [input_x(row) input_y(row)]);
    bank_point = bank_points(idx_bank,:);
    
    %% Raster data
    
    % Could have a function to sample SSURGO rasters...
    % bd(row) = model_soils(bank_point);
    
    % But I broght them into the input table manually
    bd(row)=input_table.bd(row);
    
    % Sample tree cover raster
    % requires mapping toolbox!
    fc(row) = model_landsat(bank_point,pathrow(row),input_dates(row,:),'FC');
    
    %fc(row) = input_table.fc(row); % If you don't have mapping toolbox!
    
    % At one time I also used EVI/NDVI from landsat, this is how
    % evi{row} = model_landsat(bank_point,pathrow(row),input_dates(row,:),'EVI');
        
    % SM01Y1 had inaccurate tree cover data, this assigns a value to it
    if strcmpi(input_table.Properties.RowNames(row),'SM01Y1')
        fc(row)=0.05;
    end
        
    % Land cover - no longer used
    % nlcd(row) = model_land_cover(bank_point(:,1), bank_point(:,2),'16 N');
    % mcd(row) = mcd_lookup(nlcd(row));
    
    %% 1-D Streamflow model
    % Model monthly discharge and flow depth
    
    % Monthly flow volume from Noah LSM
    %
    % runoff + baseflow, averaged over entire basin
    V_noah = input_table.noah_runoff_2{row} + input_table.noah_baseflow_2{row}; % kg/m2;
    V_noah(V_noah<0) = 0; % give negative months a value of 0
    
    % Multiply by drainage area and convert units
    V_noah = V_noah.*A(row).*1000.*1000; % kg (or L, or dm^3)
    V_noah = V_noah./1000; % m^3
    
    % Discharge calculation
    %
    Q_dates = input_table.dates_2{row}; % valid dates for observations at this site
    d=datenum(Q_dates); % used for plotting
    ndays = eomday(year(Q_dates),month(Q_dates)); % number of days in each month
    
    % number of wet days = storm frequency (PRISM) * no. days
    wet_days = f{row}.*ndays;

    % Assuming events last 1 day and each event is the same magnitude
    Q_est = V_noah./wet_days./24./60./60; % estimated event Q m3/s
    delta_t{row} = wet_days(1:numel(model_dates))*24*60*60; % assuming 1 event = 1 day, this is the total length of storm events in seconds
    
    % Monthly flow depth
    % assuming steady uniform flow
    H_est = model_normal_depth(Q_est,H(row),B(row),S(row),manning_n);
    Rh_est = H_est.*B(row)./(2.*H_est + B(row)); % monthly hydraulic radius for given monthly depths
    
    % Estimate channel roughness assuming constant Manning's n
    Cf_est = 9.81*manning_n.^2./(Rh_est.^(1/6)).^2; % Cf = gn^2/Rh^(1/36)
    
    %% Plot data up to now   
    if plt
        fig=figure(row);
        fig.Color = 'w';
        suptitle(input_table.Properties.RowNames{row})
        h1=subplot(2,2,1);
        hold on
        base_color = rgb('cloudy blue');
        dark_color = base_color./[2 2 1.5];
        
        d=datenum(Q_dates);
        
        dt =wet_days./7; % weeks
        b1=bar(d,Q_est);
        b1.EdgeColor='none';
        b1.FaceColor=base_color;
        plot(d,H_est,'k-');
        plot([min(Q_dates) max(Q_dates)],[H(row) H(row)],'k:')
        h2.XLim=h1.XLim;
        h2.Color='None';
        h2.Box = 'off';
        legend('Q','h','h_{bkf}','\Delta t','EVI','Location','north','Orientation','horizontal','FontName','Myriad Pro')
        p1 = scatter(NaN,NaN,'or');
        p2 = scatter(NaN,NaN,'.r');
    end
    
    %% Flow Modeling
    
    % clip the centerline to a reasonable length
    % expressed in node indices (=1/10 width)
    len = 160;
    low = max(1,idx_cl-len); % upstream limit: "len" indices
    hi = min(idx_cl+len/2,numel(x)); % downstream limit: "len/2" indices
    range = low:hi;
    
    % clip vectors to these limits
    x=x(range);
    y=y(range);
    s=s(range);
    iR=iR(range);
    
    if animate % create .gif filename for this site
        giffile = strcat('gifs/',input_table.Properties.RowNames{row},'.gif');
    end
    
    %% Loop over simulation months (still in loop over sites)
    
    for mnth = 1:numel(model_dates)      
        if mnth==1
            % First month of simulations...
            
            
            % simulate bankfull conditions to set bed topography
            AR0 = 0;
            [unl,~,~,~] = ...
            model_velocity_nonlinear_swe(x,y,iR,Q(row),B(row),H(row),s,Cf0(row),S(row),1,AR0,1);
            
            % get bed topography for subsequent runs
            AR0 = unl.AR;
        end
        
        % Allow flow Q and depth to vary each month, 
        % keep bed topography from the bankfull condition
        
        % Blanckaert, DeVriend, and Ottevanger's model
        %  Inputs
        %    x,y: centerline geom
        %    iR: curvature
        %    Q_est(mnth): this month's estimated Q
        %    B: bankfull width
        %    H_est(mnth): this month's estimated normal depth
        %    s: downstream coordinate
        %    Cf_est(mnth): this month's estimated friction factor
        %    S: downstream slope'
        %    1: nonlinear model (0=linear)
        %    AR0: bed topography (scour factor) modeled with mnth=1
        %    3: bed type, use AR from input (see meander_model_copy)
        %  Outputs
        %    unl: struct containing lots of data
        %    n: normal vector
        [unl,~,~,n] = ...
            model_velocity_nonlinear_swe(x,y,iR,Q_est(mnth),B(row),H_est(mnth),s,Cf_est(mnth),S(row),1,AR0,3);
        
        % Select nearest centerline point (required again 
        % because we clipped the centerline)
        idx_cl = dsearchn([x; y]', [input_x(row) input_y(row)]);
        
        % From BdV,
        % Unb = U * as/R * B/2
        % Hnb = H* A/R * B/2
        U_nb = unl.q./unl.h.*unl.asR.*B(row)/2;
        H_nb = H(row).*(unl.AR+unl.Fr2R).*B(row)/2;
        BR = B(row)*unl.iR(:); % B/R(s)
        
        % Make sure we are using the correct bank!!
        if input_table.bank(row)=='Left'
            U_bank_monthly{row}(mnth) = U_nb(idx_cl);
            H_bank_monthly{row}(mnth) = H_nb(idx_cl);
        else
            U_bank_monthly{row}(mnth) = -U_nb(idx_cl);
            H_bank_monthly{row}(mnth) = -H_nb(idx_cl);
        end
        
        %% Plot data 
        if plt
            figure(row)
            
            subplot(2,2,1)
            p1.XData = d(mnth);
            p1.YData = H_est(mnth);
            p2.XData = d(1:mnth-1);
            p2.YData = H_est(1:mnth-1);
            
            subplot(2,2,2)
            axis equal
            hold on
            pcolor([x+B(row)/2*cos(n);x-B(row)/2*cos(n)],...
                [y+B(row)/2*sin(n);y-B(row)/2*sin(n)],...
                mean(unl.h) + [-H_nb, H_nb]');
            shading interp
            colormap parula
            colorbar
            box on;
            title('Water depth (m)','FontWeight','normal','FontName','Myriad Pro')
            xlabel('UTM e')
            ylabel('UTM b')
             
            scatter(x(idx_cl),y(idx_cl),20,'bo')
            scatter(bank_point(:,1),bank_point(:,2),'ob')
            scatter(input_x(row),input_y(row),'b.');
            
            ax=gca;
            cmap1 = rgbmap('ice blue','cornflower blue', 'dark slate blue');
            cmap2 = rgbmap('army green','cement','off white');
            cmap = [cmap2; cmap1];
            colormap(cmap)
            ax.CLim = [-H(row)*3 H(row)*3];
            ax.XTick=[];
            ax.YTick=[];
            
            subplot(2,2,3)
            hold on
            axis equal
            pcolor([x+B(row)/2*cos(n);x-B(row)/2*cos(n)],...
                [y+B(row)/2*sin(n);y-B(row)/2*sin(n)],...
                mean(unl.q./unl.h)*0 + [-U_nb, U_nb]');
            shading interp
            colorbar
            axis equal
            box on
            scatter(x(idx_cl),y(idx_cl),20,'bo')
            scatter(bank_point(:,1),bank_point(:,2),'ob')
            scatter(input_x(row),input_y(row),'b.');
            
            title('\Delta U (m/s)','FontWeight','normal')
            xlabel('UTM e')
            ylabel('UTM n')
            ax=gca;
            ax.CLim = [-U(row)*1.5 U(row)*1.5];
            ax.XTick=[];
            ax.YTick=[];

            figure(row)
            subplot(2,2,4)
            p=plot(s,BR);
            p.Color=[0.6 0.6 0.6];
            hold on
            plot(s,U_nb./unl.q./unl.h,'r-',s,H_nb./unl.h,'b-')
            drawnow;
            legend('B/R','\Delta U/U','\Delta H/H','FontName','Myriad Pro')
            axis tight
            hold off
            
            if animate
                fig=figure(row);
                fig.Position=[200 150 1440 900];
                F = getframe(fig);
                im = frame2im(F);
                [aX,map] = rgb2ind(im,256);
                if mnth==1
                    imwrite(aX,map,giffile,'gif','Loopcount',inf,'DelayTime',.5)
                else
                    imwrite(aX,map,giffile,'gif','Writemode','append','DelayTime',.5)
                end
            end
        end
        drawnow;
    end
    if plt
        figure
        plot(model_dates,U_bank_monthly{row}*10,'r',model_dates,H_bank_monthly{row},'b')
        legend('10\Delta U','\Delta h','FontName','Myriad Pro')
        hold on
        plot([model_dates(1) model_dates(end)],[mean(U_bank_monthly{row}*10) mean(U_bank_monthly{row}*10)],'r:')
        plot([model_dates(1) model_dates(end)],[mean(H_bank_monthly{row}) mean(H_bank_monthly{row})],'b:')
        
    end
end
drawnow; % draw all figures after all sites/months have been calculated



%% Model parameter fitting

y = input_table.erosion_rate./100; % m/yr
y(y<0)=0; % negative rates (deposition) are interpreted as zero erosion

% Take the mean of monthly deltaU and deltaH
U_bar = nan(size(input_table,1),1);
H_bar = nan(size(U_bar));
for row=1:size(input_table,1)
      U_bar(row) = mean(U_bank_monthly{row});
      H_bar(row) = mean(H_bank_monthly{row});
end

% negative velocity excess implies deposition, which we are interpreting as
% zero erosion (we aren't attempting to model deposition at all)
U_bar(U_bar<0)=0;

% total bank height = mean depth + scour
H_bank = H_bar+H;

% Create an output table
% Some of these variables aren't used but may be interesting
output_table = table(fc, bd, A, psand, K, U_bar, H_bar, S, H_bank, y);
output_table.Properties.RowNames=input_table.Properties.RowNames;

% subsample the output table for only rows that were included in the model
output_table = output_table(input_table.gis,:);
disp(output_table);

% Setup statistical model options
opts = statset('nlinfit');
opts.MaxIter=2000;

% Erodibility k1
mdl1 = fitnlm(output_table, ...
     @(b,x) b(1) .* (exp(b(2).*log(x(:,1)) + b(3).*x(:,9) + b(4).*x(:,2).^2.5)) .* (x(:,6)), ...
     [0.1 -1 0.5 -.5],'Options',opts,'ErrorModel','proportional');

% Erodibility k2
mdl2 = fitnlm(output_table, ...
    @(b,x) b(1) .* x(:,1).^b(3).*exp(b(2).*x(:,9)).*x(:,6),...
    [1 -1 -0.5 ],'Options',opts,'ErrorModel','proportional');

% Final plots of results
figure
gscatter(predict(mdl1),mdl1.Variables.y,output_table.Properties.RowNames)
plot1;
ax=gca; ax.XScale='log';ax.YScale='log';
legend off
xlabel('Simulated')
ylabel('Observed')

figure
gscatter(predict(mdl2),mdl2.Variables.y,output_table.Properties.RowNames)
plot1;
ax=gca; ax.XScale='log';ax.YScale='log';
legend off
xlabel('Simulated')
ylabel('Observed')

figure
subplot(1,2,1)
scatter(predict(mdl1),mdl1.Variables.y,'ko','filled')
plot1;
ax=gca; ax.XScale='log';ax.YScale='log';
ax.XLim = [0.001 5];
ax.YLim = [0.001 5];
legend off
xlabel('Predicted erosion rate (K=K_1)')
ylabel('Observed erosion rate (m/yr)')
plot1

subplot(1,2,2)
scatter(predict(mdl2),mdl2.Variables.y,'ko','filled')
plot1;
ax=gca; ax.XScale='log';ax.YScale='log';
ax.XLim = [0.001 5];
ax.YLim = [0.001 5];
legend off
xlabel('Predicted erosion rate (K=K_2)')
ylabel('Observed erosion rate (m/yr)')
plot1

% Output data
output.tbl = output_table;
output.mdl1 = mdl1;
output.mdl2 = mdl2;

end