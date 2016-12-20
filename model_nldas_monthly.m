function [soil_moisture, runoff, baseflow] = model_nldas_monthly(lat,lon,date,basinid)
old_dir = cd

nldas_dir = 'K:/GIS/MODEL/input/NLDAS/NLDAS_NOAH0125_M.002/'
basin_dir = 'K:/GIS/MODEL/input/basins/'

dates = date + calmonths(0:11)';
[yrs, mos, ~] = datevec(dates);

date_strings = cellstr(num2str(yrs*100 + mos));
soil_moisture = NaN(size(date_strings));
runoff = NaN(size(date_strings));
baseflow = NaN(size(date_strings));
lai = NaN(size(date_strings));

for i = 1:numel(date_strings)
    cd(nldas_dir)
    name = strcat(nldas_dir,'NLDAS_NOAH0125_M.A',date_strings{i},'.002.grb');
  
    grbfile=ncgeodataset(name);
    lons = double(grbfile{'lon'}(:));
    lats = double(grbfile{'lat'}(:));
    R = makerefmat(lons(1),lats(1),mode(diff(lons)),mode(diff(lats)));
    [row, col] = latlon2pix(R,lat,lon);
    
    param = 'Soil_moisture_content_layer_between_two_depths_below_surface_layer_1_Month_Average';
    soil_moisture_grid = grbfile.data(param);
    % The second index is the layer, up to 4: 
    % 0-10cm, 10-40cm, 40-100cm, 100-200cm
    soil_moisture_grid = double(squeeze(soil_moisture_grid(1,1,:,:)));
    soil_moisture(i) = soil_moisture_grid(round(row),round(col)); % kg/m2 average
    
    param = 'Surface_runoff_non-infiltrating_surface_1_Month_Accumulation';
    runoff_grid = grbfile.data(param);
    runoff_grid = double(squeeze(runoff_grid(1,:,:)));
    runoff(i) = runoff_grid(round(row),round(col)); % kg/m2 accumulated 
    
    param = 'Subsurface_runoff_baseflow_surface_1_Month_Accumulation';
    baseflow_grid = grbfile.data(param);
    baseflow_grid = double(squeeze(baseflow_grid(1,:,:)));
    baseflow(i) = baseflow_grid(round(row),round(col)); % kg/m2 accumulated
    
    param = 'Leaf_area_index_0-9_surface_1_Month_Average';
    lai_grid = grbfile.data(param);
    lai_grid = double(squeeze(lai_grid(1,:,:)));
    lai(i) = lai_grid(round(row),round(col)); % 0-9 average
    
    % For runoff and baseflow, multiply by 1000 for ml
    
    h = imagesc(lai_grid);
    ax=gca;
    ax.YDir = 'normal';
    hold on
    axis equal
    scatter(col,row,'r.')
    
    figure
    h = imagesc(lons,lats,lai_grid);
    ax=gca;
    ax.YDir = 'normal';
    hold on
    axis equal
    scatter(lon,lat,'r.')
    break

%     cd(basin_dir)
%     basins = shaperead('Aug6')
%     
%     idx = find(strcmpi(num2str(basinid),{basins.BasinID}));
%     
%     basin = basins(idx);
%     bbox = basin.BoundingBox;
%     if sum(abs(bbox(:)) < 180) == 4
%         x = lon;
%         y = lat;
%     else
%         [x, y] = deg2utm(lat, lon);
%     end
%     bx = basin.X;
%     by = basin.Y;
%     
%     clipL = min(bbox(:,1)) - 0.1;
%     clipR = max(bbox(:,1)) + 0.1;
%     clipU = min(bbox(:,2)) - 0.1;
%     clipD = max(bbox(:,2)) + 0.1;
%     
%     [clipU, clipL] = map2pix(r,clipL,clipU)
%     [clipD, clipR] = map2pix(r,clipR,clipD)
%     p(round(clipU):round(clipD), round(clipL):round(clipR))
%     imagesc(p(clipU:clipD, clipL:clipR))
%     break
end

cd(old_dir)
end

