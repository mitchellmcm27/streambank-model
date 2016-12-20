function [ basin ] = model_prism_monthly(lat,lon,date,basinid)
old_dir = cd

prism_dir = 'K:/GIS/MODEL/input/PRISM/precip/monthly/'
basin_dir = 'K:/GIS/MODEL/input/basins/'

dates = date + calmonths(0:11)';
[yrs, mos, ~] = datevec(dates);

date_strings = cellstr(num2str(yrs*100 + mos));

for i = 1:numel(date_strings)
    cd(prism_dir)
    wc = strcat('*',date_strings{i},'_asc.asc');
    files = dir(wc);
    name = extractfield(files, 'name')';
    
    assert(numel(name)==1)
    
    [p,r]=arcgridread(char(name));
    
    cd(basin_dir)
    basins = shaperead('Aug6')
    
    idx = find(strcmpi(num2str(basinid),{basins.BasinID}));
    
    basin = basins(idx);
    bbox = basin.BoundingBox;
    if sum(abs(bbox(:)) < 180) == 4
        x = lon;
        y = lat;
    else
        [x, y] = deg2utm(lat, lon);
    end
    bx = basin.X;
    by = basin.Y;
    
    clipL = min(bbox(:,1)) - 0.1;
    clipR = max(bbox(:,1)) + 0.1;
    clipU = min(bbox(:,2)) - 0.1;
    clipD = max(bbox(:,2)) + 0.1;
    
    [clipU, clipL] = map2pix(r,clipL,clipU)
    [clipD, clipR] = map2pix(r,clipR,clipD)
    p(round(clipU):round(clipD), round(clipL):round(clipR))
    imagesc(p(clipU:clipD, clipL:clipR))
    break
end

cd(old_dir)
end

