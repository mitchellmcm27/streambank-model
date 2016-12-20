function [pathrow] = latlon2pathrow(lat, lon, dir)
old_dir = cd;
cd(dir)

lat = lat(:);
lon = lon(:);

shp = shaperead('wrs2_descending.shp');

c_X = cellfun(@nanmedian,{shp.X});
c_Y = cellfun(@nanmedian,{shp.Y});

nearest = dsearchn([c_X; c_Y]', [lon lat]);

pathrow = cell2mat({shp(nearest).PR})';

cd(old_dir)
end

