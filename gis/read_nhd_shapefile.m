function [line] = read_nhd_shapefile(field,target_field,dir)

lines = shaperead(dir);

for j=1:numel(lines)
    if field == getfield(lines(j),target_field)
        line = lines(j);
        return
    end    
end

