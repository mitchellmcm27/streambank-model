function [line] = read_nhd_shapefile(value,target_field,dir)

lines = m_shaperead(dir);
string(lines.(target_field));
idx = find(strcmp(string(value), string(lines.(target_field))));
line.X = lines.ncst{idx}(:,1);
line.Y = lines.ncst{idx}(:,2);

