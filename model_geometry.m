function [W,D,Q,U] = model_geometry(A,region)

A_mi = A.*0.386102;


region=1; 
if region==1
    % Northwest Florida Coastal Plain
    Wbkf_ft = 10.4*A_mi.^0.39; % ft Metcalf 2009
    Dbkf_ft = 1.64*A_mi.^0.25; % ft Metcalf 2009
    Qbkf_ft = 27.7.*A_mi.^0.71; % cfs
else
    % North Florida Coastal Plain
    Wbkf_ft = 9.2*A_mi.^0.28; % ft Metcalf 2009
    Dbkf_ft = 0.67*A_mi.^0.43; % ft Metcalf 2009
    Qbkf_ft = 7.54.*A_mi.^0.77; % cfs
end

W = Wbkf_ft .* 0.3048; % bankfull width meters
D = Dbkf_ft .* 0.3048; % bankfull mean depth meters
Q = Qbkf_ft .* 0.028316847; % m3/s
U = Q./W./D;

% Bieger et al JAWRA 2015
% SE US sand bed streams
% W = 2.22 * A.^0.363; % R2 = 0.84
% D = 0.24 * A.^0.323;

% % Faustini et al. Geomorphology 2009
% Ecoregion 03 - South Atlantic-Gulf
W = 2.4.*A.^0.36; %<-- very similar to Bieger et al

W=W(:);
D=D(:);
Q=Q(:);
end

