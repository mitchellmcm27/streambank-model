function [Q_norm_monthly,...
    U_norm_monthly,...
    P_norm_monthly,...
    PET_norm_monthly] = get_EROM(COMIDs)

dbf_dir = 'K:\GIS\v\NHDPlusV2Data\NHDPlusSA\NHDPlus03W\EROMExtension\';

Q_norm_monthly = nan(numel(COMIDs),12);
U_norm_monthly = nan(size(Q_norm_monthly));
P_norm_monthly = nan(size(Q_norm_monthly));
PET_norm_monthly = nan(size(Q_norm_monthly));

for m=1:12
    file=strcat(dbf_dir,'EROM_',num2str(m,'%02d'),'0001.dbf');
    if m==1
        [all_data,names]=dbfread(file);
        all_comids = all_data(:,find(strcmpi('COMID',names)));
        [matching] = find(ismember(cell2mat(all_comids), COMIDs));
        matching=matching(matching>0);
    end
    matching=matching(matching>0);
    [attributes] = dbfread(file,matching,{'Comid' 'Q0001E' 'V0001E' 'PPT0001' 'PET0001'});
    
    attributes = cell2mat(attributes);
    returned_ids = attributes(:,1);

    for i=1:numel(returned_ids)
        id = returned_ids(i);
        Q_norm_monthly(COMIDs==id,m) = deal(attributes(returned_ids==id,2)); % cfps
        U_norm_monthly(COMIDs==id,m) = deal(attributes(returned_ids==id,3)); % fps
        P_norm_monthly(COMIDs==id,m) = deal(attributes(returned_ids==id,4)); % mm
        PET_norm_monthly(COMIDs==id,m) = deal(attributes(returned_ids==id,5)); % mm
    end
    
end

Q_norm_monthly = Q_norm_monthly .* 0.0283168466; % ft3 per sec -> m3/s
U_norm_monthly = U_norm_monthly .* 0.3048; % ft per sec -> m/s

Q_norm_monthly = num2cell(Q_norm_monthly,2);
U_norm_monthly = num2cell(U_norm_monthly,2);
P_norm_monthly = num2cell(P_norm_monthly,2);
PET_norm_monthly = num2cell(PET_norm_monthly,2);
end

