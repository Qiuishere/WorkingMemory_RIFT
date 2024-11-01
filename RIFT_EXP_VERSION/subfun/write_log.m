function write_log(file_name,dset)
% function to write row-wise matlab table/dataset to .txt
% remove non-numeric variables, added as additional columns at the end
if verLessThan('matlab', '8.3')
    VNAMES        = dset.Properties.VarNames;
else
    VNAMES        = dset.Properties.VariableNames;
end
checknum          = NaN(numel(VNAMES),1);
for k = 1:size(dset,2)
    checknum(k,:) = ~isnumeric(dset.(VNAMES{k}));
end
listout           = find(checknum);
dsetout           = dset(:,listout);
dsetin            = dset;
dsetin(:,listout) = [];
VNAMES_new        = VNAMES([find(checknum==0)' listout']);
% initialize the file adding the header if no such file exists
% if isempty(ls(file_name))
if ~exist(file_name,'file')
    cell2csv(file_name,VNAMES_new,' ');
    fid               = fopen(file_name,'a');
    fprintf(fid,'\n');fclose(fid);
end    
if verLessThan('matlab', '8.3')
dsetin            = double(dsetin);
dsetout           = dataset2cell(dsetout)';
if ~isempty(dsetout)
dsetout           = {dsetout{2}};
end
else
dsetin            = table2array(dsetin);
dsetout           = table2cell(dsetout);
end

fid               = fopen(file_name,'a');
fprintf(fid, [repmat('%d\t',1,size(dsetin,2)) ...
              repmat('%s\t',1,size(dsetout,2)) '\n'], dsetin, dsetout{:});
fclose(fid);
