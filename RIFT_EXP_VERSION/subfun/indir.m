function d = indir(path)
%% remove current folder and parent folder entries ( '.', '..') from dir
if nargin<1
d           = dir;
else
d           = dir(path);
end
dcell       = struct2cell(d);
exclude     = strcmp(dcell(1,:),'.') | strcmp(dcell(1,:),'..') | contains(dcell(1,:),'.ini') | contains(dcell(1,:),'.DS_Store');
% contains(dcell(1,:),'.');
d(exclude)  = [];
end
    
