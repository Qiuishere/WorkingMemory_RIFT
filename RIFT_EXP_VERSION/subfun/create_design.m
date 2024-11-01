function  design  = create_design(conditions,datap,nblocks,randomize) 
% create experimental designs with balanced conditions per trials/block
default('nblocks',1)
default('randomize',1)
n          = numel(conditions);
if mod(datap/nblocks,1)>0
    error('unequal number of data per block...');
end
conditions{n+1} = 1:(datap/nblocks);
conditions{n+2} = 1:nblocks;
ni         = numel(conditions);
order      = ni:-1:1 ;
[A{order}] = ndgrid(conditions{order}) ;
% concatenate
design     = reshape(cat(ni+1,A{:}),[],ni);
% block ordering
[~,idx]    = sort(design(:,end));
design     = design(idx,:);
if randomize
% shuffle within block
tmp        = NaN(size(design));
for k = 1:nblocks
    block_design            = design(design(:,end)==k,:);
    trials_in_block         = size(block_design,1);
    tmp(design(:,end)==k,:) = block_design(randperm(trials_in_block),:);
end
design     = tmp;
end
% remove the trial counter (datap)
design(:,end-1) = [];
% add trial list
design     = horzcat(design,(1:size(design,1))');