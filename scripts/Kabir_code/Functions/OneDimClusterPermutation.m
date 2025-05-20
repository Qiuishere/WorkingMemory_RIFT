function [sig,si,cs,Tmass_true,sig_mask] = OneDimClusterPermutation(MAT,min_cluster_size,permutations)
%1D Cluster permutation test
    %Parameters
    alpha = 0.05; disp("Here");
    
    %Calculate Surrogate Distribution
    
%     %%%Within subject ||  Randomly multiply by -1 or 1
%     MAT_perms = zeros(permutations,size(MAT,1),size(MAT,2));
%     for perms = 1:permutations
%         MAT_tmp = MAT;
%         for subj = 1:size(MAT,1)          
%             bollarr = ones(1,randi(size(MAT,2)));
%             bollarr= [bollarr zeros(1,size(MAT,2)-length(bollarr))]; %Create matrix of random indices
%             bollarr= bollarr(randperm(length(bollarr)));
%             MAT_tmp(subj,logical(bollarr)) = MAT_tmp(subj,logical(bollarr))*-1;
%         end
%         MAT_perms(perms,:,:) = MAT_tmp;
%     end
     
    %%%Within subject ||  Randomly multiply by -1 or 1
    MAT_perms = zeros(permutations,size(MAT,1),size(MAT,2));
    for perms = 1:permutations
        MAT_tmp = MAT;
        MAT_tmp = MAT_tmp.*((randi([0 1],1,size(MAT,1))*2)-1)';
        MAT_perms(perms,:,:) = MAT_tmp;
    end
    
    %Do t-tests and extract clusters
    Tmass_perm = zeros(1,permutations);
    for perms = 1:permutations
        MAT_PM = squeeze(MAT_perms(perms,:,:));
        [h,~,~,stats] = ttest(MAT_PM);
        f = find(diff([0,h,0]==1)); %Find consecutive significant clusters
        si = f(1:2:end-1);  % Find start indices of clusters
        cs = f(2:2:end)-si;  % Find size of clusters
        cs(cs<min_cluster_size) = [];
        si(cs<min_cluster_size) = [];
                
        %Find the fatest cluster with respect to T-Mass
        Tmass_temp = [];
        for clstr = 1:length(cs)
            Tmass_temp(clstr) = sum(stats.tstat(si(clstr):(si(clstr)+cs(clstr)-1)));
        end
        
        [~,idx] = max(Tmass_temp);
        Tmass_perm(perms) = sum(stats.tstat(si(idx):(si(idx)+cs(idx)-1)));%Calculate t-mass of fattest cluster
    end
    
    [h,~,~,stats] = ttest(MAT,0,'Alpha',0.05);
    f = find(diff([0,h,0]==1)); %Find consecutive significant clusters
    si = f(1:2:end-1);  % Find start indices of clusters
    cs = f(2:2:end)-si;  % Find size of clusters
    
    cs(cs<min_cluster_size) = [];
    si(cs<min_cluster_size) = [];
    
    Tmass_true = [];
    for clstr = 1:length(cs)
        Tmass_true(clstr) = sum(stats.tstat(si(clstr):(si(clstr)+cs(clstr)-1)));%Calculate t-mass of significant clusters
    end

    sig = [];
    for clstr = 1:length(Tmass_true)
        sig(clstr) = Tmass_true(clstr)>prctile(Tmass_perm,100-(alpha*100));%Compare true clusters to surrogate distribution
    end

    sig_mask = zeros(1,size(MAT,2));
    
    for clstr = 1:length(Tmass_true)
        if sig(clstr)
            sig_mask(si(clstr):(si(clstr)+cs(clstr)-1)) = 1;
        end
    end