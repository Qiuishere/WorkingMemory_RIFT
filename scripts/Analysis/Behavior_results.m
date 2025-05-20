addpath('/project/3018085.01/scripts/subfun')
Start_up
subjects = datainfo();

stim = prm.fac.stimuli;
% read in behaviour data
for thesub = [1,2,3,5,6]

    load(subjects(thesub).behav);
    Behav = T;
    for thestim = 1:2
        trialid = strcmp(T.stimulus, stim{thestim}) ;
        meanError(thesub,thestim) = mean(T.Error(trialid)); % only calculate acc for RS trials
        STD(thesub,thestim) = std(T.Error(trialid));
    end

end


bar([1,3], meanError); hold on
errorbar([1,3], meanError, STD)

    
%% ===========plot===============
MarkerSize = 5;
Col.pink = [255 153 153] /255;
Col.green = [91 153 153] /255;
Col.black = [0 0 0];
Col.grey = [.4 .4 .4];
Col.white = [1 1 1];


todraw = meanError;
nbars = size(todraw,2);
groupwidth = min(0.8, nbars/(nbars+1.5));

x = [1-0.6 1+0.6];
figure('OuterPosition',[0 0 350 600]);
hold on;
%x = [0.5 1 ];x = x(1:nbars);
for ii = 1:nbars
    [Group_mu(ii),Group_sigma(ii), Group_MUCI(ii,:), Group_SigmaCI(ii,:)] = normfit(todraw(:,ii),0.05);
    cilength(ii) = Group_MUCI(ii,2) - Group_ave(ii);

    barhd(ii) = bar(x(ii),Group_mu(ii) );
    scathd = plot(x(ii),todraw(:,ii),'o','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.black,'MarkerSize',MarkerSize);
    errorbar(x(ii), Group_mu(ii), cilength(ii),  'k','linestyle','none','linewidth',1.2,'alignvertexcenters','on');

end
%plot(ones(size(todraw,2))*x,[alPSE{todraw(1)} alPSE{todraw(2)}],'o','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.black,'MarkerSize',MarkerSize);
%plot(ones(size(alPSE{todraw(1)}))*x,[alPSE{todraw(1)} ],'o','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.black,'MarkerSize',MarkerSize);
plot(x,todraw,':','LineWidth',0.5 ,'color',[.7 .7 .7]); %k: black;  -: solid line ; d:diamond
set(gca, 'YDir', 'reverse')

set(gca,'XTick',x, 'XTickLabel', stim,'fontsize',12);
barhd(1).FaceColor= Col.grey;
barhd(2).FaceColor= Col.white;
barhd(1).EdgeColor= Col.grey;
barhd(2).LineWidth = 1;

%title(ExpName(expid) ,'fontsize',16); %% specify the figure name

ylabel('Error (degree)','fontsize',12);
%set(gca,'ylim',[-0.12 0.09],'yTick',-0.12:0.03:0.09,'fontsize',12);
box off
axis ij