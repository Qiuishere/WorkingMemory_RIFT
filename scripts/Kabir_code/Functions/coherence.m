function [coh] = coherence(signalX,signalY)

    n = size(signalX,2); %number of trials
    %signalY = repmat(signalY,n,1);

    %Hilbert
    signalXh = hilbert(signalX);
    signalYh = hilbert(signalY');

    %Magnitude
    mX = (abs(signalXh))';
    mY = abs(signalYh);

%     %Phase
%     pX = angle(signalXh);
%     pY = angle(signalYh);

    %phase_diff = pX-pY;
    %phase_diff = angle(signalXh'-signalYh);
    phase_diff = angle(signalXh')-angle(signalYh);

    coh = [];
    for t = 1:length(signalY)
       %coh(t) = (abs(mX(t)*mY(t)*exp(1)^(1i*phase_diff(t)))^2)/((mX(t)^2)*(mY(t)^2)); 
       coh(t) = (  abs(sum(mX(:,t).*mY(:,t).*exp(1).^(1i*phase_diff(:,t)))/n)  .^2)  /  (sum(  (mX(:,t).^2).*(mY(:,t).^2)  )/n);
       %coh(t) = abs(sum(signalXh(:,t).*signalYh(:,t).*exp(1).^(1i*phase_diff(:,t)))/n) / ((sum((signalXh(:,t).^2))/n).*(sum((signalYh(:,t).^2))/n));
    end



