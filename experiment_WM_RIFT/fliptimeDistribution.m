dur = diff(fliptimes(1,:));

%scatter(ones(length(dur),1), dur)
figure;

%subplot(1,3,1)


scatter(1*ones(121,1), dur(144:264)); hold on % target
scatter(2*ones(101,1), dur(264:364))          % mask
scatter(3*ones(301,1), dur(400:700))          % delay
