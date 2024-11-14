dt = 1 / 1440; % in second
tag_freq = 71;%[2 4 8 1];
timax = dt:dt:2;
tag_sigs = cos(2*pi*tag_freq*timax) / 2 + 0.5;

figure
subplot(211)
plot(tag_sigs)
xlim([0 200])

subplot(212)
ft = abs(fft(tag_sigs));
hz = linspace(0,1440/2,floor(length(timax)/2)+1);
plot(hz,ft(1:length(hz)))