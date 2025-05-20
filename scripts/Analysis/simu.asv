flat = ones(200,1);
fs = 100;  % Sampling frequency (assumed to be 100 Hz)
T = 1/fs;  % Sampling period
N = 100;   % Number of data points
t = (0:N-1) * T; % Time vector
f = 10;
amplitude = 10;


for amplitude = 1:20
sinwave = amplitude*sin(2*pi*f*t)';
coswave = amplitude*cos(2*pi*f*t)';
simu_1 = [flat; sinwave; flat];
simu_2 = [flat; coswave; flat];

% simu_1 = [ sinwave; ];
% simu_2 = [coswave; ];

% figure
% plot(1:length(simu_1), simu_1); hold on
% plot(1:length(simu_2), simu_2);

fourier = fft(simu_1);
phase = angle(fourier);
amp = abs(fourier);

% figure
% plot(1:length(simu_1), phase); hold on
% plot(1:length(simu_1), amp);


fourier2 = fft(simu_2);
phase2 = angle(fourier2);
amp2 = abs(fourier2);


% 10 hz is at the 31 location is the signal is 3s, 51 if the signal is 5s
phases(amplitude, 1) = phase(51);
phases(amplitude, 2)  = phase2(51);

end



plot(1:size(phases,1), phases)


figure
plot(1:length(simu_1), phase); hold on
plot(1:length(simu_1), phase2);
