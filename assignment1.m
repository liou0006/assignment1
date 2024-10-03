clear

% Parameters
fs = 44100; % Sampling frequency in Hz
f0 = 1000; % Frequency of the main peak in Hz
T0 = 0.5; % Duration in seconds
t = 0:1/fs:T0-1/fs; % Time vector for T0 seconds
A1 = 1; % Amplitude for synthetic signal

% Synthetic signal
signal = A1*sin(2*pi*f0*t);

% Paramters for the recorded voice
Fs = fs; %% sample freq
nBits = 16; 
nChannels = 2; 
ID = -1;       % default audio input device 
recObj = audiorecorder(Fs,nBits,nChannels,ID);
disp("Begin speaking.")
recDuration = 0.5; 
recordblocking(recObj,recDuration);
disp("End of recording.")
play(recObj);
audioData = getaudiodata(recObj);

% Extract the first channel of the recorded audio (getaudiodata gives you
% a 2D array. We only want one
audioData_y = audioData(:, 1);

% Ensure both signals have the same length
min_length = min(length(t), length(audioData_y));
signal = signal(1:min_length);
audioData_y = audioData_y(1:min_length);

% Plot the time-domain signals
figure(1);
clf
h1 = plot(t(1:min_length), signal);
hold on
h2 = plot(t(1:min_length), audioData_y);
xlim([0.25 0.27])
xlabel('Time (s)')
ylabel('Magnitude (V)')
title('Magnitude Spectrum in Time Domain')
legend([h1, h2], {'1000 Hz', 'Recorded audio'}) % Use plot handles for legend
grid on
saveas(gcf, 'Time Domain.png');

% Compute the magnitude spectrum
n = length(signal);
Y2 = fft(signal);
Y = Y2(1:n/2+1); % Take only the positive frequencies
Y = Y / n; % Normalize the FFT output

n1 = length(audioData_y);
Y3 = fft(audioData_y);
Y1 = Y3(1:n1/2+1); % Take only the positive frequencies
Y1 = Y1 / n1; % Normalize the FFT output

% Frequency vector
f = (0:n/2)*(fs/n);

abs1 = 20*log10(2*abs(Y1));

% Plot the magnitude spectrum
figure(2);
clf
h1=plot(f, 20*log10(2*abs(Y)));
hold on
h2 = plot(f, abs1);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude Spectrum in Frequency Domain')
xlim([0 3000]) % Limit the x-axis to 3000 Hz
legend([h1, h2], {'1000 Hz', 'Recorded audio'}) % Use plot handles for legend
grid on
saveas(gcf, 'Frequency Domain.png');

