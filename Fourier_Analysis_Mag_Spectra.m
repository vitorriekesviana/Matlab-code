clc;
clear;
% Load file
file_name=input('Enter path location:')
data = McsHDF5.McsData(file_name);
% Split the file_name by the file separator
parts = strsplit(file_name, '\');
% Extract the relevant parts (5th and 7th parts)
part1 = parts{6};
part2 = parts{8};
% Concatenate the parts with a '/' separator
desired_path = [part1, '/', part2];
% Display the result
disp(desired_path);
% Initialize full matrix
Matrix_full = data.Recording{1, 1}.AnalogStream{1, 1}.ChannelData;
time_micro_s = double(data.Recording{1, 1}.Duration);
time_s = double(time_micro_s / 10^6);
num_samples_full = size(Matrix_full, 2);
num_channels_full = size(Matrix_full, 1);
sampling_frequency = num_samples_full / round(time_s);
disp('The sampling frequency is')
disp(round(sampling_frequency))
% Break Full matrix in a wanted time interval
tempo_initial = 120;
tempo_final = 130;
x = round(double(tempo_initial * sampling_frequency));
y = round(double(tempo_final * sampling_frequency));
Data_analysis = Matrix_full(:, x:y);
% Parameters for new/reduced matrix
num_samples_reduced = size(Data_analysis, 2);
num_channels_reduced = size(Data_analysis, 1);
Fs = double(sampling_frequency);
T = 1 / Fs;                   % Sampling period (s)
L = num_samples_reduced;     % Length of signal
t = (0:L-1) * T;              % Time vector
% Check if L is even
if mod(L, 2) == 1
   % If odd, change to the nearest even number
   L = L - 1;
end
% Filter Signal
f_low = 1; % Lower cutoff frequency in Hz
f_high = 4; % Upper cutoff frequency in Hz
order = 2; % Filter order
Data_analysis1 = zeros(size(Data_analysis)); % Initialize Data_analysis1
for i = 1:size(Data_analysis, 1)
   Wn = [f_low, f_high] / (Fs / 2);
   % Create Butterworth band-pass filter
   [b, a] = butter(order, Wn, 'bandpass');
   % Apply the filter
   Data_analysis1(i, :) = filtfilt(b, a, Data_analysis(i, :));
end
% Plot the magnitude spectrum of all channels separately
figure;
hold on;
n_colors = num_channels_reduced; % Number of colors needed
colors = hsv(n_colors);
for i = 1:num_channels_reduced
   % Perform Fourier Transform
   X = fft(Data_analysis1(i, :));
   P2 = abs(X / L);            % Two-sided spectrum
   P1 = P2(1:L/2 + 1);         % Single-sided spectrum
   P1(2:end-1) = 2 * P1(2:end-1);
   f = Fs * (0:(L/2)) / L;     % Frequency domain
   % Plot the magnitude spectrum
   plot(f, P1, 'Color', colors(i, :));
end
% Formatting plot
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of All Channels');
xlim([0, 10]); % Adjust x-axis for better visibility
hold off;
% Compute the magnitude spectrum for each channel and average them
magnitude_spectra = zeros(num_channels_reduced, L/2 + 1);
for i = 1:num_channels_reduced
   % Perform Fourier Transform"G:\Marilyn_MEA\MEA Recordings\Data Analysis\h5 data raw (Matlab)\N102312\39070\39070_E1_div7\2023-12-01T17-23-52McsRecording.h5"
   X = fft(Data_analysis1(i, :));
   P2 = abs(X / L);            % Two-sided spectrum
   P1 = P2(1:L/2 + 1);         % Single-sided spectrum
   P1(2:end-1) = 2 * P1(2:end-1);
   magnitude_spectra(i, :) = P1;
end
% parentpectrum
avg_magnitude_spectrum = mean(magnitude_spectra, 1);
% Smooth the average magnitude spectrum using a moving average filter
windowSize = 10; % Adjust this value for more or less smoothing
avg_magnitude_spectrum_smoothed = movmean(avg_magnitude_spectrum, windowSize);
% Convert the smoothed spectrum to dB
avg_magnitude_spectrum_smoothed_dB = 20 * log10(avg_magnitude_spectrum_smoothed);
% Frequency domain
f = Fs * (0:(L/2)) / L;
% Plot the average magnitude spectrum
figure;
plot(f, avg_magnitude_spectrum_smoothed_dB, 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Average Magnitude Spectrum');
xlim([0, 10]);
ylim([60,165]);
% Save the matrix to a text file
%Save graph file as txt format
full_file_name = strcat(desired_path, '.txt');
disp(['Saving to: ', full_file_name]);
writematrix(avg_magnitude_spectrum_smoothed_dB, full_file_name)

