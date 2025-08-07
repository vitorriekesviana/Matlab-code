% CONTROL SCRIPT: Visualize frequency content of simulated sine waves
% Author: Vitor Riekes Viana
% ---------------------------------------
% This script simulates 60 sine waves (each with a unique frequency between 1–60 Hz)
% and allows the user to visualize selected signals in time and frequency domains.


% Ask the user how many frequency signals they want to visualize
visualize = input('How many frequency signals do you want to visualize? (1-60): ');

% Validate input
if visualize < 1 || visualize > 60
    error('Please enter a value between 1 and 60.');
end

% Define parameters
num_channels = 60;          % Total number of available signals
num_samples = 60;           % Number of samples per signal
Fs = 200;                   % Sampling frequency (Hz)
T = 1/Fs;                   % Sampling period (s)
L = num_samples;            % Length of each signal
t = (0:L-1)*T;              % Time vector

% Frequency range (example: from 1 Hz to 60 Hz)
frequency_range = linspace(1, 60, num_channels);

% Preallocate matrix to store sine waves
matrix_control = zeros(num_channels, num_samples);

% Generate sine waves with different frequencies for each channel
for i = 1:num_channels
    f = frequency_range(i);
    matrix_control(i, :) = sin(2 * pi * f * t);
end

% Display the frequency of the selected channels
disp('Frequencies of selected channels:');
disp(frequency_range(1:visualize));

% Display the generated signal matrix (optional)
% disp(matrix_control);

%% Plot magnitude spectrum of selected channels
figure;
hold on;
n_colors = visualize;
colors = hsv(n_colors); % Generate distinct colors

for i = 1:visualize
    % Perform FFT
    X = fft(matrix_control(i, :));
    P2 = abs(X/L);            % Two-sided spectrum
    P1 = P2(1:L/2+1);         % Single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);

    % Frequency domain
    f = Fs*(0:(L/2))/L;

    % Plot magnitude spectrum
    plot(f, P1, 'Color', colors(i, :), 'DisplayName', ...
        ['Channel ' num2str(i) ' - ' num2str(frequency_range(i), '%.2f') ' Hz']);
end

xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(['Magnitude Spectrum of First ' num2str(visualize) ' Channels']);
legend;
xlim([0, 80]);
hold off;

%% Compute and plot average magnitude spectrum
magnitude_spectra = zeros(visualize, L/2+1);
for i = 1:visualize
    X = fft(matrix_control(i, :));
    P2 = abs(X/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    magnitude_spectra(i, :) = P1;
end

% Average magnitude spectrum (in dB)
avg_magnitude_spectrum = mean(magnitude_spectra, 1);
avg_magnitude_spectrum_dB = 20 * log10(avg_magnitude_spectrum);

% Frequency domain for plot
f = Fs*(0:(L/2))/L;

figure;
plot(f, avg_magnitude_spectrum_dB, 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(['Average Magnitude Spectrum of First ' num2str(visualize) ' Channels']);
xlim([0, 80]);

%% Plot the time-domain sine waves
figure;
hold on;
for i = 1:visualize
    plot(t, matrix_control(i, :), 'Color', colors(i, :), ...
        'DisplayName', ['Channel ' num2str(i) ' - ' num2str(frequency_range(i), '%.2f') ' Hz']);
end
xlabel('Time (s)');
ylabel('Amplitude');
title(['Time-Domain Sine Waves of First ' num2str(visualize) ' Channels']);
legend;
hold off;

