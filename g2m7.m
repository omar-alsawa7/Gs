clc; clear; close all;

SectionNumber = 7; % section number
GroupNumber = 2; % Group number

%% 1. Plot the original signal x(t) for -5 < t < 5 
% Define time range 
SamplingFreq = 16000; 
SamplingPeriod = 1/SamplingFreq; 
time = -5:SamplingPeriod:5; 

% Define x(t) as |t| for -1 <= t <= 1, and 0 elsewhere
signal = zeros(size(time));
signal(time >= -1 & time <= 1) = abs(time(time >= -1 & time <= 1));

figure;
subplot(4,1,1);
plot(time, signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('x(t) = |t| for -1 <= t <= 1 and 0 elsewhere');
grid on;

%% 2. Plot the Time-Shifted signal shifted(t)
% shifted = x(t - 0.2 * s)
shift = -0.2 * SectionNumber;
shifted = zeros(size(time));
shifted_time = time + shift;
shifted(shifted_time >= -1 & shifted_time <= 1) = abs(shifted_time(shifted_time >= -1 & shifted_time <= 1));

subplot(4,1,2);
plot(time, shifted, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('shifted(t)');
grid on;

%% 3. Plot the Time-Scaled signal scaled(t)
alpha = 1/GroupNumber;
scaled = zeros(size(time));
scaled_time = time * alpha + shift;
scaled(scaled_time >= -1 & scaled_time <= 1) = abs(scaled_time(scaled_time >= -1 & scaled_time <= 1));

subplot(4,1,3);
plot(time, scaled, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled(t)');
grid on;

%% 4. Plot the Amplitude-Scaled signal scaledAmplitude(t)
Amplitude = 0.1 * GroupNumber * SectionNumber;
scaledAmplitude = Amplitude * scaled;

subplot(4,1,4);
plot(time, scaledAmplitude, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaledAmplitude(t)');
grid on;

%% 5. Plot the filter system h(t)
filterSystem = (1/sqrt(2*pi)) * exp(-time.^2 / 2);
figure;
plot(time, filterSystem, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('h(t)');
grid on;

%% 6. Compute and plot the magnitude spectrum of the input signal: |inputSignal(w)|
frequency = linspace(-SamplingFreq/2, SamplingFreq/2, length(signal)); 
angularVelocity = 2 * pi * frequency; 
numberOfPeriods = SamplingFreq;
inputSignal = fft(signal) / numberOfPeriods; 
figure;
plot(angularVelocity, abs(fftshift(inputSignal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|inputSignal(w)|');
xlim([-50000 50000]);
grid on;

%% 7. Compute and plot the magnitude spectrum of the input signal: |inputSignalAmplitude(w)|
inputSignalAmplitude = fft(scaledAmplitude) / numberOfPeriods; 
figure;
plot(angularVelocity, abs(fftshift(inputSignalAmplitude)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|inputSignalAmplitude(w)|');
xlim([-50000 50000]);
grid on;

%% 8. Compute and plot the magnitude spectrum of the filtering system: |filterSystem(w)|
filterSystemFrequency = fft(filterSystem) / numberOfPeriods; 
figure;
plot(angularVelocity, fftshift(abs(filterSystemFrequency)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|filterSystem(w)|');
xlim([-50000 50000]);
grid on;

%% 9. Compute and plot the magnitude spectrum of the output signal |outputSignal(w)|
outputSignal = inputSignal .* filterSystemFrequency; 
figure;
plot(angularVelocity, abs(fftshift(outputSignal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|outputSignal(w)|');
xlim([-50000 50000]);
grid on;

%% 10. Compute and plot the time-domain output signal outputTime(t)
outputTime = (ifft(outputSignal)) * numberOfPeriods;
figure;
plot(time, real(outputTime), 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('outputTime(t)');
grid on;

%% 11. Plot scaledAmplitude(t) and outputTime(t) on the same graph
figure;
plot(time, real(outputTime), 'b', 'LineWidth', 1);
hold on;
plot(time, scaledAmplitude, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaledAmplitude(t) & outputTime(t)');
legend('outputTime(t)', 'scaledAmplitude(t)');
grid on;
