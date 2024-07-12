clc; clear; close all;

SectionNum = 7; % section number
GroupNum = 2; % Group number

%% 1. Plot the original signal x(t) for -5 < t < 5 
% Define time range 
SamplingFreq = 16000; 
SamplingPer = 1/SamplingFreq; 
time_vec = -5:SamplingPer:5; 

% Define x(t) as |t| for -1 <= t <= 1, and 0 elsewhere
sig = zeros(size(time_vec));
sig(time_vec >= -1 & time_vec <= 1) = abs(time_vec(time_vec >= -1 & time_vec <= 1));

figure;
subplot(4,1,1);
plot(time_vec, sig, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('x(t) = |t| for -1 <= t <= 1 and 0 elsewhere');
grid on;

%% 2. Plot the Time-Shifted signal shifted(t)
% shifted= x(t - 0.2 * s)
shift_amt = -0.2 * SectionNum;
shifted_sig = zeros(size(time_vec));
shifted_time_vec = time_vec + shift_amt;
shifted_sig(shifted_time_vec >= -1 & shifted_time_vec <= 1) = abs(shifted_time_vec(shifted_time_vec >= -1 & shifted_time_vec <= 1));

subplot(4,1,2);
plot(time_vec, shifted_sig, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('shifted(t)');
grid on;

%% 3. Plot the Time-Scaled signal scaled_sig(t)
scaling_factor = 1/GroupNum;
scaled_sig = zeros(size(time_vec));
scaled_time_vec = time_vec * scaling_factor + shift_amt;
scaled_sig(scaled_time_vec >= -1 & scaled_time_vec <= 1) = abs(scaled_time_vec(scaled_time_vec >= -1 & scaled_time_vec <= 1));

subplot(4,1,3);
plot(time_vec, scaled_sig, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_sig(t)');
grid on;

%% 4. Plot the Amplitude-Scaled signal scaled_Amp_sig(t)
Amp_factor = 0.1 * GroupNum * SectionNum;
scaled_Amp_sig = Amp_factor * scaled_sig;

subplot(4,1,4);
plot(time_vec, scaled_Amp_sig, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_Amp\_sig(t)');
grid on;

%% 5. Plot the filter system filter_sys(t)
filter_sys = (1/sqrt(2*pi)) * exp(-time_vec.^2 / 2);
figure;
plot(time_vec, filter_sys, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('filter\_sys(t)');
grid on;

%% 6. Compute and plot the magnitude spectrum of the input signal: |input_signal(w)|
freq_vec = linspace(-SamplingFreq/2, SamplingFreq/2, length(sig)); 
ang_freq_vec = 2 * pi * freq_vec; 
num_periods = SamplingFreq;
input_signal = fft(sig) / num_periods; 
figure;
plot(ang_freq_vec, abs(fftshift(input_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|input\_signal(w)|');
xlim([-50000 50000]);
grid on;

%% 7. Compute and plot the magnitude spectrum of the input signal: |input_Amp_sig(w)|
input_Amp_sig = fft(scaled_Amp_sig) / num_periods; 
figure;
plot(ang_freq_vec, abs(fftshift(input_Amp_sig)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|input\_Amp\_sig(w)|');
xlim([-50000 50000]);
grid on;

%% 8. Compute and plot the magnitude spectrum of the filtering system: |filter_sys(w)|
filter_sys_freq = fft(filter_sys) / num_periods; 
figure;
plot(ang_freq_vec, fftshift(abs(filter_sys_freq)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|filter\_sys(w)|');
xlim([-50000 50000]);
grid on;

%% 9. Compute and plot the magnitude spectrum of the output signal |output_signal(w)|
output_signal = input_signal .* filter_sys_freq; 
figure;
plot(ang_freq_vec, abs(fftshift(output_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|output\_signal(w)|');
xlim([-50000 50000]);
grid on;

%% 10. Compute and plot the time-domain output signal output_time(t)
output_time = (ifft(output_signal)) * num_periods;
figure;
plot(time_vec, real(output_time), 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('output\_time(t)');
grid on;

%% 11. Plot scaled_Amp_sig(t) and output_time(t) on the same graph
figure;
plot(time_vec, real(output_time), 'b', 'LineWidth', 1);
hold on;
plot(time_vec, scaled_Amp_sig, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_Amp\_sig(t) & output\_time(t)');
legend('output\_time(t)', 'scaled\_Amp\_sig(t)');
grid on;
