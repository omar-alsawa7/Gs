clc; clear; close all;

SectionNum = 3; % section number
GroupNum = 1; % Group number

%% 1. Plot the original signal signal_t(t) for -5 < t < 5 
% Define time range 
SamplingFreq = 16000; 
SamplingPeriod = 1/SamplingFreq; 
time_vec = -5:SamplingPeriod:5; 

% Define signal_t(t) as |t| for -1 <= t <= 1, and 0 elsewhere
signal_t = zeros(size(time_vec));
signal_t(time_vec >= -1 & time_vec <= 1) = abs(time_vec(time_vec >= -1 & time_vec <= 1));

figure;
subplot(4,1,1);
plot(time_vec, signal_t, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('signal\_t(t) = |t| for -1 <= t <= 1 and 0 elsewhere');
grid on;

%% 2. Plot the Time-Shifted signal shifted_signal(t)
% shifted_signal= signal_t(t - 0.2 * SectionNum)
shift_amt = -0.2 * SectionNum;
shifted_signal = zeros(size(time_vec));
shifted_time = time_vec + shift_amt;
shifted_signal(shifted_time >= -1 & shifted_time <= 1) = abs(shifted_time(shifted_time >= -1 & shifted_time <= 1));

subplot(4,1,2);
plot(time_vec, shifted_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('shifted\_signal(t)');
grid on;

%% 3. Plot the Time-Scaled signal scaled_signal(t)
scaling_factor = GroupNum + 1;
scaled_signal = zeros(size(time_vec));
scaled_time = time_vec * scaling_factor + shift_amt;
scaled_signal(scaled_time >= -1 & scaled_time <= 1) = abs(scaled_time(scaled_time >= -1 & scaled_time <= 1));

subplot(4,1,3);
plot(time_vec, scaled_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_signal(t)');
grid on;

%% 4. Plot the Amplitude-Scaled signal scaled_amp_signal(t)
Amp_factor = 0.1 * GroupNum * SectionNum;
scaled_amp_signal = Amp_factor * scaled_signal;

subplot(4,1,4);
plot(time_vec, scaled_amp_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_amp\_signal(t)');
grid on;

%% 5. Plot the filter system filter_sys(t)
filter_sys = (1/sqrt(2*pi)) * exp(-time_vec.^2 / 2);
figure;
plot(time_vec, filter_sys, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('filter\_sys(t)');
grid on;

%% 6. Compute and plot the magnitude spectrum of the input signal: |Input_signal(w)|
freq_vec = linspace(-SamplingFreq/2, SamplingFreq/2, length(signal_t)); 
ang_freq_vec = 2 * pi * freq_vec; 
num_periods = SamplingFreq;
Input_signal = fft(signal_t) / num_periods; 
figure;
plot(ang_freq_vec, abs(fftshift(Input_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Input\_signal(w)|');
xlim([-50000 50000]);
grid on;

%% 7. Compute and plot the magnitude spectrum of the input signal: |Amp_scaled_signal(w)|
Amp_scaled_signal = fft(scaled_amp_signal) / num_periods; 
figure;
plot(ang_freq_vec, abs(fftshift(Amp_scaled_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Amp\_scaled\_signal(w)|');
xlim([-50000 50000]);
grid on;

%% 8. Compute and plot the magnitude spectrum of the filtering system: |Filter_sys(w)|
Filter_sys = fft(filter_sys) / num_periods; 
figure;
plot(ang_freq_vec, fftshift(abs(Filter_sys)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Filter\_sys(w)|');
xlim([-50000 50000]);
grid on;

%% 9. Compute and plot the magnitude spectrum of the output signal |Output_signal(w)|
Output_signal = Input_signal .* Filter_sys; 
figure;
plot(ang_freq_vec, abs(fftshift(Output_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Output\_signal(w)|');
xlim([-50000 50000]);
grid on;

%% 10. Compute and plot the time-domain output signal output_time(t)
output_time = (ifft(Output_signal)) * num_periods;
figure;
plot(time_vec, real(output_time), 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('output\_time(t)');
grid on;

%% 11. Plot scaled_amp_signal(t) and output_time(t) on the same graph
figure;
plot(time_vec, real(output_time), 'b', 'LineWidth', 1);
hold on;
plot(time_vec, scaled_amp_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_amp\_signal(t) & output\_time(t)');
legend('output\_time(t)', 'scaled\_amp\_signal(t)');
grid on;
