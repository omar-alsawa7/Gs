clc; clear; close all;

SectionNum = 3; 
GroupNum = 1;

SamplingFreq = 16000; 
SamplingPeriod = 1/SamplingFreq; 
time_vec = -5:SamplingPeriod:5; 

signal_t = zeros(size(time_vec));
signal_t(time_vec >= -1 & time_vec <= 1) = abs(time_vec(time_vec >= -1 & time_vec <= 1));

figure;
subplot(4,1,1);
plot(time_vec, signal_t, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('signal\_t(t) = |t| for -1 <= t <= 1 and 0 elsewhere');
grid on;

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

Amp_factor = 0.1 * GroupNum * SectionNum;
scaled_amp_signal = Amp_factor * scaled_signal;

subplot(4,1,4);
plot(time_vec, scaled_amp_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_amp\_signal(t)');
grid on;

filter_sys = (1/sqrt(2*pi)) * exp(-time_vec.^2 / 2);
figure;
plot(time_vec, filter_sys, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('filter\_sys(t)');
grid on;

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

Amp_scaled_signal = fft(scaled_amp_signal) / num_periods; 
figure;
plot(ang_freq_vec, abs(fftshift(Amp_scaled_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Amp\_scaled\_signal(w)|');
xlim([-50000 50000]);
grid on;

Filter_sys = fft(filter_sys) / num_periods; 
figure;
plot(ang_freq_vec, fftshift(abs(Filter_sys)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Filter\_sys(w)|');
xlim([-50000 50000]);
grid on;

Output_signal = Input_signal .* Filter_sys; 
figure;
plot(ang_freq_vec, abs(fftshift(Output_signal)), 'r', 'LineWidth', 1);
xlabel('w rad/sec');
ylabel('Magnitude');
title('|Output\_signal(w)|');
xlim([-50000 50000]);
grid on;

output_time = (ifft(Output_signal)) * num_periods;
figure;
plot(time_vec, real(output_time), 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('output\_time(t)');
grid on;

figure;
plot(time_vec, real(output_time), 'b', 'LineWidth', 1);
hold on;
plot(time_vec, scaled_amp_signal, 'r', 'LineWidth', 1);
xlabel('t(sec)');
ylabel('Amplitude');
title('scaled\_amp\_signal(t) & output\_time(t)');
legend('output\_time(t)', 'scaled\_amp\_signal(t)');
grid on;
