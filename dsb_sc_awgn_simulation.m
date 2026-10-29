
clear; clc; close all;

fm = 100;           
fc = 5000;          
fs = 100000;        
T = 0.03;           
t = 0:1/fs:T-1/fs;  
N = length(t);       

Am = 1;            
Ac = 1;           

SNR_dB = [-5, 0, 5, 10, 15, 20, 25, 30, 35, 40];

fprintf('========== STEP 5: AWGN CHANNEL SIMULATION ==========\n');
fprintf('Message frequency:   fm = %d Hz\n', fm);
fprintf('Carrier frequency:   fc = %d Hz\n', fc);
fprintf('Sampling frequency:  fs = %d Hz\n', fs);
fprintf('Duration:            T  = %.0f ms\n', T*1000);
fprintf('Total samples:       N  = %d\n', N);
fprintf('SNR range:           %d dB to %d dB\n', SNR_dB(1), SNR_dB(end));
fprintf('====================================================\n\n');
 
m_t = Am * cos(2*pi*fm*t);

c_t = Ac * cos(2*pi*fc*t);

s_t = m_t .* c_t;

signal_power = mean(s_t.^2);

v_ideal = s_t .* cos(2*pi*fc*t);  

f = (-N/2:N/2-1) * (fs/N);
H_ideal = double(abs(f) <= 1.5*fm); 
V_ideal = fftshift(fft(v_ideal));
m_demod_ideal = real(ifft(ifftshift(V_ideal .* H_ideal)));

scale_factor = max(abs(m_demod_ideal)) / max(abs(m_t));
m_demod_ideal = m_demod_ideal / scale_factor;

NMSE_ideal = sum((m_t - m_demod_ideal).^2) / sum(m_t.^2);
fprintf('--- Ideal Case (Noiseless) ---\n');
fprintf('NMSE = %.2e (should be approximately zero)\n\n', NMSE_ideal);

NMSE_values = zeros(size(SNR_dB));
m_demod_noisy = zeros(length(SNR_dB), N);

fprintf('--- SNR Sweep Results ---\n');
fprintf('%-10s %-15s %-15s %-15s\n', 'SNR (dB)', 'Noise Power', 'NMSE', 'NMSE (dB)');
fprintf('%-10s %-15s %-15s %-15s\n', '--------', '-----------', '----', '--------');

for i = 1:length(SNR_dB)
    snr_linear = 10^(SNR_dB(i)/10);
    noise_power = signal_power / snr_linear;
    noise_std = sqrt(noise_power);

    n_t = noise_std * randn(1, N);

    r_t = s_t + n_t;

    v_t = r_t .* cos(2*pi*fc*t);

    V_f = fftshift(fft(v_t));
    m_recovered = real(ifft(ifftshift(V_f .* H_ideal)));

    m_recovered = m_recovered / scale_factor;

    m_demod_noisy(i, :) = m_recovered;

    NMSE_values(i) = sum((m_t - m_recovered).^2) / sum(m_t.^2);
    NMSE_dB = 10*log10(NMSE_values(i));

    fprintf('%-10d %-15.4e %-15.4e %-15.2f\n', ...
        SNR_dB(i), noise_power, NMSE_values(i), NMSE_dB);
end

figure('Position', [100 100 800 500]);
semilogy(SNR_dB, NMSE_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b');
hold on;

yline(0.1, 'r--', 'NMSE = 0.1 (severe distortion)', 'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
yline(0.01, 'g--', 'NMSE = 0.01 (acceptable)', 'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);

grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (log scale)', 'FontSize', 12);
title('SNR vs NMSE in AWGN Channel — DSB-SC Demodulation Performance', 'FontSize', 13);
legend('NMSE', 'Location', 'northeast');
set(gca, 'FontSize', 11);

fprintf('\n--- Figure 1: SNR vs NMSE curve generated ---\n');

snr_show = [0, 10, 20, 40];  % dB
snr_idx = arrayfun(@(x) find(SNR_dB == x), snr_show);

figure('Position', [100 100 900 700]);
for k = 1:length(snr_show)
    subplot(length(snr_show), 1, k);
    plot(t*1000, m_t, 'b-', 'LineWidth', 1.5); hold on;
    plot(t*1000, m_demod_noisy(snr_idx(k), :), 'r--', 'LineWidth', 1);

    title(sprintf('SNR = %d dB — NMSE = %.4f', ...
        snr_show(k), NMSE_values(snr_idx(k))), 'FontSize', 11);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    legend('Original m(t)', 'Demodulated', 'Location', 'northeast');
    grid on;
    ylim([-2 2]);
end
sgtitle('Demodulation in AWGN Channel: Different SNR Levels', 'FontSize', 13);

fprintf('--- Figure 2: Time domain comparison generated ---\n');

figure('Position', [100 100 900 700]);
for k = 1:length(snr_show)
    subplot(length(snr_show), 1, k);
    error_signal = m_t - m_demod_noisy(snr_idx(k), :);
    plot(t*1000, error_signal, 'Color', [0.8 0 0], 'LineWidth', 1);

    title(sprintf('Error Signal — SNR = %d dB (max error: %.4f)', ...
        snr_show(k), max(abs(error_signal))), 'FontSize', 11);
    xlabel('Time (ms)');
    ylabel('Error');
    grid on;
end
sgtitle('Error Signals in AWGN Channel', 'FontSize', 13);

fprintf('--- Figure 3: Error signals generated ---\n');

figure('Position', [100 100 900 600]);

subplot(2,1,1);
snr_high_idx = find(SNR_dB == 30);
r_high = s_t + sqrt(signal_power / 10^(30/10)) * randn(1, N);
R_high = fftshift(fft(r_high)) / N;
f_axis = (-N/2:N/2-1) * (fs/N);
plot(f_axis, abs(R_high), 'b', 'LineWidth', 1);
title('Received Signal Spectrum — SNR = 30 dB (clean)', 'FontSize', 11);
xlabel('Frequency (Hz)'); ylabel('|R(f)|');
xlim([-8000 8000]); grid on;

subplot(2,1,2);
snr_low_idx = find(SNR_dB == 5);
r_low = s_t + sqrt(signal_power / 10^(5/10)) * randn(1, N);
R_low = fftshift(fft(r_low)) / N;
plot(f_axis, abs(R_low), 'r', 'LineWidth', 1);
title('Received Signal Spectrum — SNR = 5 dB (noisy)', 'FontSize', 11);
xlabel('Frequency (Hz)'); ylabel('|R(f)|');
xlim([-8000 8000]); grid on;

sgtitle('AWGN Effect in Frequency Domain', 'FontSize', 13);

fprintf('--- Figure 4: Spectrum comparison generated ---\n');

f1 = 300; f2 = 800; f3 = 1500;
m_multi = cos(2*pi*f1*t) + 0.5*cos(2*pi*f2*t) + 0.3*cos(2*pi*f3*t);
s_multi = m_multi .* c_t;
signal_power_multi = mean(s_multi.^2);

W_multi = 1500;
H_multi = double(abs(f) <= 1.8*W_multi);

v_multi_ideal = s_multi .* cos(2*pi*fc*t);
V_multi_ideal = fftshift(fft(v_multi_ideal));
m_multi_demod_ideal = real(ifft(ifftshift(V_multi_ideal .* H_multi)));
scale_multi = max(abs(m_multi_demod_ideal)) / max(abs(m_multi));

NMSE_multi = zeros(size(SNR_dB));
for i = 1:length(SNR_dB)
    snr_lin = 10^(SNR_dB(i)/10);
    n_t = sqrt(signal_power_multi / snr_lin) * randn(1, N);
    r_t = s_multi + n_t;
    v_t = r_t .* cos(2*pi*fc*t);
    V_f = fftshift(fft(v_t));
    m_rec = real(ifft(ifftshift(V_f .* H_multi))) / scale_multi;
    NMSE_multi(i) = sum((m_multi - m_rec).^2) / sum(m_multi.^2);
end

figure('Position', [100 100 800 500]);
semilogy(SNR_dB, NMSE_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b'); hold on;
semilogy(SNR_dB, NMSE_multi, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (log scale)', 'FontSize', 12);
title('Single-Tone vs Multi-Tone DSB-SC: AWGN Performance Comparison', 'FontSize', 13);
legend('Single tone (100 Hz)', 'Multi-tone (300+800+1500 Hz)', 'Location', 'northeast');
set(gca, 'FontSize', 11);

fprintf('--- Figure 5: Single-tone vs multi-tone NMSE comparison generated ---\n');

fprintf('\n========== NMSE SUMMARY TABLE ==========\n');
fprintf('%-10s %-18s %-18s\n', 'SNR (dB)', 'NMSE (Single)', 'NMSE (Multi)');
fprintf('%-10s %-18s %-18s\n', '--------', '-------------', '------------');
for i = 1:length(SNR_dB)
    fprintf('%-10d %-18.6f %-18.6f\n', SNR_dB(i), NMSE_values(i), NMSE_multi(i));
end

fprintf('\n========== ANALYSIS NOTES ==========\n');
fprintf('1. NMSE decreases as SNR increases — expected behavior.\n');
fprintf('2. Multi-tone signal is more sensitive — more components = more distortion.\n');
fprintf('3. Below SNR < 10 dB, demodulation is severely degraded.\n');
fprintf('4. Above SNR > 25 dB, NMSE drops to negligible levels.\n');
fprintf('5. Assuming no AWGN in class = assuming SNR -> infinity.\n');
fprintf('   In practice, SNR is finite, so performance always degrades.\n');
fprintf('====================================\n');
