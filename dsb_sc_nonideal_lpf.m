
clear; clc; close all;

fm = 100;        
fc = 5000;     
fs = 100000;        
T = 0.03;         
t = 0:1/fs:T-1/fs;
N = length(t);

Am = 1; Ac = 1;

fprintf('========== NON-IDEAL LPF SIMULATION — DSB-SC ==========\n');
fprintf('Message frequency:   fm = %d Hz\n', fm);
fprintf('Carrier frequency:   fc = %d Hz\n', fc);
fprintf('Sampling frequency:  fs = %d Hz\n', fs);
fprintf('Duration:            T  = %.0f ms\n', T*1000);
fprintf('======================================================\n\n');

m_t = Am * cos(2*pi*fm*t);
c_t = Ac * cos(2*pi*fc*t);
s_t = m_t .* c_t;

f_axis = (-N/2:N/2-1) * (fs/N);

H_ideal = double(abs(f_axis) <= 1.5*fm);

v_ideal = s_t .* cos(2*pi*fc*t);
V_ideal = fftshift(fft(v_ideal));
m_ideal = real(ifft(ifftshift(V_ideal .* H_ideal)));
scale_factor = max(abs(m_ideal)) / max(abs(m_t));
m_ideal_norm = m_ideal / scale_factor;

NMSE_ideal = sum((m_t - m_ideal_norm).^2) / sum(m_t.^2);
fprintf('Ideal brick-wall LPF:  NMSE = %.2e\n\n', NMSE_ideal);

filter_orders = [2, 4, 6];
cutoff_norm = 1.5 * fm / (fs/2);  

NMSE_butter = zeros(size(filter_orders));
m_demod_butter = zeros(length(filter_orders), N);
gd_values = zeros(size(filter_orders));

fprintf('--- Butterworth Filter Sweep ---\n');
fprintf('Cutoff frequency: %.0f Hz (1.5 × fm)\n', 1.5*fm);
fprintf('Normalized cutoff: %.6f\n\n', cutoff_norm);
fprintf('%-10s %-15s %-15s %-15s\n', 'Order', 'Group Delay', 'NMSE', 'NMSE (dB)');
fprintf('%-10s %-15s %-15s %-15s\n', '-----', '-----------', '----', '--------');

for i = 1:length(filter_orders)
    n_order = filter_orders(i);
    
    [b, a] = butter(n_order, cutoff_norm);
    
    [gd_vec, f_vec] = grpdelay(b, a, 512, fs);
    [~, gd_idx] = min(abs(f_vec - fm));
    delay_samples = round(gd_vec(gd_idx));
    gd_values(i) = delay_samples;
    
    v_t = s_t .* cos(2*pi*fc*t);
    
    m_filtered = filter(b, a, v_t);
    m_compensated = [m_filtered(delay_samples+1:end), zeros(1, delay_samples)];
    
    m_compensated = m_compensated / scale_factor;
    m_demod_butter(i, :) = m_compensated;
    
    NMSE_butter(i) = sum((m_t - m_compensated).^2) / sum(m_t.^2);
    
    fprintf('%-10d %-15s %-15.4e %-15.2f\n', ...
        n_order, sprintf('%d samples', delay_samples), ...
        NMSE_butter(i), 10*log10(NMSE_butter(i) + eps));
end

figure('Position', [100 100 800 500]);
bar(filter_orders, NMSE_butter, 0.5, 'FaceColor', [0.1 0.5 0.9]);

for i = 1:length(filter_orders)
    text(filter_orders(i), NMSE_butter(i) + 0.01, ...
        sprintf('%.2e', NMSE_butter(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

hold on;
bar_highlight = bar(4, NMSE_butter(2), 0.5, 'FaceColor', [0.9 0.3 0.1]);

xlabel('Butterworth Filter Order', 'FontSize', 12);
ylabel('NMSE', 'FontSize', 12);
title('Non-Ideal LPF: Demodulation Error vs Filter Order', 'FontSize', 13);
legend('Butterworth LPF', 'Order 4 (used in combined analysis)', 'Location', 'northwest');
set(gca, 'FontSize', 11, 'XTick', filter_orders);
grid on;

fprintf('\n--- Figure 1: NMSE vs Filter Order generated ---\n');

figure('Position', [100 100 900 600]);

f_plot = linspace(0, 500, 1000);

subplot(2,1,1);
plot(f_axis(f_axis >= 0 & f_axis <= 500), ...
    H_ideal(f_axis >= 0 & f_axis <= 500), 'k-', 'LineWidth', 2);
hold on;

colors = [0.2 0.7 0.2; 0.1 0.5 0.9; 0.9 0.3 0.1];
line_styles = {'-', '--', ':'};

for i = 1:length(filter_orders)
    [b, a] = butter(filter_orders(i), cutoff_norm);
    [H_resp, f_resp] = freqz(b, a, 1024, fs);
    plot(f_resp, abs(H_resp), line_styles{i}, 'Color', colors(i,:), 'LineWidth', 2);
end

xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('|H(f)|', 'FontSize', 11);
title('Magnitude Response: Ideal vs Butterworth Filters', 'FontSize', 12);
legend('Ideal brick-wall', 'Butter order 2', 'Butter order 4', 'Butter order 6', ...
    'Location', 'northeast');
xlim([0 500]); grid on;
xline(1.5*fm, 'r--', 'Cutoff = 150 Hz', 'LineWidth', 1);

subplot(2,1,2);
for i = 1:length(filter_orders)
    [b, a] = butter(filter_orders(i), cutoff_norm);
    [gd, f_gd] = grpdelay(b, a, 512, fs);
    plot(f_gd, gd/fs*1000, line_styles{i}, 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end

xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Group Delay (ms)', 'FontSize', 11);
title('Group Delay: Higher Order = More Delay Distortion', 'FontSize', 12);
legend('Order 2', 'Order 4', 'Order 6', 'Location', 'northeast');
xlim([0 500]); grid on;

sgtitle('Butterworth Filter Characteristics', 'FontSize', 13);

fprintf('--- Figure 2: Frequency response comparison generated ---\n');

figure('Position', [100 100 900 700]);

subplot(length(filter_orders)+1, 1, 1);
plot(t*1000, m_t, 'b-', 'LineWidth', 1.5); hold on;
plot(t*1000, m_ideal_norm, 'r--', 'LineWidth', 1);
title(sprintf('Ideal Brick-Wall LPF — NMSE = %.2e', NMSE_ideal), 'FontSize', 11);
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Original m(t)', 'Demodulated', 'Location', 'northeast');
grid on; ylim([-1.5 1.5]);

for i = 1:length(filter_orders)
    subplot(length(filter_orders)+1, 1, i+1);
    plot(t*1000, m_t, 'b-', 'LineWidth', 1.5); hold on;
    plot(t*1000, m_demod_butter(i, :), 'r--', 'LineWidth', 1);
    
    title(sprintf('Butterworth Order %d — NMSE = %.4e, Group Delay = %d samples', ...
        filter_orders(i), NMSE_butter(i), gd_values(i)), 'FontSize', 11);
    xlabel('Time (ms)'); ylabel('Amplitude');
    grid on; ylim([-1.5 1.5]);
end
sgtitle('Non-Ideal LPF: Time-Domain Demodulation Comparison', 'FontSize', 13);

fprintf('--- Figure 3: Time domain comparison generated ---\n');

figure('Position', [100 100 900 500]);
for i = 1:length(filter_orders)
    subplot(length(filter_orders), 1, i);
    error_signal = m_t - m_demod_butter(i, :);
    plot(t*1000, error_signal, 'Color', [0.8 0 0], 'LineWidth', 1);
    
    title(sprintf('Error Signal — Butterworth Order %d (max error: %.4f)', ...
        filter_orders(i), max(abs(error_signal))), 'FontSize', 11);
    xlabel('Time (ms)'); ylabel('Error');
    grid on;
end
sgtitle('Non-Ideal LPF: Error Signals', 'FontSize', 13);

fprintf('--- Figure 4: Error signals generated ---\n');

fprintf('\n============================================================\n');
fprintf('               NON-IDEAL LPF SUMMARY\n');
fprintf('============================================================\n');
fprintf('Key findings:\n');
fprintf('  1. Even a 2nd-order Butterworth causes NMSE = %.2e\n', NMSE_butter(1));
fprintf('  2. Order 4 (used in combined analysis): NMSE = %.2e\n', NMSE_butter(2));
fprintf('  3. LPF distortion >> phase error >> AWGN at moderate conditions\n');
fprintf('  4. Higher order = sharper cutoff but more group delay distortion\n');
fprintf('  5. Error is DETERMINISTIC — unlike AWGN, it cannot be reduced\n');
fprintf('     by increasing transmit power.\n');
fprintf('\nIn class: ideal brick-wall LPF, |H(f)|=1 below cutoff, 0 above.\n');
fprintf('In practice: finite rolloff, passband ripple, group delay.\n');
fprintf('This is the DOMINANT impairment in the paper''s analysis.\n');
fprintf('============================================================\n');