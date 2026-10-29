clear; clc; close all;

fm = 100;          
fc = 5000;        
fs = 100000;       
T = 0.03;           
t = 0:1/fs:T-1/fs;
N = length(t);

Am = 1; Ac = 1;

fprintf('========== PHASE ERROR SIMULATION — DSB-SC ==========\n');
fprintf('Message frequency:   fm = %d Hz\n', fm);
fprintf('Carrier frequency:   fc = %d Hz\n', fc);
fprintf('Sampling frequency:  fs = %d Hz\n', fs);
fprintf('Duration:            T  = %.0f ms\n', T*1000);
fprintf('====================================================\n\n');

m_t = Am * cos(2*pi*fm*t);
c_t = Ac * cos(2*pi*fc*t);
s_t = m_t .* c_t;

f_axis = (-N/2:N/2-1) * (fs/N);
H_ideal = double(abs(f_axis) <= 1.5*fm);

v_ideal = s_t .* cos(2*pi*fc*t);
V_ideal = fftshift(fft(v_ideal));
m_ideal = real(ifft(ifftshift(V_ideal .* H_ideal)));
scale_factor = max(abs(m_ideal)) / max(abs(m_t));

NMSE_ideal = sum((m_t - m_ideal/scale_factor).^2) / sum(m_t.^2);
fprintf('Ideal case (dphi = 0):  NMSE = %.2e\n\n', NMSE_ideal);

phi_deg = 0:5:90;     
phi_rad = phi_deg * pi / 180;
NMSE_phase = zeros(size(phi_deg));
m_demod_phase = zeros(length(phi_deg), N);

fprintf('--- Phase Error Sweep ---\n');
fprintf('%-15s %-15s %-15s %-15s\n', 'Phase (deg)', 'cos(dphi)', 'NMSE', 'NMSE (dB)');
fprintf('%-15s %-15s %-15s %-15s\n', '-----------', '---------', '----', '--------');

for i = 1:length(phi_deg)
    lo_t = cos(2*pi*fc*t + phi_rad(i));
    
    v_t = s_t .* lo_t;
    
    V_f = fftshift(fft(v_t));
    m_recovered = real(ifft(ifftshift(V_f .* H_ideal)));
    
    m_recovered = m_recovered / scale_factor;
    m_demod_phase(i, :) = m_recovered;
    
    NMSE_phase(i) = sum((m_t - m_recovered).^2) / sum(m_t.^2);
    
    fprintf('%-15d %-15.4f %-15.4e %-15.2f\n', ...
        phi_deg(i), cos(phi_rad(i)), NMSE_phase(i), 10*log10(NMSE_phase(i) + eps));
end

NMSE_theory = (1 - cos(phi_rad)).^2;

fprintf('\n--- Theory vs Simulation ---\n');
fprintf('%-15s %-15s %-15s %-15s\n', 'Phase (deg)', 'NMSE (sim)', 'NMSE (theory)', 'Match?');
fprintf('%-15s %-15s %-15s %-15s\n', '-----------', '----------', '-------------', '------');
for i = 1:length(phi_deg)
    match = abs(NMSE_phase(i) - NMSE_theory(i)) / max(NMSE_theory(i), eps) < 0.05;
    fprintf('%-15d %-15.4e %-15.4e %-15s\n', ...
        phi_deg(i), NMSE_phase(i), NMSE_theory(i), string(match));
end

figure('Position', [100 100 800 500]);
plot(phi_deg, NMSE_phase, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b'); hold on;
plot(phi_deg, NMSE_theory, 'r--', 'LineWidth', 1.5);

idx_30 = find(phi_deg == 30);
if ~isempty(idx_30)
    plot(30, NMSE_phase(idx_30), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    text(32, NMSE_phase(idx_30)*1.3, sprintf('30°: NMSE = %.2e', NMSE_phase(idx_30)), ...
        'FontSize', 10, 'Color', 'r');
end

grid on;
xlabel('Phase Error \Delta\phi (degrees)', 'FontSize', 12);
ylabel('NMSE', 'FontSize', 12);
title('Phase Synchronization Error vs Demodulation Quality — DSB-SC', 'FontSize', 13);
legend('Simulation', 'Theory: (1 - cos\Delta\phi)^2', '30° (used in combined analysis)', ...
    'Location', 'northwest');
set(gca, 'FontSize', 11);

fprintf('\n--- Figure 1: NMSE vs Phase Error generated ---\n');

phi_show = [0, 15, 30, 45, 60, 90];
phi_show_idx = arrayfun(@(x) find(phi_deg == x), phi_show);

figure('Position', [100 100 900 800]);
for k = 1:length(phi_show)
    subplot(length(phi_show), 1, k);
    plot(t*1000, m_t, 'b-', 'LineWidth', 1.5); hold on;
    plot(t*1000, m_demod_phase(phi_show_idx(k), :), 'r--', 'LineWidth', 1);
    
    title(sprintf('\\Delta\\phi = %d° — NMSE = %.4e, cos(\\Delta\\phi) = %.3f', ...
        phi_show(k), NMSE_phase(phi_show_idx(k)), cos(phi_show(k)*pi/180)), ...
        'FontSize', 10);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    if k == 1
        legend('Original m(t)', 'Demodulated', 'Location', 'northeast');
    end
    grid on;
    ylim([-1.5 1.5]);
end
sgtitle('Effect of Phase Error on DSB-SC Demodulation', 'FontSize', 13);

fprintf('--- Figure 2: Time domain comparison generated ---\n');

phi_detail = [15, 30, 60];
phi_detail_idx = arrayfun(@(x) find(phi_deg == x), phi_detail);

figure('Position', [100 100 900 500]);
for k = 1:length(phi_detail)
    subplot(length(phi_detail), 1, k);
    error_signal = m_t - m_demod_phase(phi_detail_idx(k), :);
    plot(t*1000, error_signal, 'Color', [0.8 0 0], 'LineWidth', 1);
    
    title(sprintf('Error Signal — \\Delta\\phi = %d° (max error: %.4f)', ...
        phi_detail(k), max(abs(error_signal))), 'FontSize', 11);
    xlabel('Time (ms)');
    ylabel('Error');
    grid on;
end
sgtitle('Phase Error: Error Signals at Different Angles', 'FontSize', 13);

fprintf('--- Figure 3: Error signals generated ---\n');

fprintf('\n============================================================\n');
fprintf('                   PHASE ERROR SUMMARY\n');
fprintf('============================================================\n');
fprintf('Key finding: Phase error causes amplitude scaling by cos(dphi).\n');
fprintf('  - At 30°: NMSE = %.2e — moderate degradation\n', NMSE_phase(idx_30));
fprintf('  - At 60°: output amplitude halved (cos60° = 0.5)\n');
fprintf('  - At 90°: complete signal loss (cos90° = 0)\n');
fprintf('  - Simulation matches theory: NMSE = (1 - cos(dphi))^2\n');
fprintf('\nIn class, we assume dphi = 0. In practice, PLLs achieve\n');
fprintf('lock with residual phase error — this script shows the cost.\n');
fprintf('============================================================\n');