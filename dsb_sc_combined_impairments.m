
clear; clc; close all;

fm = 100;           
fc = 5000;          
fs = 100000;       
T = 0.03;           
t = 0:1/fs:T-1/fs;
N = length(t);

Am = 1; Ac = 1;

m_t = Am * cos(2*pi*fm*t);
c_t = Ac * cos(2*pi*fc*t);

s_t = m_t .* c_t;
signal_power = mean(s_t.^2);

f_axis = (-N/2:N/2-1) * (fs/N);

fprintf('============================================================\n');
fprintf('   COMBINED IMPAIRMENTS ANALYSIS — DSB-SC DEMODULATION\n');
fprintf('============================================================\n');
fprintf('Parameters: fm=%d Hz, fc=%d Hz, fs=%d Hz, T=%.0f ms\n\n', fm, fc, fs, T*1000);

delta_phi_deg = 30;
delta_phi = delta_phi_deg * pi / 180;  

butter_order = 4;
cutoff_norm = 1.5 * fm / (fs/2);  
[b_butt, a_butt] = butter(butter_order, cutoff_norm);

SNR_dB_combined = 15;
snr_lin = 10^(SNR_dB_combined/10);
noise_std = sqrt(signal_power / snr_lin);

H_ideal = double(abs(f_axis) <= 1.5*fm);

[gd_vec, f_vec] = grpdelay(b_butt, a_butt, 512, fs);
[~, gd_idx] = min(abs(f_vec - fm));
delay_samples = round(gd_vec(gd_idx));

fprintf('Impairment levels:\n');
fprintf('  I1 — Phase error:     delta_phi = %d degrees\n', delta_phi_deg);
fprintf('  I2 — Non-ideal LPF:   Butterworth order %d, fc = %.0f Hz\n', butter_order, 1.5*fm);
fprintf('       Group delay at fm = %d samples (%.2f ms)\n', delay_samples, delay_samples/fs*1000);
fprintf('  I3 — AWGN channel:    SNR = %d dB\n\n', SNR_dB_combined);

v_ideal = s_t .* cos(2*pi*fc*t);    
V_ideal = fftshift(fft(v_ideal));
m_ideal = real(ifft(ifftshift(V_ideal .* H_ideal)));  
scale_factor = max(abs(m_ideal)) / max(abs(m_t));

fprintf('Scale factor (from ideal case): %.6f\n', scale_factor);
fprintf('  (All scenarios use this same scale factor)\n\n');

scenario_names = {
    'Ideal (no impairments)',           ...  
    'Phase error only',                 ...  
    'Non-ideal LPF only',              ...  
    'AWGN only',                        ...  
    'Phase error + LPF',               ...  
    'Phase error + AWGN',              ... 
    'LPF + AWGN',                       ...  
    'All three combined'                ...  
};

% [use_phase, use_butter, use_awgn]
scenario_flags = [
    0 0 0;   
    1 0 0;   
    0 1 0;   
    0 0 1;   
    1 1 0;   
    1 0 1;   
    0 1 1;   
    1 1 1;   
];

n_scenarios = size(scenario_flags, 1);
NMSE_scenarios = zeros(1, n_scenarios);
m_demod_scenarios = zeros(n_scenarios, N);

fprintf('--- Running %d scenarios ---\n\n', n_scenarios);

for s = 1:n_scenarios
    [m_hat, nmse_val] = run_demod(s_t, m_t, t, fc, delta_phi, ...
        scenario_flags(s,1), H_ideal, f_axis, fs, N, ...
        b_butt, a_butt, scenario_flags(s,2), ...
        noise_std, scenario_flags(s,3), ...
        scale_factor, delay_samples);

    NMSE_scenarios(s) = nmse_val;
    m_demod_scenarios(s, :) = m_hat;

    fprintf('  Scenario %d: %-30s  NMSE = %.6e  (%.2f dB)\n', ...
        s, scenario_names{s}, nmse_val, 10*log10(nmse_val + eps));
end

figure('Position', [100 100 1000 550]);

bar_colors = [
    0.2 0.7 0.2;   
    0.9 0.5 0.1;   
    0.1 0.5 0.9;   
    0.8 0.2 0.2;    
    0.7 0.4 0.7;   
    0.9 0.6 0.3; 
    0.3 0.6 0.7;   
    0.1 0.1 0.1;  
];

b = bar(1:n_scenarios, log10(NMSE_scenarios + eps), 'FaceColor', 'flat');
for k = 1:n_scenarios
    b.CData(k,:) = bar_colors(k,:);
end

set(gca, 'XTickLabel', scenario_names, 'XTickLabelRotation', 35, 'FontSize', 9);
ylabel('log_{10}(NMSE)', 'FontSize', 12);
title('Combined Impairments: NMSE Comparison Across All Scenarios', 'FontSize', 13);
grid on;

for k = 1:n_scenarios
    text(k, log10(NMSE_scenarios(k) + eps) + 0.15, ...
        sprintf('%.1e', NMSE_scenarios(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end

fprintf('\n--- Figure 1: NMSE bar chart generated ---\n');

show_idx = [1, 2, 7, 8]; 
show_labels = scenario_names(show_idx);

figure('Position', [100 100 900 700]);
for k = 1:length(show_idx)
    subplot(length(show_idx), 1, k);
    plot(t*1000, m_t, 'b-', 'LineWidth', 1.5); hold on;
    plot(t*1000, m_demod_scenarios(show_idx(k), :), 'r--', 'LineWidth', 1);

    title(sprintf('%s — NMSE = %.4e', show_labels{k}, NMSE_scenarios(show_idx(k))), ...
        'FontSize', 11);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    legend('Original m(t)', 'Demodulated', 'Location', 'northeast');
    grid on;
    ylim([-1.5 1.5]);
end
sgtitle('Combined Impairments: Time-Domain Comparison', 'FontSize', 13);

fprintf('--- Figure 2: Time domain comparison generated ---\n');

nmse_phase = NMSE_scenarios(2);
nmse_lpf   = NMSE_scenarios(3);
nmse_awgn  = NMSE_scenarios(4);

predicted_phase_lpf  = nmse_phase + nmse_lpf;
predicted_phase_awgn = nmse_phase + nmse_awgn;
predicted_lpf_awgn   = nmse_lpf + nmse_awgn;
predicted_all        = nmse_phase + nmse_lpf + nmse_awgn;

actual_phase_lpf  = NMSE_scenarios(5);
actual_phase_awgn = NMSE_scenarios(6);
actual_lpf_awgn   = NMSE_scenarios(7);
actual_all        = NMSE_scenarios(8);

combo_names = {'Phase+LPF', 'Phase+AWGN', 'LPF+AWGN', 'All Three'};
predicted_vals = [predicted_phase_lpf, predicted_phase_awgn, predicted_lpf_awgn, predicted_all];
actual_vals    = [actual_phase_lpf, actual_phase_awgn, actual_lpf_awgn, actual_all];

figure('Position', [100 100 900 500]);
x = 1:4;
bar_data = [predicted_vals; actual_vals]';
bh = bar(x, bar_data, 'grouped');
bh(1).FaceColor = [0.3 0.6 0.9];  
bh(2).FaceColor = [0.9 0.3 0.3];  

set(gca, 'XTickLabel', combo_names, 'FontSize', 10);
ylabel('NMSE', 'FontSize', 12);
title('Additivity Test: Predicted (Sum of Individuals) vs Actual Combined NMSE', 'FontSize', 13);
legend('Predicted (sum)', 'Actual (measured)', 'Location', 'northwest');
grid on;

for k = 1:4
    ratio = actual_vals(k) / predicted_vals(k);
    y_pos = max(predicted_vals(k), actual_vals(k));
    text(k, y_pos * 1.05, sprintf('Ratio: %.2f', ratio), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

fprintf('--- Figure 3: Additivity analysis generated ---\n');

figure('Position', [100 100 800 500]);

individual_nmse = [nmse_phase, nmse_lpf, nmse_awgn];
individual_labels = {'Phase Error (30°)', 'Non-ideal LPF (Butter-4)', 'AWGN (15 dB)'};
total_individual = sum(individual_nmse);
contributions_pct = (individual_nmse / total_individual) * 100;

pie_h = pie(individual_nmse, {
    sprintf('Phase Error\n%.2e (%.1f%%)', nmse_phase, contributions_pct(1)), ...
    sprintf('Non-ideal LPF\n%.2e (%.1f%%)', nmse_lpf, contributions_pct(2)), ...
    sprintf('AWGN\n%.2e (%.1f%%)', nmse_awgn, contributions_pct(3))});

for k = 2:2:length(pie_h)
    pie_h(k).FontSize = 8;
end

colormap([0.9 0.5 0.1; 0.1 0.5 0.9; 0.8 0.2 0.2]);

th = title(sprintf('Individual Impairment Contributions to Total NMSE\n(Sum = %.2e, Actual combined = %.2e)', ...
    total_individual, actual_all), 'FontSize', 12);
set(th, 'Units', 'normalized', 'Position', [0.5 1.12 0]);

fprintf('--- Figure 4: Dominance pie chart generated ---\n');

SNR_sweep = [0, 5, 10, 15, 20, 25, 30];
n_snr = length(SNR_sweep);

NMSE_awgn_only  = zeros(1, n_snr);
NMSE_phase_awgn = zeros(1, n_snr);
NMSE_lpf_awgn   = zeros(1, n_snr);
NMSE_all_three  = zeros(1, n_snr);

for i = 1:n_snr
    snr_i = 10^(SNR_sweep(i)/10);
    ns_i = sqrt(signal_power / snr_i);

    [~, NMSE_awgn_only(i)] = run_demod(s_t, m_t, t, fc, 0, ...
        false, H_ideal, f_axis, fs, N, b_butt, a_butt, false, ns_i, true, ...
        scale_factor, delay_samples);

    [~, NMSE_phase_awgn(i)] = run_demod(s_t, m_t, t, fc, delta_phi, ...
        true, H_ideal, f_axis, fs, N, b_butt, a_butt, false, ns_i, true, ...
        scale_factor, delay_samples);

    [~, NMSE_lpf_awgn(i)] = run_demod(s_t, m_t, t, fc, 0, ...
        false, H_ideal, f_axis, fs, N, b_butt, a_butt, true, ns_i, true, ...
        scale_factor, delay_samples);

    [~, NMSE_all_three(i)] = run_demod(s_t, m_t, t, fc, delta_phi, ...
        true, H_ideal, f_axis, fs, N, b_butt, a_butt, true, ns_i, true, ...
        scale_factor, delay_samples);
end

figure('Position', [100 100 900 550]);
semilogy(SNR_sweep, NMSE_awgn_only, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b'); hold on;
semilogy(SNR_sweep, NMSE_phase_awgn, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');
semilogy(SNR_sweep, NMSE_lpf_awgn, 'g-d', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');
semilogy(SNR_sweep, NMSE_all_three, 'k-^', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'k');

grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('NMSE (log scale)', 'FontSize', 12);
title('DSB-SC Demodulation Performance: Single vs Combined Impairments', 'FontSize', 13);
legend('AWGN only', ...
    sprintf('Phase error (%d°) + AWGN', delta_phi_deg), ...
    sprintf('Non-ideal LPF (order %d) + AWGN', butter_order), ...
    'All three combined', ...
    'Location', 'northeast');
set(gca, 'FontSize', 11);

fprintf('--- Figure 5: SNR sweep with combined impairments generated ---\n');

fprintf('\n============================================================\n');
fprintf('                   COMBINED IMPAIRMENTS SUMMARY\n');
fprintf('============================================================\n\n');

fprintf('%-35s %-15s %-12s\n', 'Scenario', 'NMSE', 'NMSE (dB)');
fprintf('%-35s %-15s %-12s\n', repmat('-',1,35), repmat('-',1,15), repmat('-',1,12));
for s = 1:n_scenarios
    fprintf('%-35s %-15.4e %-12.2f\n', ...
        scenario_names{s}, NMSE_scenarios(s), 10*log10(NMSE_scenarios(s) + eps));
end

fprintf('\n--- ADDITIVITY TEST ---\n');
fprintf('%-20s %-15s %-15s %-10s\n', 'Combination', 'Predicted', 'Actual', 'Ratio');
fprintf('%-20s %-15s %-15s %-10s\n', repmat('-',1,20), repmat('-',1,15), repmat('-',1,15), repmat('-',1,10));
for k = 1:4
    ratio = actual_vals(k) / predicted_vals(k);
    fprintf('%-20s %-15.4e %-15.4e %-10.2f\n', ...
        combo_names{k}, predicted_vals(k), actual_vals(k), ratio);
end

fprintf('\n--- DOMINANT IMPAIRMENT ---\n');
[~, dom_idx] = max(individual_nmse);
fprintf('At the tested levels (phase=%d°, Butter-%d, SNR=%d dB):\n', ...
    delta_phi_deg, butter_order, SNR_dB_combined);
fprintf('Dominant impairment: %s (%.1f%% of sum)\n', ...
    individual_labels{dom_idx}, contributions_pct(dom_idx));

fprintf('\n--- ENGINEERING INSIGHT ---\n');
fprintf('When all three impairments are present:\n');
fprintf('  - Fix the dominant impairment first for maximum improvement.\n');
fprintf('  - At low SNR, AWGN dominates — improving sync or filter has little effect.\n');
fprintf('  - At high SNR, deterministic impairments (phase, filter) become the bottleneck.\n');
fprintf('  - This crossover illustrates why real system design involves trade-offs.\n');
fprintf('============================================================\n');

function [m_hat, nmse_val] = run_demod(s_in, m_ref, t, fc, dphi, ...
        use_phase, H_lp, ~, fs, Ns, b_filt, a_filt, use_butter, ...
        n_std, use_awgn, ext_scale, gd_samples)
    %#ok<*INUSD>  % All arguments are used inside conditional branches
   
    if use_awgn
        noise = n_std * randn(1, Ns);
        r_t = s_in + noise;
    else
        r_t = s_in;
    end
    if use_phase
        lo_t = cos(2*pi*fc*t + dphi);
    else
        lo_t = cos(2*pi*fc*t);
    end
    v_t = r_t .* lo_t;

    if use_butter
        m_raw = filter(b_filt, a_filt, v_t);
        m_hat = [m_raw(gd_samples+1:end), zeros(1, gd_samples)];
    else
        V_f = fftshift(fft(v_t));
        m_hat = real(ifft(ifftshift(V_f .* H_lp)));
    end

    if ext_scale > 0
        m_hat = m_hat / ext_scale;
    end
    nmse_val = sum((m_ref - m_hat).^2) / sum(m_ref.^2);
end
