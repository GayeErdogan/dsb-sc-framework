# DSB-SC Demodulation: Bridging the Theory-Practice Gap

A MATLAB simulation framework that systematically tests the hidden assumptions behind DSB-SC (Double Sideband Suppressed Carrier) demodulation. This repository accompanies the paper:

> **"From Textbook to Testbench: A MATLAB Framework for Exploring the Hidden Assumptions of DSB-SC Demodulation"**  
> Gaye Erdogan, TOBB University of Economics and Technology, 2026

## Overview

Every communications textbook derives DSB-SC demodulation under seven implicit assumptions. This framework breaks those assumptions one by one and measures the resulting degradation using a unified metric: **Normalized Mean Square Error (NMSE)**.

The key finding is the **dominance hierarchy**: at moderate operating conditions, the lowpass filter (LPF) accounts for 91.6% of total degradation, phase error contributes 8.3%, and AWGN contributes less than 0.1% — contradicting the emphasis of the standard curriculum.

## The Seven Assumptions

| # | What the textbook assumes | Tested? |
|---|--------------------------|---------|
| A1 | Carrier frequency >> message bandwidth | Qualitatively (Code 1) |
| A2 | Receiver oscillator has same frequency as transmitter | Qualitatively |
| A3 | Receiver oscillator has same phase as carrier | Quantitatively (Code 2) |
| A4 | Mixer is perfectly linear | Qualitatively |
| A5 | Lowpass filter is a perfect brick wall | Quantitatively (Code 3) |
| A6 | Channel adds no noise | Quantitatively (Code 1) |
| A7 | Channel has flat response and linear phase | Qualitatively |

## Files

| File | Description | Paper Section |
|------|-------------|---------------|
| `dsb_sc_awgn_simulation.m` | AWGN channel simulation; SNR sweep from -5 to 40 dB; single-tone vs. multi-tone comparison | Section IV-C5 |
| `dsb_sc_phase_error.m` | Phase error sweep 0-90 degrees; theory vs. simulation validation: NMSE = (1 - cos(delta_phi))^2 | Section IV-C3 |
| `dsb_sc_nonideal_lpf.m` | Butterworth filter orders 2, 4, 6; group delay compensation; frequency response analysis | Section IV-C4 |
| `dsb_sc_combined_impairments.m` | All 8 on/off combinations of phase error, LPF, AWGN; additivity test; SNR crossover sweep | Section V |
| `dsb_sc_sensitivity_analysis.m` | Parameter sweep: 5 phase angles x 3 filter orders; dominance heatmap; robustness verification | Section V-E |

## Requirements

- **MATLAB R2020a or later** (tested on R2026b)
- **Signal Processing Toolbox** (for `butter`, `filter`, `grpdelay`, `sosfilt`, `freqz`, `zp2sos`)

### Checking for the Signal Processing Toolbox

```matlab
ver('signal')
```

If the toolbox is not installed, add it via MATLAB's Add-On Explorer (Home > Add-Ons > Get Add-Ons).

## Quick Start

1. Clone the repository:
   ```
   git clone https://github.com/YOUR_USERNAME/dsb-sc-framework.git
   ```

2. Open MATLAB and navigate to the repository folder.

3. Run the scripts in order:
   ```matlab
   dsb_sc_awgn_simulation        % Step 1: AWGN effects
   dsb_sc_phase_error            % Step 2: Phase synchronization error
   dsb_sc_nonideal_lpf           % Step 3: Non-ideal lowpass filter
   dsb_sc_combined_impairments   % Step 4: All combinations
   dsb_sc_sensitivity_analysis   % Step 5: Robustness check
   ```

4. Each script generates figures and prints a summary table to the command window.

## Simulation Parameters

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Message frequency | f_m | 100 Hz |
| Carrier frequency | f_c | 5000 Hz |
| Sampling frequency | f_s | 100,000 Hz |
| Duration | T | 30 ms |
| Message amplitude | A_m | 1 |
| Carrier amplitude | A_c | 1 |

## Key Results

### Individual Impairments (at default operating point)

| Impairment | NMSE | Contribution |
|------------|------|-------------|
| Non-ideal LPF (Butterworth order 4) | 1.97 x 10^-1 | 91.6% |
| Phase error (30 degrees) | 1.79 x 10^-2 | 8.3% |
| AWGN (15 dB SNR) | 6.04 x 10^-5 | < 0.1% |

### Combined Effects

- Impairment effects are approximately **additive** (ratio R = 0.82 to 1.13)
- LPF remains dominant in **87% of 15 tested configurations**
- The **crossover effect** shows that increasing transmit power cannot reduce errors caused by deterministic impairments (filter distortion, phase offset)

## License

This project is licensed under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. You are free to share and adapt the material, provided you give appropriate credit.

## Citation

If you use this framework in your research or teaching, please cite:

```
Erdogan, G. (2026). From Textbook to Testbench: A MATLAB Framework for Exploring 
the Hidden Assumptions of DSB-SC Demodulation. TechRxiv. 
DOI: [to be added after upload]
```

## Author

**Gaye Erdogan**  
Department of Electrical and Electronics Engineering  
TOBB University of Economics and Technology, Ankara, Turkey

## Acknowledgment

The author used Claude (Anthropic) for assistance with manuscript drafting, code review, and language editing. All research design, simulations, analysis, and interpretation of results were performed by the author.
