clear; clc; close all;

%% -------------------- Parameters --------------------
N_total = 100000000; % total bits
chunk_bits = 100;    % latency in bits
Rb = 1000;           % bit rate (bps)
Fs = 200e3;          % sampling rate

f1a = 10e3;          % signal A: bit 1
f0a = 12e3;          % signal A: bit 0
f1b = 6e3;           % signal B: bit 1
f0b = 8e3;           % signal B: bit 0

%% -------------------- Filters --------------------
bp_a = make_bandpass(Fs, 9e3, 13e3);
bp_b = make_bandpass(Fs, 5e3, 9e3);

fprintf('Task 4 running...\n');
start_t = tic;
[err_a, err_b] = fsk_task4_streaming(N_total, chunk_bits, Fs, Rb, ...
    f1a, f0a, f1b, f0b, bp_a, bp_b);
elapsed = toc(start_t);

ber_a = err_a / N_total;
ber_b = err_b / N_total;

fprintf('Task 4 errors A = %d, total = %d, BER = %.6f\n', err_a, N_total, ber_a);
fprintf('Task 4 errors B = %d, total = %d, BER = %.6f\n', err_b, N_total, ber_b);
fprintf('Task 4 elapsed time = %.2f s\n', elapsed);

%% ================== Local functions ==================
function [s, t] = fsk_modulate(bits, Fs, Rb, f1, f0)
    Ns = round(Fs / Rb);
    bits_ups = repelem(bits, Ns);
    freq = f0 + (f1 - f0) * bits_ups;
    phase = 2 * pi * cumsum(freq) / Fs;
    s = cos(phase);
    t = (0:numel(s)-1)' / Fs;
end

function bits_hat = fsk_demod_discriminator(s, Fs, Rb, f1, f0, bp, n_bits)
    Ns = round(Fs / Rb);
    if nargin < 7
        n_bits = floor(numel(s) / Ns);
    end

    x1 = filtfilt(bp, s);
    x2 = sign(x1);
    phi = unwrap(angle(hilbert(x2)));
    f_inst = [0; diff(phi)] * Fs / (2 * pi);

    lp = make_lowpass(Fs, Rb);
    x3 = filtfilt(lp, f_inst);

    x3 = x3(1:n_bits * Ns);
    X = reshape(x3, Ns, n_bits);
    x_bit = mean(X, 1)';

    threshold = (f0 + f1) / 2;
    bits_hat = x_bit < threshold;
end

function bp = make_bandpass(Fs, f_low, f_high)
    bp = designfilt('bandpassiir', ...
        'FilterOrder', 6, ...
        'HalfPowerFrequency1', f_low, ...
        'HalfPowerFrequency2', f_high, ...
        'SampleRate', Fs);
end

function lp = make_lowpass(Fs, Rb)
    lp = designfilt('lowpassiir', ...
        'FilterOrder', 6, ...
        'HalfPowerFrequency', Rb, ...
        'SampleRate', Fs);
end

function [err_a, err_b] = fsk_task4_streaming(N_total, chunk_bits, Fs, Rb, ...
    f1a, f0a, f1b, f0b, bp_a, bp_b)
    err_a = 0;
    err_b = 0;

    for start_idx = 1:chunk_bits:N_total
        n_bits = min(chunk_bits, N_total - start_idx + 1);
        bits_a = randi([0 1], n_bits, 1);
        bits_b = randi([0 1], n_bits, 1);

        [s_a, ~] = fsk_modulate(bits_a, Fs, Rb, f1a, f0a);
        [s_b, ~] = fsk_modulate(bits_b, Fs, Rb, f1b, f0b);
        s_mix = s_a + s_b;

        bits_a_hat = fsk_demod_discriminator(s_mix, Fs, Rb, f1a, f0a, bp_a, n_bits);
        bits_b_hat = fsk_demod_discriminator(s_mix, Fs, Rb, f1b, f0b, bp_b, n_bits);

        err_a = err_a + sum(bits_a_hat ~= bits_a);
        err_b = err_b + sum(bits_b_hat ~= bits_b);
    end
end
