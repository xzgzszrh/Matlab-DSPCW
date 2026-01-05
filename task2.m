function task2(p, bits, s)
% Task 2: detect FSK signal from Task 1
bp = FskHelpers.make_bandpass(p.Fs, 9e3, 11e3);
FskHelpers.plot_filter_response(bp, p.Fs, '任务2：带通滤波器频率响应（10kHz）');

lp = FskHelpers.make_lowpass(p.Fs, p.Rb);
FskHelpers.plot_filter_response(lp, p.Fs, '任务2：低通滤波器频率响应');

% Energy detection after 10 kHz extraction
s_noisy = FskHelpers.add_awgn(s, p.SNRdB);
if isfield(p, 'filterMode')
    mode = lower(p.filterMode);
else
    mode = 'both';
end
use_ideal = strcmp(mode, 'filtfilt') || strcmp(mode, 'both');
use_real = strcmp(mode, 'filter') || strcmp(mode, 'both');
if ~use_ideal && ~use_real
    use_ideal = true;
    use_real = true;
end

if use_ideal
    [bits_hat_ideal, ~, env_ideal, lp_ideal] = FskHelpers.fsk_demod_energy( ...
        s_noisy, p.Fs, p.Rb, bp, lp, true, p.N);
end
if use_real
    [bits_hat_real, x_bit_real, env_real, lp_real] = FskHelpers.fsk_demod_energy( ...
        s_noisy, p.Fs, p.Rb, bp, lp, false, p.N);
    [bits_real_ref, bits_hat_real] = align_bits(bits, bits_hat_real, x_bit_real);
else
    bits_real_ref = bits;
end

Ns = round(p.Fs / p.Rb);
n_plot_bits = min(30, p.N);
idx = 1:(n_plot_bits * Ns);
t = (0:numel(s_noisy)-1)' / p.Fs;

figure;
plot(t(idx), s_noisy(idx));
grid on;
xlabel('时间(秒)'); ylabel('幅度');
title('任务2：含噪FSK信号（前3位）');

figure;
if use_ideal && use_real
    subplot(2, 1, 1);
    plot(t(idx), env_ideal(idx), 'b'); hold on;
    plot(t(idx), env_real(idx), 'r--');
    grid on;
    xlabel('时间(秒)'); ylabel('包络');
    legend('零相位(理想)', '因果(现实)');
    title('任务2：包络对比（前3位）');

    subplot(2, 1, 2);
    plot(t(idx), lp_ideal(idx), 'b'); hold on;
    plot(t(idx), lp_real(idx), 'r--');
    grid on;
    xlabel('时间(秒)'); ylabel('低通输出');
    legend('零相位(理想)', '因果(现实)');
    title('任务2：低通输出对比（前3位）');
elseif use_ideal
    plot(t(idx), env_ideal(idx), 'b');
    grid on;
    xlabel('时间(秒)'); ylabel('包络');
    title('任务2：包络（零相位，前3位）');
elseif use_real
    plot(t(idx), env_real(idx), 'r');
    grid on;
    xlabel('时间(秒)'); ylabel('包络');
    title('任务2：包络（因果滤波，前3位）');
end

if use_real
    FskHelpers.plot_bits(bits_hat_real, p.Rb, '任务2：检测数据序列波形（因果滤波）');
    FskHelpers.plot_bits_compare(bits_real_ref, bits_hat_real, p.Rb, '任务2：检测比特（因果滤波，延时已对齐）');
    bits_str = char(bits_hat_real' + '0');
    fprintf('任务2检测数据序列(filter): ...%s...\n', bits_str);
end
if use_ideal && ~use_real
    FskHelpers.plot_bits(bits_hat_ideal, p.Rb, '任务2：检测数据序列波形（零相位）');
    FskHelpers.plot_bits_compare(bits, bits_hat_ideal, p.Rb, '任务2：检测比特（零相位）');
    bits_str = char(bits_hat_ideal' + '0');
    fprintf('任务2检测数据序列(filtfilt): ...%s...\n', bits_str);
end

if use_ideal
    ber_ideal = mean(bits_hat_ideal ~= bits);
    fprintf('任务2 BER(filtfilt) = %.6f (%d / %d errors)\n', ber_ideal, sum(bits_hat_ideal ~= bits), p.N);
end
if use_real
    ber_real = mean(bits_hat_real ~= bits_real_ref);
    fprintf('任务2 BER(filter)   = %.6f (%d / %d errors)\n', ber_real, sum(bits_hat_real ~= bits_real_ref), numel(bits_real_ref));
end

% BER vs SNR curve
snr_list = 0:2:16;
if use_ideal
    ber_list_ideal = zeros(size(snr_list));
end
if use_real
    ber_list_real = zeros(size(snr_list));
end
frames = 20;

for k = 1:numel(snr_list)
    if use_ideal
        err_ideal = 0;
        total_bits_ideal = 0;
    end
    if use_real
        err_real = 0;
        total_bits_real = 0;
    end
    for n = 1:frames
        bits_k = randi([0 1], p.N, 1);
        [s_k, ~] = FskHelpers.fsk_modulate(bits_k, p.Fs, p.Rb, p.f1a, p.f0a);
        s_k = FskHelpers.add_awgn(s_k, snr_list(k));
        if use_ideal
            bits_hat_k_ideal = FskHelpers.fsk_demod_energy(s_k, p.Fs, p.Rb, bp, lp, true, p.N);
            err_ideal = err_ideal + sum(bits_hat_k_ideal ~= bits_k);
            total_bits_ideal = total_bits_ideal + p.N;
        end
        if use_real
            [bits_hat_k_real, x_bit_k_real] = FskHelpers.fsk_demod_energy(s_k, p.Fs, p.Rb, bp, lp, false, p.N);
            [bits_k_ref, bits_hat_k_real] = align_bits(bits_k, bits_hat_k_real, x_bit_k_real);
            err_real = err_real + sum(bits_hat_k_real ~= bits_k_ref);
            total_bits_real = total_bits_real + numel(bits_k_ref);
        end
    end
    if use_ideal
        ber_list_ideal(k) = err_ideal / total_bits_ideal;
    end
    if use_real
        ber_list_real(k) = err_real / total_bits_real;
    end
end

figure;
hold on;
if use_ideal
    semilogy(snr_list, ber_list_ideal, '-o', 'LineWidth', 1.2);
end
if use_real
    semilogy(snr_list, ber_list_real, '-s', 'LineWidth', 1.2);
end
grid on;
xlabel('信噪比(dB)'); ylabel('误码率');
if use_ideal && use_real
    legend('零相位(理想)', '因果(现实)');
elseif use_ideal
    legend('零相位(理想)');
else
    legend('因果(现实)');
end
title('任务2：误码率-信噪比曲线');
end

function [bits_ref, bits_hat_adj] = align_bits(bits_ref, bits_hat, metric)
    ref = double(bits_ref(:));
    sig = double(metric(:));
    ref = ref - mean(ref);
    sig = sig - mean(sig);
    maxlag = min(10, numel(ref) - 1);
    if maxlag < 1
        bits_hat_adj = bits_hat;
        return;
    end
    [c, lags] = xcorr(sig, ref, maxlag);
    [~, idx] = max(c);
    lag = lags(idx);
    if lag > 0
        bits_hat_adj = bits_hat(lag+1:end);
        bits_ref = bits_ref(1:end-lag);
    elseif lag < 0
        lag = -lag;
        bits_hat_adj = bits_hat(1:end-lag);
        bits_ref = bits_ref(lag+1:end);
    else
        bits_hat_adj = bits_hat;
    end
end
