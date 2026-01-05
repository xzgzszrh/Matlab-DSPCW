function task3(p)
% 任务3：两路FSK叠加在一起，再分别把它们捞出来
bits_a = randi([0 1], p.N, 1);
bits_b = randi([0 1], p.N, 1);

% 先看两路随机比特
FskHelpers.plot_bits(bits_a, p.Rb, '任务3：信号A随机比特');
FskHelpers.plot_bits(bits_b, p.Rb, '任务3：信号B随机比特');

% 各自调制成FSK
[s_a, t_a] = FskHelpers.fsk_modulate(bits_a, p.Fs, p.Rb, p.f1a, p.f0a);
[s_b, t_b] = FskHelpers.fsk_modulate(bits_b, p.Fs, p.Rb, p.f1b, p.f0b);

% 单独看看两路的FSK长啥样
FskHelpers.plot_fsk(t_a, s_a, '任务3：FSK信号A（10/12kHz）');
FskHelpers.plot_fsk(t_b, s_b, '任务3：FSK信号B（6/8kHz）');

% 叠加并加噪
s_mix = s_a + s_b;
s_mix = FskHelpers.add_awgn(s_mix, p.SNRdB);
FskHelpers.plot_fsk(t_a, s_mix, '任务3：叠加后的FSK信号');

% 分别给A/B做带通，最后统一低通
bp_a = FskHelpers.make_bandpass(p.Fs, 9e3, 13e3);
bp_b = FskHelpers.make_bandpass(p.Fs, 5e3, 9e3);
bp_a_energy = FskHelpers.make_bandpass_fir(p.Fs, 9e3, 11e3);
bp_b_energy = FskHelpers.make_bandpass_fir(p.Fs, 5e3, 7e3);
lp_energy = FskHelpers.make_lowpass_fir(p.Fs, p.Rb);

if isfield(p, 'filterMode')
    mode = lower(p.filterMode);
else
    mode = 'both';
end
% 零相位 vs 因果滤波，方便对比效果
use_ideal = strcmp(mode, 'filtfilt') || strcmp(mode, 'both');
use_real = strcmp(mode, 'filter') || strcmp(mode, 'both');
if ~use_ideal && ~use_real
    use_ideal = true;
    use_real = true;
end

if use_ideal
    FskHelpers.plot_filter_response(bp_a, p.Fs, '任务3：信号A带通滤波器频率响应');
    FskHelpers.plot_filter_response(bp_b, p.Fs, '任务3：信号B带通滤波器频率响应');
end
if use_real
    FskHelpers.plot_filter_response(bp_a_energy, p.Fs, '任务3：信号A能量检测带通响应（10kHz）');
    FskHelpers.plot_filter_response(bp_b_energy, p.Fs, '任务3：信号B能量检测带通响应（6kHz）');
end

if use_ideal
    % 理想版本：用鉴频器直接判比特
    bits_a_hat_ideal = FskHelpers.fsk_demod_discriminator(s_mix, p.Fs, p.Rb, p.f1a, p.f0a, bp_a, p.N, [], true);
    bits_b_hat_ideal = FskHelpers.fsk_demod_discriminator(s_mix, p.Fs, p.Rb, p.f1b, p.f0b, bp_b, p.N, [], true);
end
if use_real
    % 现实版本：能量检测 + 因果滤波
    [bits_a_hat_real, x_a_real] = FskHelpers.fsk_demod_energy(s_mix, p.Fs, p.Rb, bp_a_energy, lp_energy, false, p.N);
    [bits_b_hat_real, x_b_real] = FskHelpers.fsk_demod_energy(s_mix, p.Fs, p.Rb, bp_b_energy, lp_energy, false, p.N);
    % 对齐一下延迟，不然误码率会看起来偏大
    [bits_a_ref, bits_a_hat_real] = align_bits(bits_a, bits_a_hat_real, x_a_real);
    [bits_b_ref, bits_b_hat_real] = align_bits(bits_b, bits_b_hat_real, x_b_real);
else
    bits_a_ref = bits_a;
    bits_b_ref = bits_b;
end

if use_real
    % 因果滤波的可视化
    FskHelpers.plot_bits(bits_a_hat_real, p.Rb, '任务3：信号A检测数据序列波形（因果滤波）');
    FskHelpers.plot_bits(bits_b_hat_real, p.Rb, '任务3：信号B检测数据序列波形（因果滤波）');
    FskHelpers.plot_bits_compare(bits_a_ref, bits_a_hat_real, p.Rb, '任务3：信号A检测比特（因果滤波，延时已对齐）');
    FskHelpers.plot_bits_compare(bits_b_ref, bits_b_hat_real, p.Rb, '任务3：信号B检测比特（因果滤波，延时已对齐）');
elseif use_ideal
    % 零相位滤波的可视化
    FskHelpers.plot_bits(bits_a_hat_ideal, p.Rb, '任务3：信号A检测数据序列波形（零相位）');
    FskHelpers.plot_bits(bits_b_hat_ideal, p.Rb, '任务3：信号B检测数据序列波形（零相位）');
    FskHelpers.plot_bits_compare(bits_a, bits_a_hat_ideal, p.Rb, '任务3：信号A检测比特（零相位）');
    FskHelpers.plot_bits_compare(bits_b, bits_b_hat_ideal, p.Rb, '任务3：信号B检测比特（零相位）');
end

% 打印一下检测到的比特，眼睛扫一眼就知道对不对
if use_real
    bits_a_str = char(bits_a_hat_real' + '0');
    bits_b_str = char(bits_b_hat_real' + '0');
else
    bits_a_str = char(bits_a_hat_ideal' + '0');
    bits_b_str = char(bits_b_hat_ideal' + '0');
end
fprintf('任务3检测数据序列A: ...%s...\n', bits_a_str);
fprintf('任务3检测数据序列B: ...%s...\n', bits_b_str);

if use_ideal
    % 误码率统计（理想零相位）
    ber_a_ideal = mean(bits_a_hat_ideal ~= bits_a);
    ber_b_ideal = mean(bits_b_hat_ideal ~= bits_b);
    fprintf('任务3 BER (A, filtfilt) = %.6f (%d / %d errors)\n', ber_a_ideal, sum(bits_a_hat_ideal ~= bits_a), p.N);
    fprintf('任务3 BER (B, filtfilt) = %.6f (%d / %d errors)\n', ber_b_ideal, sum(bits_b_hat_ideal ~= bits_b), p.N);
end
if use_real
    % 误码率统计（因果滤波）
    ber_a_real = mean(bits_a_hat_real ~= bits_a_ref);
    ber_b_real = mean(bits_b_hat_real ~= bits_b_ref);
    fprintf('任务3 BER (A, filter)   = %.6f (%d / %d errors)\n', ber_a_real, sum(bits_a_hat_real ~= bits_a_ref), numel(bits_a_ref));
    fprintf('任务3 BER (B, filter)   = %.6f (%d / %d errors)\n', ber_b_real, sum(bits_b_hat_real ~= bits_b_ref), numel(bits_b_ref));
end

% 误码率-信噪比 曲线（和任务2一致：零相位 vs 因果）
snr_list = -6:2:12;
frames = 30;
if use_ideal
    ber_ideal_curve = zeros(size(snr_list));
end
if use_real
    ber_real_curve = zeros(size(snr_list));
end

for k = 1:numel(snr_list)
    if use_ideal
        err_a_ideal = 0; err_b_ideal = 0; total_ideal = 0;
    end
    if use_real
        err_a_real = 0; err_b_real = 0; total_a_real = 0; total_b_real = 0;
    end
    for n = 1:frames
        bits_a_k = randi([0 1], p.N, 1);
        bits_b_k = randi([0 1], p.N, 1);
        [s_a_k, ~] = FskHelpers.fsk_modulate(bits_a_k, p.Fs, p.Rb, p.f1a, p.f0a);
        [s_b_k, ~] = FskHelpers.fsk_modulate(bits_b_k, p.Fs, p.Rb, p.f1b, p.f0b);
        s_mix_k = s_a_k + s_b_k;
        s_mix_k = FskHelpers.add_awgn(s_mix_k, snr_list(k));

        if use_ideal
            bits_a_hat_k = FskHelpers.fsk_demod_discriminator(s_mix_k, p.Fs, p.Rb, p.f1a, p.f0a, bp_a, p.N, [], true);
            bits_b_hat_k = FskHelpers.fsk_demod_discriminator(s_mix_k, p.Fs, p.Rb, p.f1b, p.f0b, bp_b, p.N, [], true);
            err_a_ideal = err_a_ideal + sum(bits_a_hat_k ~= bits_a_k);
            err_b_ideal = err_b_ideal + sum(bits_b_hat_k ~= bits_b_k);
            total_ideal = total_ideal + p.N;
        end
        if use_real
            [bits_a_hat_k, x_a_k] = FskHelpers.fsk_demod_energy(s_mix_k, p.Fs, p.Rb, bp_a_energy, lp_energy, false, p.N);
            [bits_b_hat_k, x_b_k] = FskHelpers.fsk_demod_energy(s_mix_k, p.Fs, p.Rb, bp_b_energy, lp_energy, false, p.N);
            [bits_a_ref_k, bits_a_hat_k] = align_bits(bits_a_k, bits_a_hat_k, x_a_k);
            [bits_b_ref_k, bits_b_hat_k] = align_bits(bits_b_k, bits_b_hat_k, x_b_k);
            err_a_real = err_a_real + sum(bits_a_hat_k ~= bits_a_ref_k);
            err_b_real = err_b_real + sum(bits_b_hat_k ~= bits_b_ref_k);
            total_a_real = total_a_real + numel(bits_a_ref_k);
            total_b_real = total_b_real + numel(bits_b_ref_k);
        end
    end
    if use_ideal
        ber_ideal_curve(k) = 0.5 * (err_a_ideal + err_b_ideal) / total_ideal;
    end
    if use_real
        ber_real_curve(k) = 0.5 * (err_a_real / total_a_real + err_b_real / total_b_real);
    end
end

figure;
hold on;
if use_ideal
    semilogy(snr_list, ber_ideal_curve, '-o', 'LineWidth', 1.2);
end
if use_real
    semilogy(snr_list, ber_real_curve, '-s', 'LineWidth', 1.2);
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
title('任务3：误码率-信噪比曲线');
end

function [bits_ref, bits_hat_adj] = align_bits(bits_ref, bits_hat, metric)
    % 跟task2一样，用相关性把延迟对齐掉
    ref = double(bits_ref(:));
    sig = double(metric(:));
    % 去掉均值，相关性更稳
    ref = ref - mean(ref);
    sig = sig - mean(sig);
    % 只在小范围里找延迟，省时间
    maxlag = min(10, numel(ref) - 1);
    if maxlag < 1
        bits_hat_adj = bits_hat;
        return;
    end
    % 找相关性最大的lag
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
