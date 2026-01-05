function task3(p)
% Task 3: two FSK signals, mix, filter, detect
bits_a = randi([0 1], p.N, 1);
bits_b = randi([0 1], p.N, 1);

FskHelpers.plot_bits(bits_a, p.Rb, '任务3：信号A随机比特');
FskHelpers.plot_bits(bits_b, p.Rb, '任务3：信号B随机比特');

[s_a, t_a] = FskHelpers.fsk_modulate(bits_a, p.Fs, p.Rb, p.f1a, p.f0a);
[s_b, t_b] = FskHelpers.fsk_modulate(bits_b, p.Fs, p.Rb, p.f1b, p.f0b);

FskHelpers.plot_fsk(t_a, s_a, '任务3：FSK信号A（10/12kHz）');
FskHelpers.plot_fsk(t_b, s_b, '任务3：FSK信号B（6/8kHz）');

s_mix = s_a + s_b;
s_mix = FskHelpers.add_awgn(s_mix, p.SNRdB);
FskHelpers.plot_fsk(t_a, s_mix, '任务3：叠加后的FSK信号');

bp_a = FskHelpers.make_bandpass(p.Fs, 9e3, 13e3);
bp_b = FskHelpers.make_bandpass(p.Fs, 5e3, 9e3);
lp = FskHelpers.make_lowpass(p.Fs, p.Rb);
FskHelpers.plot_filter_response(bp_a, p.Fs, '任务3：信号A带通滤波器频率响应');
FskHelpers.plot_filter_response(bp_b, p.Fs, '任务3：信号B带通滤波器频率响应');

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
    bits_a_hat_ideal = FskHelpers.fsk_demod_discriminator(s_mix, p.Fs, p.Rb, p.f1a, p.f0a, bp_a, p.N, [], true);
    bits_b_hat_ideal = FskHelpers.fsk_demod_discriminator(s_mix, p.Fs, p.Rb, p.f1b, p.f0b, bp_b, p.N, [], true);
end
if use_real
    [bits_a_hat_real, x_a_real] = FskHelpers.fsk_demod_energy(s_mix, p.Fs, p.Rb, bp_a, lp, false, p.N);
    [bits_b_hat_real, x_b_real] = FskHelpers.fsk_demod_energy(s_mix, p.Fs, p.Rb, bp_b, lp, false, p.N);
    [bits_a_ref, bits_a_hat_real] = align_bits(bits_a, bits_a_hat_real, x_a_real);
    [bits_b_ref, bits_b_hat_real] = align_bits(bits_b, bits_b_hat_real, x_b_real);
else
    bits_a_ref = bits_a;
    bits_b_ref = bits_b;
end

if use_real
    FskHelpers.plot_bits(bits_a_hat_real, p.Rb, '任务3：信号A检测数据序列波形（因果滤波）');
    FskHelpers.plot_bits(bits_b_hat_real, p.Rb, '任务3：信号B检测数据序列波形（因果滤波）');
    FskHelpers.plot_bits_compare(bits_a_ref, bits_a_hat_real, p.Rb, '任务3：信号A检测比特（因果滤波，延时已对齐）');
    FskHelpers.plot_bits_compare(bits_b_ref, bits_b_hat_real, p.Rb, '任务3：信号B检测比特（因果滤波，延时已对齐）');
elseif use_ideal
    FskHelpers.plot_bits(bits_a_hat_ideal, p.Rb, '任务3：信号A检测数据序列波形（零相位）');
    FskHelpers.plot_bits(bits_b_hat_ideal, p.Rb, '任务3：信号B检测数据序列波形（零相位）');
    FskHelpers.plot_bits_compare(bits_a, bits_a_hat_ideal, p.Rb, '任务3：信号A检测比特（零相位）');
    FskHelpers.plot_bits_compare(bits_b, bits_b_hat_ideal, p.Rb, '任务3：信号B检测比特（零相位）');
end

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
    ber_a_ideal = mean(bits_a_hat_ideal ~= bits_a);
    ber_b_ideal = mean(bits_b_hat_ideal ~= bits_b);
    fprintf('任务3 BER (A, filtfilt) = %.6f (%d / %d errors)\n', ber_a_ideal, sum(bits_a_hat_ideal ~= bits_a), p.N);
    fprintf('任务3 BER (B, filtfilt) = %.6f (%d / %d errors)\n', ber_b_ideal, sum(bits_b_hat_ideal ~= bits_b), p.N);
end
if use_real
    ber_a_real = mean(bits_a_hat_real ~= bits_a_ref);
    ber_b_real = mean(bits_b_hat_real ~= bits_b_ref);
    fprintf('任务3 BER (A, filter)   = %.6f (%d / %d errors)\n', ber_a_real, sum(bits_a_hat_real ~= bits_a_ref), numel(bits_a_ref));
    fprintf('任务3 BER (B, filter)   = %.6f (%d / %d errors)\n', ber_b_real, sum(bits_b_hat_real ~= bits_b_ref), numel(bits_b_ref));
end
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
