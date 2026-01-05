function task4(p)
% 任务4：流式检测，用100比特的“延迟”一块一块处理
bp_a = FskHelpers.make_bandpass(p.Fs, 9e3, 13e3);
bp_b = FskHelpers.make_bandpass(p.Fs, 5e3, 9e3);

fprintf('任务4运行中...\n');

if isfield(p, 'filterMode')
    mode = lower(p.filterMode);
else
    mode = 'both';
end
% 任务4只关心零相位/因果的选择，不做双通道对比
use_zero_phase = strcmp(mode, 'filtfilt') || strcmp(mode, 'both');

% 计时一下，看看跑多久
start_t = tic;
[err_a, err_b] = FskHelpers.fsk_task4_streaming(p.N_total, p.chunk_bits, ...
    p.Fs, p.Rb, p.f1a, p.f0a, p.f1b, p.f0b, bp_a, bp_b, p.SNRdB, use_zero_phase);
elapsed = toc(start_t);

% 把误码率算出来
ber_a = err_a / p.N_total;
ber_b = err_b / p.N_total;

fprintf('任务4错误数 A = %d, 总数 = %d, BER = %.6f\n', err_a, p.N_total, ber_a);
fprintf('任务4错误数 B = %d, 总数 = %d, BER = %.6f\n', err_b, p.N_total, ber_b);
fprintf('任务4耗时 = %.2f s\n', elapsed);
end
