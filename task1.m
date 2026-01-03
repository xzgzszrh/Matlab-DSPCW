function [bits, s] = task1(p)
% Task 1: generate bits and FSK waveform
bits = randi([0 1], p.N, 1);
FskHelpers.plot_bits(bits, p.Rb, '任务1：随机数据比特（100位）');

[s, ~] = FskHelpers.fsk_modulate(bits, p.Fs, p.Rb, p.f1a, p.f0a);
FskHelpers.plot_fsk_colored(bits, p.Fs, p.Rb, s, '任务1：FSK信号（红=0，蓝=1，前20位）', 20);
end
