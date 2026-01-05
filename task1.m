function [bits, s] = task1(p)
% 任务1：先随手造一段随机比特，再做FSK调制
bits = randi([0 1], p.N, 1);
% 先看看原始比特长啥样
FskHelpers.plot_bits(bits, p.Rb, '任务1：随机数据比特（100位）');

% 做FSK调制，后面任务要用到这个信号
[s, ~] = FskHelpers.fsk_modulate(bits, p.Fs, p.Rb, p.f1a, p.f0a);
% 只画前20位，颜色区分0/1，图看着更直观
FskHelpers.plot_fsk_colored(bits, p.Fs, p.Rb, s, '任务1：FSK信号（红=0，蓝=1，前20位）', 20);
end
