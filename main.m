clear; clc; close all;

% 先把环境清一清，避免上次的变量串进来

% 这里放的是几组“全局参数”，后面所有任务都会用到
params.N = 100;            % 任务1-3的比特数，先用100个练手
params.Rb = 1000;          % 比特率(bps)，直接当每秒1000比特
params.Fs = 200e3;         % 采样率，给个比较稳的200 kHz

params.f1a = 10e3;         % 信号A：比特1走10 kHz
params.f0a = 12e3;         % 信号A：比特0走12 kHz
params.f1b = 6e3;          % 信号B：比特1走6 kHz
params.f0b = 8e3;          % 信号B：比特0走8 kHz

params.SNRdB = 2;         % 噪声强度(dB)，0就是不加噪

% 任务4的总比特数，正常应该上到1e8，但那得跑一小时，先缩小规模
params.N_total = 100000;
params.chunk_bits = 100;    % 任务4分块处理的“延迟”，一块100比特
params.exportPdf = false;   % 要不要把所有图导成一个PDF
params.exportPdfPath = 'results/figures.pdf';
params.filterMode = 'both';  % 支持 'filtfilt' | 'filter' | 'both'

% 下面是按顺序跑四个任务
[bits1, s1] = task1(params);
task2(params, bits1, s1);
task3(params);
task4(params);

if params.exportPdf
    % 需要的话把图一次性打包成PDF，省得手动导出
    FskHelpers.export_figures_pdf(params.exportPdfPath);
end
