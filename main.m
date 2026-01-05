clear; clc; close all;

% Shared parameters
params.N = 100;            % number of bits for Task 1-3
params.Rb = 1000;          % bit rate (bps)
params.Fs = 200e3;         % sampling rate

params.f1a = 10e3;         % signal A: bit 1
params.f0a = 12e3;         % signal A: bit 0
params.f1b = 6e3;          % signal B: bit 1
params.f0b = 8e3;          % signal B: bit 0

params.SNRdB = 10;          % noise level (dB), 0 means no noise⚠️

params.N_total = 100000; % Task 4 total bits 正常应该是1e8 但是需要一个小时才能算完
params.chunk_bits = 100;    % Task 4 latency
params.exportPdf = false;    % true to export all figures to a single PDF
params.exportPdfPath = 'results/figures.pdf';
params.filterMode = 'both';  % 'filtfilt' | 'filter' | 'both'

% Call each task from the main script
[bits1, s1] = task1(params);
task2(params, bits1, s1);
task3(params);
task4(params);

if params.exportPdf
    FskHelpers.export_figures_pdf(params.exportPdfPath);
end
