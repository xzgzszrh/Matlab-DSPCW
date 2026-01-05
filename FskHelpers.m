classdef FskHelpers
    % FSK相关的工具箱，任务里需要的零碎函数都放这儿
    methods(Static)
        function [s, t] = fsk_modulate(bits, Fs, Rb, f1, f0)
            % 把比特拉成采样点，再用频率切换的方式做FSK
            Ns = round(Fs / Rb);
            bits_ups = repelem(bits, Ns);
            freq = f0 + (f1 - f0) * bits_ups;
            phase = 2 * pi * cumsum(freq) / Fs;
            s = cos(phase);
            t = (0:numel(s)-1)' / Fs;
        end

        function [bits_hat, x_bit] = fsk_demod_discriminator(s, Fs, Rb, f1, f0, bp, n_bits, lp, use_zero_phase)
            % 用鉴频器做检测：带通 → 限幅 → 希尔伯特 → 低通 → 判决
            Ns = round(Fs / Rb);
            if nargin < 7
                % 没给bit数就按信号长度来估一个
                n_bits = floor(numel(s) / Ns);
            end
            if nargin < 8 || isempty(lp)
                % 没有低通就现场造一个
                lp = FskHelpers.make_lowpass(Fs, Rb);
            end
            if nargin < 9 || isempty(use_zero_phase)
                % 默认走零相位
                use_zero_phase = true;
            end

            if use_zero_phase
                % filtfilt是零相位，看起来更理想
                filt_fn = @filtfilt;
            else
                % filter是因果，现实里就是这个味
                filt_fn = @filter;
            end

            x1 = filt_fn(bp, s);
            % 限幅一下，鉴频器更稳
            x2 = sign(x1);
            phi = unwrap(angle(hilbert(x2)));
            f_inst = [0; diff(phi)] * Fs / (2 * pi);

            x3 = filt_fn(lp, f_inst);

            % 按bit长度切片，再取平均当作判决统计量
            x3 = x3(1:n_bits * Ns);
            X = reshape(x3, Ns, n_bits);
            x_bit = mean(X, 1)';

            % 阈值就用两频率的中点
            threshold = (f0 + f1) / 2;
            bits_hat = x_bit < threshold;
        end

        function bp = make_bandpass(Fs, f_low, f_high)
            % 带通滤波器，IIR就够用
            bp = designfilt('bandpassiir', ...
                'FilterOrder', 6, ...
                'HalfPowerFrequency1', f_low, ...
                'HalfPowerFrequency2', f_high, ...
                'SampleRate', Fs);
        end

        function lp = make_lowpass(Fs, Rb)
            % 低通滤波器，截止频率按比特率来
            lp = designfilt('lowpassiir', ...
                'FilterOrder', 6, ...
                'HalfPowerFrequency', Rb, ...
                'SampleRate', Fs);
        end

        function plot_bits(bits, Rb, plot_title)
            % 画阶梯图，看比特序列最直观
            N = numel(bits);
            t_bit = (0:N-1)' / Rb;
            figure;
            stairs(t_bit, bits, 'LineWidth', 1.5);
            ylim([-0.2 1.2]); grid on;
            xlabel('时间(秒)'); ylabel('比特');
            title(plot_title);
        end

        function plot_bits_compare(bits, bits_hat, Rb, plot_title)
            % 原始/检测比特叠在一起，方便肉眼比对
            N = numel(bits);
            t_bit = (0:N-1)' / Rb;
            figure;
            stairs(t_bit, bits_hat, 'LineWidth', 1.5); hold on;
            stairs(t_bit, bits, '--', 'LineWidth', 1.2);
            ylim([-0.2 1.2]); grid on;
            xlabel('时间(秒)'); ylabel('比特');
            legend('检测', '原始');
            title(plot_title);
        end

        function plot_fsk(t, s, plot_title)
            % 普通波形图
            figure;
            plot(t, s);
            grid on;
            xlabel('时间(秒)'); ylabel('幅度');
            title(plot_title);
        end

        function plot_fsk_colored(bits, Fs, Rb, s, plot_title, n_plot_bits)
            % 按比特上色画FSK，1是蓝，0是红
            if nargin < 6 || isempty(n_plot_bits)
                n_plot_bits = min(20, numel(bits));
            end
            Ns = round(Fs / Rb);
            n_plot_bits = min(n_plot_bits, numel(bits));
            n_samples = n_plot_bits * Ns;
            t = (0:n_samples-1)' / Fs;

            figure;
            hold on;
            for k = 1:n_plot_bits
                idx = (k-1) * Ns + (1:Ns);
                % 颜色区分比特，图一眼就能看懂
                if bits(k) == 1
                    plot(t(idx), s(idx), 'b');
                else
                    plot(t(idx), s(idx), 'r');
                end
            end
            grid on;
            xlabel('时间(秒)'); ylabel('幅度');
            title(plot_title);
            legend('比特=1', '比特=0');
        end

        function plot_filter_response(filt_obj, Fs, plot_title)
            % 画滤波器频响，顺手把坐标名改成中文
            figure;
            freqz(filt_obj, 4096, Fs);
            title(plot_title);
            ax = findall(gcf, 'Type', 'axes');
            for k = 1:numel(ax)
                xlabel(ax(k), '频率(Hz)');
                ylab = get(get(ax(k), 'YLabel'), 'String');
                % MATLAB默认标签是英文，这里统一成中文
                if ischar(ylab) && contains(ylab, 'Phase')
                    ylabel(ax(k), '相位(弧度)');
                else
                    ylabel(ax(k), '幅度(dB)');
                end
            end
        end

        function s_noisy = add_awgn(s, snr_db)
            % 按SNR加高斯白噪声，snr<=0就当不加
            if snr_db <= 0
                s_noisy = s;
                return;
            end
            sig_power = mean(abs(s).^2);
            noise_power = sig_power / (10^(snr_db/10));
            noise = sqrt(noise_power) * randn(size(s));
            s_noisy = s + noise;
        end

        function export_figures_pdf(pdf_path)
            % 把所有figure打包成一个PDF，省得一个个导
            figs = findall(0, 'Type', 'figure');
            if isempty(figs)
                fprintf('没有图可以导出。\n');
                return;
            end
            [~, order] = sort([figs.Number]);
            figs = figs(order);

            out_dir = fileparts(pdf_path);
            if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                % 目录不存在就先创建
                mkdir(out_dir);
            end

            if exist(pdf_path, 'file')
                % 老文件先删掉，避免追加混乱
                delete(pdf_path);
            end

            for k = 1:numel(figs)
                figure(figs(k));
                drawnow;
                try
                    % 优先用exportgraphics，矢量质量更好
                    if k == 1
                        exportgraphics(figs(k), pdf_path, 'ContentType', 'vector');
                    else
                        exportgraphics(figs(k), pdf_path, 'ContentType', 'vector', 'Append', true);
                    end
                catch
                    % 兼容一下老版本MATLAB
                    if k == 1
                        print(figs(k), pdf_path, '-dpdf', '-bestfit');
                    else
                        print(figs(k), pdf_path, '-dpdf', '-bestfit', '-append');
                    end
                end
            end

            fprintf('已导出%d张图到%s\n', numel(figs), pdf_path);
        end

        function [bits_hat, x_bit, env, lp_out] = fsk_demod_energy(s, Fs, Rb, bp, lp, use_zero_phase, n_bits)
            % 能量检测：带通→包络→低通→按bit取平均判决
            Ns = round(Fs / Rb);
            if nargin < 7 || isempty(n_bits)
                % 没给bit数就按长度估算
                n_bits = floor(numel(s) / Ns);
            end
            if nargin < 6 || isempty(use_zero_phase)
                % 默认零相位
                use_zero_phase = true;
            end

            if use_zero_phase
                % 零相位滤波
                filt_fn = @filtfilt;
            else
                % 因果滤波
                filt_fn = @filter;
            end

            x1 = filt_fn(bp, s);
            env = abs(hilbert(x1));
            lp_out = filt_fn(lp, env);

            % 切成bit块，算均值当统计量
            lp_out = lp_out(1:n_bits * Ns);
            X = reshape(lp_out, Ns, n_bits);
            x_bit = mean(X, 1)';
            % 阈值取最大最小的中点，简单好用
            threshold = (max(x_bit) + min(x_bit)) / 2;
            bits_hat = x_bit > threshold;
        end

        function [err_a, err_b] = fsk_task4_streaming(N_total, chunk_bits, Fs, Rb, ...
            f1a, f0a, f1b, f0b, bp_a, bp_b, snr_db, use_zero_phase)
            % 流式跑任务4：一块一块生成、调制、检测并累计误码
            err_a = 0;
            err_b = 0;
            processed = 0;
            % 每1%输出一次进度，跑长任务不容易迷路
            report_every = max(chunk_bits, floor(N_total / 100));
            start_t = tic;
            lp = FskHelpers.make_lowpass(Fs, Rb);
            if nargin < 13 || isempty(use_zero_phase)
                % 默认零相位
                use_zero_phase = true;
            end

            for start_idx = 1:chunk_bits:N_total
                % 每一块单独生成A/B比特
                n_bits = min(chunk_bits, N_total - start_idx + 1);
                bits_a = randi([0 1], n_bits, 1);
                bits_b = randi([0 1], n_bits, 1);

                % 先调制再叠加，然后加噪
                [s_a, ~] = FskHelpers.fsk_modulate(bits_a, Fs, Rb, f1a, f0a);
                [s_b, ~] = FskHelpers.fsk_modulate(bits_b, Fs, Rb, f1b, f0b);
                s_mix = s_a + s_b;
                s_mix = FskHelpers.add_awgn(s_mix, snr_db);

                % 两路分别检测
                bits_a_hat = FskHelpers.fsk_demod_discriminator(s_mix, Fs, Rb, f1a, f0a, bp_a, n_bits, lp, use_zero_phase);
                bits_b_hat = FskHelpers.fsk_demod_discriminator(s_mix, Fs, Rb, f1b, f0b, bp_b, n_bits, lp, use_zero_phase);

                % 累计误码
                err_a = err_a + sum(bits_a_hat ~= bits_a);
                err_b = err_b + sum(bits_b_hat ~= bits_b);

                processed = processed + n_bits;
                if processed == N_total || processed >= report_every
                    report_every = report_every + max(chunk_bits, floor(N_total / 100));
                    elapsed = toc(start_t);
                    pct = 100 * processed / N_total;
                    eta = elapsed * (N_total - processed) / max(processed, 1);
                    % 打印进度和预计剩余时间
                    fprintf('任务4进度: %.1f%% | 已用: %.1fs | 预计剩余: %.1fs\n', pct, elapsed, eta);
                end
            end
        end
    end
end
