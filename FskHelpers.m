classdef FskHelpers
    % Helper functions for FSK tasks
    methods(Static)
        function [s, t] = fsk_modulate(bits, Fs, Rb, f1, f0)
            Ns = round(Fs / Rb);
            bits_ups = repelem(bits, Ns);
            freq = f0 + (f1 - f0) * bits_ups;
            phase = 2 * pi * cumsum(freq) / Fs;
            s = cos(phase);
            t = (0:numel(s)-1)' / Fs;
        end

        function bits_hat = fsk_demod_discriminator(s, Fs, Rb, f1, f0, bp, n_bits, lp, use_zero_phase)
            Ns = round(Fs / Rb);
            if nargin < 7
                n_bits = floor(numel(s) / Ns);
            end
            if nargin < 8 || isempty(lp)
                lp = FskHelpers.make_lowpass(Fs, Rb);
            end
            if nargin < 9 || isempty(use_zero_phase)
                use_zero_phase = true;
            end

            if use_zero_phase
                filt_fn = @filtfilt;
            else
                filt_fn = @filter;
            end

            x1 = filt_fn(bp, s);
            x2 = sign(x1);  % limiter to stabilize discriminator
            phi = unwrap(angle(hilbert(x2)));
            f_inst = [0; diff(phi)] * Fs / (2 * pi);

            x3 = filt_fn(lp, f_inst);

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

        function plot_bits(bits, Rb, plot_title)
            N = numel(bits);
            t_bit = (0:N-1)' / Rb;
            figure;
            stairs(t_bit, bits, 'LineWidth', 1.5);
            ylim([-0.2 1.2]); grid on;
            xlabel('Time (s)'); ylabel('Bit');
            title(plot_title);
        end

        function plot_bits_compare(bits, bits_hat, Rb, plot_title)
            N = numel(bits);
            t_bit = (0:N-1)' / Rb;
            figure;
            stairs(t_bit, bits_hat, 'LineWidth', 1.5); hold on;
            stairs(t_bit, bits, '--', 'LineWidth', 1.2);
            ylim([-0.2 1.2]); grid on;
            xlabel('Time (s)'); ylabel('Bit');
            legend('Detected', 'Original');
            title(plot_title);
        end

        function plot_fsk(t, s, plot_title)
            figure;
            plot(t, s);
            grid on;
            xlabel('Time (s)'); ylabel('Amplitude');
            title(plot_title);
        end

        function plot_fsk_colored(bits, Fs, Rb, s, plot_title, n_plot_bits)
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
                if bits(k) == 1
                    plot(t(idx), s(idx), 'b');
                else
                    plot(t(idx), s(idx), 'r');
                end
            end
            grid on;
            xlabel('Time (s)'); ylabel('Amplitude');
            title(plot_title);
            legend('Bit=1', 'Bit=0');
        end

        function plot_filter_response(filt_obj, Fs, plot_title)
            figure;
            freqz(filt_obj, 4096, Fs);
            title(plot_title);
        end

        function s_noisy = add_awgn(s, snr_db)
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
            figs = findall(0, 'Type', 'figure');
            if isempty(figs)
                fprintf('No figures to export.\n');
                return;
            end
            [~, order] = sort([figs.Number]);
            figs = figs(order);

            out_dir = fileparts(pdf_path);
            if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end

            if exist(pdf_path, 'file')
                delete(pdf_path);
            end

            for k = 1:numel(figs)
                figure(figs(k));
                drawnow;
                try
                    if k == 1
                        exportgraphics(figs(k), pdf_path, 'ContentType', 'vector');
                    else
                        exportgraphics(figs(k), pdf_path, 'ContentType', 'vector', 'Append', true);
                    end
                catch
                    if k == 1
                        print(figs(k), pdf_path, '-dpdf', '-bestfit');
                    else
                        print(figs(k), pdf_path, '-dpdf', '-bestfit', '-append');
                    end
                end
            end

            fprintf('Exported %d figures to %s\n', numel(figs), pdf_path);
        end

        function [bits_hat, x_bit, env, lp_out] = fsk_demod_energy(s, Fs, Rb, bp, lp, use_zero_phase, n_bits)
            Ns = round(Fs / Rb);
            if nargin < 7 || isempty(n_bits)
                n_bits = floor(numel(s) / Ns);
            end
            if nargin < 6 || isempty(use_zero_phase)
                use_zero_phase = true;
            end

            if use_zero_phase
                filt_fn = @filtfilt;
            else
                filt_fn = @filter;
            end

            x1 = filt_fn(bp, s);
            env = abs(hilbert(x1));
            lp_out = filt_fn(lp, env);

            lp_out = lp_out(1:n_bits * Ns);
            X = reshape(lp_out, Ns, n_bits);
            x_bit = mean(X, 1)';
            threshold = (max(x_bit) + min(x_bit)) / 2;
            bits_hat = x_bit > threshold;
        end

        function [err_a, err_b] = fsk_task4_streaming(N_total, chunk_bits, Fs, Rb, ...
            f1a, f0a, f1b, f0b, bp_a, bp_b, snr_db, use_zero_phase)
            err_a = 0;
            err_b = 0;
            processed = 0;
            report_every = max(chunk_bits, floor(N_total / 100));
            start_t = tic;
            lp = FskHelpers.make_lowpass(Fs, Rb);
            if nargin < 13 || isempty(use_zero_phase)
                use_zero_phase = true;
            end

            for start_idx = 1:chunk_bits:N_total
                n_bits = min(chunk_bits, N_total - start_idx + 1);
                bits_a = randi([0 1], n_bits, 1);
                bits_b = randi([0 1], n_bits, 1);

                [s_a, ~] = FskHelpers.fsk_modulate(bits_a, Fs, Rb, f1a, f0a);
                [s_b, ~] = FskHelpers.fsk_modulate(bits_b, Fs, Rb, f1b, f0b);
                s_mix = s_a + s_b;
                s_mix = FskHelpers.add_awgn(s_mix, snr_db);

                bits_a_hat = FskHelpers.fsk_demod_discriminator(s_mix, Fs, Rb, f1a, f0a, bp_a, n_bits, lp, use_zero_phase);
                bits_b_hat = FskHelpers.fsk_demod_discriminator(s_mix, Fs, Rb, f1b, f0b, bp_b, n_bits, lp, use_zero_phase);

                err_a = err_a + sum(bits_a_hat ~= bits_a);
                err_b = err_b + sum(bits_b_hat ~= bits_b);

                processed = processed + n_bits;
                if processed == N_total || processed >= report_every
                    report_every = report_every + max(chunk_bits, floor(N_total / 100));
                    elapsed = toc(start_t);
                    pct = 100 * processed / N_total;
                    eta = elapsed * (N_total - processed) / max(processed, 1);
                    fprintf('任务4进度: %.1f%% | 已用: %.1fs | 预计剩余: %.1fs\n', pct, elapsed, eta);
                end
            end
        end
    end
end
