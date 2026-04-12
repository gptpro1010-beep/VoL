clearvars;
close all;
clc;

%% Parameters
noise = 0.01;

rhos   = [1.00 4.00 16.00];
speeds = [0.01 0.04 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3j_nnd4_dist';

envFile = fullfile('data', 'data_4j_nnd4Dist_csr', ...
    'csr_envelope_nnd4_dist_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4j_part_2_nnd4Dist_csr_test');
if ~isfolder(outDir), mkdir(outDir); end

%% Load CSR envelope
envDat   = readmatrix(envFile);

t_s_env = envDat(:,1);
F4_low  = envDat(:,2);
F4_high = envDat(:,3);
F4_mid  = 0.5*(F4_low + F4_high);

%% Colors
cc = parula(6);
line_env = cc(1,:);
line_obs = cc(4,:);

%% Loop over parameters
for iRho = 1:numel(rhos)

    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        %% --- Load averaged NND4 file ---
        fpath = fullfile(dataDir, ...
            sprintf('nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                    rho, s, noise));

        if ~isfile(fpath)
            warning('Missing file: %s', fpath);
            continue;
        end

        dat = readmatrix(fpath);

        if isempty(dat) || size(dat,2) < 5
            warning('Bad file format: %s', fpath);
            continue;
        end

        % columns: snap_idx  t_s  mean_F4  std_F4  mean_d4
        snap_all = dat(:,1);
        t_s_all  = dat(:,2);
        F4_all   = dat(:,3);

        %% --- Loop over snapshots ---
        for ik = 1:numel(kvals)

            snap_idx = kvals(ik);

            mask = (snap_all == snap_idx);

            if ~any(mask)
                warning('No data for snap_idx = %d (rho=%.2f, s=%.2f)', ...
                        snap_idx, rho, s);
                continue;
            end

            %% Extract data
            t_s_obs = t_s_all(mask);
            F4_obs  = F4_all(mask);

            % --- Sort (critical) ---
            [t_s_obs, idx] = sort(t_s_obs);
            F4_obs = F4_obs(idx);

            % --- Interpolate onto CSR grid ---
            F4_obs = interp1(t_s_obs, F4_obs, t_s_env, 'linear', 'extrap');

            % Replace potential NaNs (rare but safe)
            F4_obs(isnan(F4_obs)) = 0;

            %% --- Deviation metrics ---
            Dmid = trapz(t_s_env, abs(F4_obs - F4_mid));

            d = zeros(size(F4_obs));

            above = F4_obs > F4_high;
            below = F4_obs < F4_low;

            d(above) = F4_obs(above) - F4_high(above);
            d(below) = F4_low(below) - F4_obs(below);

            Denv = trapz(t_s_env, d);

            %% --- Plot ---
            figure('Units','pixels','Position',[200 200 520 420]);

            ax = axes('Units','pixels','Position',[80 70 380 300]);
            hold(ax,'on');

            % CSR envelope
            plot(ax, t_s_env, F4_low,  '--', 'Color', line_env, 'LineWidth', 1.5);
            plot(ax, t_s_env, F4_high, '--', 'Color', line_env, 'LineWidth', 1.5);

            % Observed
            plot(ax, t_s_env, F4_obs, '-', 'Color', line_obs, 'LineWidth', 2.0);

            box(ax,'on');
            grid(ax,'on');
            axis(ax,'square');

            xlabel('$t_s$', 'Interpreter','latex','FontSize',16);
            ylabel('$H_4(t_s)$', 'Interpreter','latex','FontSize',16);

            ax.FontSize = 16;
            ax.TickLabelInterpreter = 'latex';

            xlim(ax,[t_s_env(1) t_s_env(end)]);
            ylim(ax,[0 1]);

            %% --- Annotations ---
            txt1 = sprintf('$D_{mid} = %.4f$\\n$D_{env} = %.4f$', Dmid, Denv);

            text(ax,0.95,0.75,txt1, ...
                'Units','normalized', ...
                'VerticalAlignment','top', ...
                'HorizontalAlignment','right', ...
                'Interpreter','latex', ...
                'FontSize',16, ...
                'BackgroundColor','w', ...
                'EdgeColor','k');

            txt2 = sprintf('$\\rho = %.2f$\\n$s = %.2f$\\n$T/T_{cross} = %d$', ...
                           rho, s, 5*snap_idx);

            text(ax,0.95,0.05,txt2, ...
                'Units','normalized', ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right', ...
                'Interpreter','latex', ...
                'FontSize',16, ...
                'BackgroundColor','w', ...
                'EdgeColor','k');

            %% --- Save ---
            outname = sprintf( ...
                'csr_test_nnd4_dist_rho_%.2f_speed_%.2f_snap_%05d.png', ...
                rho, s, snap_idx);

            print(gcf, fullfile(outDir, outname), '-dpng', '-r600');
            close(gcf);

        end
    end
end

disp('NND4 CSR test plots generated.');
