clearvars; close all; clc;

%% Parameters
noise = 0.01;
q     = 0.10;   % quantile used for irregularity threshold

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDirNnd  = '../part_3d_properties_data/part_3f_nnd_dist';
dataDirNnd4 = '../part_3d_properties_data/part_3j_nnd4_dist';

envNndFile = fullfile('data', 'data_4f_nndDist_csr', ...
    'csr_envelope_nnd_dist_periodic_N_100_B_5000.txt');

envNnd4File = fullfile('data', 'data_4j_nnd4Dist_csr', ...
    'csr_envelope_nnd4_dist_periodic_N_100_B_5000.txt');

%% CSR thresholds
envNnd = readmatrix(envNndFile);
r_env_nnd = envNnd(:,1);
F_nnd_mid = 0.5 * (envNnd(:,2) + envNnd(:,3));
idx_nnd = find(F_nnd_mid >= q, 1, 'first');
if isempty(idx_nnd)
    error('NND CSR mean CDF never reaches q = %.2f', q);
end
r_q_CSR_nnd = r_env_nnd(idx_nnd);

envNnd4 = readmatrix(envNnd4File);
r_env_nnd4 = envNnd4(:,1);
F_nnd4_mid = 0.5 * (envNnd4(:,2) + envNnd4(:,3));
idx_nnd4 = find(F_nnd4_mid >= q, 1, 'first');
if isempty(idx_nnd4)
    error('NND4 CSR mean CDF never reaches q = %.2f', q);
end
r_q_CSR_nnd4 = r_env_nnd4(idx_nnd4);

%% Output matrices
Tcross_nnd  = nan(numel(speeds), numel(rhos));
Tcross_nnd4 = nan(numel(speeds), numel(rhos));

%% Loop over parameters
for iRho = 1:numel(rhos)
    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)
        s = speeds(iSpeed);

        fileNnd = fullfile(dataDirNnd, sprintf( ...
            'nnd_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', rho, s, noise));
        fileNnd4 = fullfile(dataDirNnd4, sprintf( ...
            'nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', rho, s, noise));

        if isfile(fileNnd)
            datNnd = readmatrix(fileNnd);
            snapN = datNnd(:,1);
            rN = datNnd(:,2);
            FN = datNnd(:,3);

            foundNnd = false;
            for ik = 1:numel(kvals)
                snap_idx = kvals(ik);
                m = abs(snapN - snap_idx) < 1e-9;
                if ~any(m), continue; end

                r_obs = rN(m);
                F_obs = FN(m);
                [r_obs, idx] = sort(r_obs);
                F_obs = F_obs(idx);

                iq = find(F_obs >= q, 1, 'first');
                if isempty(iq), continue; end

                if r_obs(iq) <= r_q_CSR_nnd
                    Tcross_nnd(iSpeed, iRho) = 5 * snap_idx;
                    foundNnd = true;
                    break;
                end
            end
            if ~foundNnd
                Tcross_nnd(iSpeed, iRho) = 5 * kvals(end);
            end
        end

        if isfile(fileNnd4)
            datNnd4 = readmatrix(fileNnd4);
            snapN4 = datNnd4(:,1);
            rN4 = datNnd4(:,2);
            FN4 = datNnd4(:,3);

            foundNnd4 = false;
            for ik = 1:numel(kvals)
                snap_idx = kvals(ik);
                m = abs(snapN4 - snap_idx) < 1e-9;
                if ~any(m), continue; end

                r_obs = rN4(m);
                F_obs = FN4(m);
                [r_obs, idx] = sort(r_obs);
                F_obs = F_obs(idx);

                iq = find(F_obs >= q, 1, 'first');
                if isempty(iq), continue; end

                if r_obs(iq) <= r_q_CSR_nnd4
                    Tcross_nnd4(iSpeed, iRho) = 5 * snap_idx;
                    foundNnd4 = true;
                    break;
                end
            end
            if ~foundNnd4
                Tcross_nnd4(iSpeed, iRho) = 5 * kvals(end);
            end
        end
    end
end

%% Difference: positive means NND4 takes longer than NND
DeltaT = Tcross_nnd4 - Tcross_nnd;

%% Plot
figure('Units','pixels','Position',[200 200 760 620]);
ax = axes('Units','pixels','Position',[130 100 430 430]);

imagesc(DeltaT);
set(ax,'YDir','normal');

colormap(ax, turbo);
cb = colorbar(ax);
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;
cb.Label.String = '$\Delta T = T_{\mathrm{irr}}^{\mathrm{NND4}} - T_{\mathrm{irr}}^{\mathrm{NND}}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

ax.FontSize = 16;
ax.XTick = 1:numel(rhos);
ax.XTickLabel = compose('%g', rhos);
ax.YTick = 1:numel(speeds);
ax.YTickLabel = compose('%.2f', speeds);
ax.TickLabelInterpreter = 'latex';

xlabel(ax,'$\rho$','Interpreter','latex','FontSize',16);
ylabel(ax,'$s$','Interpreter','latex','FontSize',16);
title(ax, sprintf('Difference in time to irregularity (q = %.2f)', q), ...
    'Interpreter','latex','FontSize',16);

box(ax,'on');

%% Save
outDir = fullfile('figs','figs_4j2_part_4_nnd_vs_nnd4Dist_time_to_get_irregular');
if ~isfolder(outDir), mkdir(outDir); end

fname = fullfile(outDir, ...
    sprintf('nnd_vs_nnd4_time_to_irregularity_diff_q_%.2f.png', q));
print(gcf, fname, '-dpng', '-r600');


disp('NND vs NND4 time-to-irregularity difference plot generated.');
