clearvars; close all; clc;

%% Parameters
noise = 0.01;
q     = 0.10;   % NND quantile

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3f_nnd_dist';

envFile = fullfile('data', 'data_4f_nndDist_csr', ...
    'csr_envelope_nnd_dist_periodic_N_100_B_5000.txt');

envDat = readmatrix(envFile);
r_env  = envDat(:,1);
F_low  = envDat(:,2);
F_high = envDat(:,3);
F_mid  = 0.5 * (F_low + F_high);

% CSR mean quantile r_q
idx_csr = find(F_mid >= q, 1, 'first');
if isempty(idx_csr)
    error('CSR mean CDF never reaches q = %.2f', q);
end
r_q_CSR = r_env(idx_csr);

%% Output matrix: transition time
Tcross = nan(numel(speeds), numel(rhos));

%% Loop over parameters
for iRho = 1:numel(rhos)
    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)
        s = speeds(iSpeed);

        fullSubDir = fullfile(dataDir, ...
            sprintf('nnd_cdf_avg_rho_%.2f_speed_%.2f_noise_%.2f', ...
                    rho, s, noise));

        fname = sprintf( ...
            'nnd_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);
        
        fpath = fullfile(dataDir, fname);
        
        if ~isfile(fpath)
            continue;
        end
        
        dat = readmatrix(fpath);
        
        if isempty(dat) || size(dat,2) < 3
            continue;
        end
        
        snap_all = dat(:,1);
        r_all    = dat(:,2);
        F_all    = dat(:,3);
        
        found = false;
        
        for ik = 1:numel(kvals)
        
            snap_idx = kvals(ik);
        
            mask = abs(snap_all - snap_idx) < 1e-9;
        
            if ~any(mask)
                continue;
            end
        
            r_obs = r_all(mask);
            F_obs = F_all(mask);
        
            [r_obs, idx] = sort(r_obs);
            F_obs = F_obs(idx);
        
            idx_obs = find(F_obs >= q, 1, 'first');
            if isempty(idx_obs)
                continue;
            end
        
            r_q_obs = r_obs(idx_obs);
        
            if r_q_obs <= r_q_CSR
                Tcross(iSpeed, iRho) = 5 * snap_idx;
                found = true;
                break;
            end
        end
        
        if ~found
            Tcross(iSpeed, iRho) = 5 * kvals(end);
        end
    end
end


%% Plot: time-to-irregularity phase diagram
figure('Units','pixels','Position',[200 200 700 600]);
ax = axes('Units','pixels','Position',[120 100 420 420]);

imagesc(Tcross);
set(ax,'YDir','normal');

colormap(ax, parula);
cb = colorbar(ax);
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

ax.FontSize = 16;
ax.XTick = 1:numel(rhos);
ax.XTickLabel = compose('%g', rhos);
ax.YTick = 1:numel(speeds);
ax.YTickLabel = compose('%.2f', speeds);
ax.TickLabelInterpreter = 'latex';

xlabel(ax,'$\rho$','Interpreter','latex','FontSize',16);
ylabel(ax,'$s$','Interpreter','latex','FontSize',16);
title(ax, sprintf('Time to loss regularity (q = %.1f)', q), ...
      'Interpreter','latex','FontSize',16);

box(ax,'on');

%% Save
outDir = fullfile('figs','figs_4f_part_4_nndDist_time_to_get_irregular');
if ~isfolder(outDir), mkdir(outDir); end

fname = fullfile(outDir, sprintf('nnd_time_to_irregularity_q_%.2f.png', q));
print(gcf, fname, '-dpng', '-r600');

disp('Time-to-irregularity phase diagram generated.');