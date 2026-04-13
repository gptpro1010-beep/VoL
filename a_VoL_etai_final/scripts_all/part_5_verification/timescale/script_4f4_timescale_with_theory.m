clearvars; close all; clc;

%% ================= PARAMETERS =================
noise = 0.01;
q     = 0.10;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3f_nnd_dist';

envFile = fullfile('data', 'data_4f_nndDist_csr', ...
    'csr_envelope_nnd_dist_periodic_N_100_B_5000.txt');

%% ================= LOAD CSR =================
envDat = readmatrix(envFile);
r_env  = envDat(:,1);
F_low  = envDat(:,2);
F_high = envDat(:,3);
F_mid  = 0.5 * (F_low + F_high);

idx_csr = find(F_mid >= q, 1, 'first');
if isempty(idx_csr)
    error('CSR mean CDF never reaches q');
end
r_q_CSR = r_env(idx_csr);

%% ================= COMPUTE Tcross =================
Tcross = nan(numel(speeds), numel(rhos));

for iRho = 1:numel(rhos)
    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)
        s = speeds(iSpeed);

        fname = sprintf( ...
            'nnd_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            continue;
        end

        dat = readmatrix(fpath);

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
            Tcross(iSpeed, iRho) = NaN;   % ignore capped values
        end
    end
end

%% ================= THEORY =================
eta_eff = 2*pi*noise;

T_theory = nan(size(Tcross));

% ===== Anchor point for fitting tau_c =====
rho_anchor   = 2;
speed_anchor = 0.04;

[~, iRho_fit]   = min(abs(rhos - rho_anchor));
[~, iSpeed_fit] = min(abs(speeds - speed_anchor));

T_anchor = Tcross(iSpeed_fit, iRho_fit);

if isnan(T_anchor)
    error('Anchor point invalid (NaN)');
end

s_anchor = speeds(iSpeed_fit);
rho_fit  = rhos(iRho_fit);

d0_fit = 1 / sqrt(rho_fit);

% ===== Fit tau_c =====
tau_c = T_anchor * (s_anchor^2) * (eta_eff^2) / (d0_fit^2);

%% ===== Generate theory curves =====
for iRho = 1:numel(rhos)

    rho = rhos(iRho);
    d0  = 1 / sqrt(rho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        T_theory(iSpeed, iRho) = tau_c * (d0^2) / (s^2 * eta_eff^2);
    end
end

%% ================= PLOT =================
figure('Units','pixels','Position',[200 200 800 500]);
ax = axes; hold(ax,'on');

cc = parula(numel(rhos));
markers = {'o','s','d','^','v'};

for iRho = 1:numel(rhos)

    color = cc(iRho,:);

    % ===== Simulation (lines) =====
    plot(speeds, Tcross(:,iRho), '-', ...
        'Color', color, ...
        'LineWidth', 2);

    % ===== Theory (markers) =====
    plot(speeds, T_theory(:,iRho), markers{iRho}, ...
        'Color', color, ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5);
end

ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$s$','Interpreter','latex','FontSize',18);
ylabel('$T_{\mathrm{cross}}$','Interpreter','latex','FontSize',18);

title('Timescale: Simulation vs Theory', ...
    'Interpreter','latex','FontSize',16);

box on;

%% ================= LEGEND =================
leg = cell(numel(rhos),1);
for i = 1:numel(rhos)
    leg{i} = sprintf('$\\rho = %g$', rhos(i));
end
legend(leg, 'Interpreter','latex','Location','northwest');

%% ================= SAVE =================
outDir = fullfile('figs','figs_5_verification_timescale');
if ~isfolder(outDir), mkdir(outDir); end

fname = fullfile(outDir, ...
    sprintf('timescale_vs_theory_q_%.2f.png', q));

print(gcf, fname, '-dpng', '-r600');

disp('Done: timescale vs theory plot');