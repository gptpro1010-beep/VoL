clearvars; close all; clc;

%% Parameters
noise = 0.01;
q     = 0.1;   % NND4 quantile

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3j_nnd4_dist';

envFile = fullfile('data', 'data_4j_nnd4Dist_csr', ...
    'csr_envelope_nnd4_dist_periodic_N_100_B_5000.txt');

envDat = readmatrix(envFile);
r_env  = envDat(:,1);
F4_low  = envDat(:,2);
F4_high = envDat(:,3);
F4_mid  = 0.5 * (F4_low + F4_high);

% CSR mean quantile r_q
idx_csr = find(F4_mid >= q, 1, 'first');
if isempty(idx_csr)
    error('CSR mean CDF never reaches q = %.2f', q);
end
r_q_CSR = r_env(idx_csr);

%% Output
outDir = fullfile('figs','figs_4j_part_3_nnd4_phase_like');
if ~isfolder(outDir), mkdir(outDir); end

%% Colors (parula, discrete)
cc = parula(6);
col_regular   = cc(2,:);
col_clustered = cc(4,:);
cmap = [col_clustered; col_regular];

for ik = 1:numel(kvals)
    snap_idx = kvals(ik);

    CLASS = nan(numel(speeds), numel(rhos));  % +1 regular, -1 clustered

    for iRho = 1:numel(rhos)
        rho = rhos(iRho);

        for iSpeed = 1:numel(speeds)
            s = speeds(iSpeed);

            fname = sprintf( ...
                'nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                rho, s, noise);

            fpath = fullfile(dataDir, fname);

            if ~isfile(fpath)
                continue;
            end

            dat = readmatrix(fpath);

            % columns: snap_idx  r  mean_F4  std_F4  mean_d4
            snap_all = dat(:,1);
            r_all    = dat(:,2);
            F4_all   = dat(:,3);

            mask = abs(snap_all - snap_idx) < 1e-9;

            if ~any(mask)
                continue;
            end

            r_obs = r_all(mask);
            F4_obs = F4_all(mask);

            % sort (important)
            [r_obs, idx] = sort(r_obs);
            F4_obs = F4_obs(idx);

            % observed quantile
            idx_obs = find(F4_obs >= q, 1, 'first');
            if isempty(idx_obs), continue; end
            r_q_obs = r_obs(idx_obs);

            % classification
            if r_q_obs > r_q_CSR
                CLASS(iSpeed, iRho) = +1;   % regular
            else
                CLASS(iSpeed, iRho) = -1;   % clustered
            end
        end
    end

    %% Convert to plotting matrix
    IMG = CLASS;
    IMG(IMG == -1) = 0;   % clustered
    IMG(IMG == +1) = 1;   % regular

    %% Plot
    figure('Units','pixels','Position',[200 200 600 600]);
    ax = axes('Units','pixels','Position',[100 100 400 400]);

    imagesc(IMG);
    set(ax,'YDir','normal');

    colormap(ax, cmap);
    clim([0 1]);

    ax.FontSize = 16;
    ax.XTick = 1:numel(rhos);
    ax.XTickLabel = compose('%g', rhos);
    ax.YTick = 1:numel(speeds);
    ax.YTickLabel = compose('%.2f', speeds);
    ax.TickLabelInterpreter = 'latex';
    grid(ax,'off');
    box(ax,'on');

    xlabel(ax,'$\rho$','Interpreter','latex','FontSize',16);
    ylabel(ax,'$s$','Interpreter','latex','FontSize',16);
    title(ax, sprintf('NND4 phase (q = %.1f) at $T/T_{cross} = %d$', ...
        q, 5*snap_idx), ...
        'Interpreter','latex','FontSize',16);

    cb = colorbar(ax);
    cb.Ticks = [0 1];
    cb.TickLabels = {'Clustered','Regular'};
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;

    %% Save
    fname = fullfile(outDir, ...
        sprintf('nnd4_phase_rho_speed_Tcross_%06d.png', 5*snap_idx));
    print(gcf, fname, '-dpng', '-r600');
    close(gcf);

    fprintf('Done NND4 phase t = %d\n', 5*snap_idx);
end

disp('NND4 regular/clustered phase diagrams generated.');
