clearvars; close all; clc;

%% Parameters
noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3d_ied';
envFile = fullfile('data', 'data_4d_ied_csr','csr_envelope_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4d_part_4_ied_phase_like');
if ~isfolder(outDir), mkdir(outDir); end

%% Load CSR envelope
envDat  = readmatrix(envFile);
t_s_env = envDat(:,1);
H_low   = envDat(:,2);
H_high  = envDat(:,3);
H_mid   = 0.5*(H_low + H_high);

%% Grid
[RHO, S] = meshgrid(rhos, speeds);

%% Loop over tvals
for ik = 1:numel(kvals)
    snap_idx = kvals(ik);

    DENV = nan(numel(speeds), numel(rhos));
    DMID = nan(numel(speeds), numel(rhos));

    for iRho = 1:numel(rhos)
        rho = rhos(iRho);

        for iSpeed = 1:numel(speeds)
            s = speeds(iSpeed);

            fullSubDir = fullfile(dataDir, ...
                sprintf('ied_avg_rho_%.2f_speed_%.2f_noise_%.2f', rho, s, noise));

            fpath = fullfile(fullSubDir, sprintf('k_%05d.txt', snap_idx));

            dat = readmatrix(fpath);

            t_s_obs = dat(:,1);
            H_obs   = dat(:,2);

            % --- interpolate CSR onto observation grid ---
            H_low_i  = interp1(t_s_env, H_low,  t_s_obs, 'linear', 'extrap');
            H_high_i = interp1(t_s_env, H_high, t_s_obs, 'linear', 'extrap');
            H_mid_i  = 0.5*(H_low_i + H_high_i);
            
            % --- Dmid ---
            DMID(iSpeed, iRho) = trapz(t_s_obs, abs(H_obs - H_mid_i));
            
            % --- Denv ---
            d = zeros(size(H_obs));
            
            above = H_obs > H_high_i;
            below = H_obs < H_low_i;
            
            d(above) = H_obs(above) - H_high_i(above);
            d(below) = H_low_i(below) - H_obs(below);
            
            DENV(iSpeed, iRho) = trapz(t_s_obs, d);
        end
    end

    %% Figure 1: Phase diagram for Denv
    figure('Units','pixels','Position',[200 200 800 800]);
    axMain = axes('Units','pixels','Position',[100 100 600 600]);
    hold(axMain,'on');

    imagesc(axMain, DENV);
    set(axMain,'YDir','normal');
    axis(axMain,'tight');
    colormap(axMain, parula(256));

    cb = colorbar(axMain);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;

    axMain.FontSize = 16;
    ax.XTick = 1:numel(rhos);
    ax.XTickLabel = compose('%g', rhos);
    ax.YTick = 1:numel(speeds);
    ax.YTickLabel = compose('%.2f', speeds);
    axMain.TickLabelInterpreter = 'latex';
    
    xlabel(axMain, '$\rho$', 'Interpreter','latex','FontSize',16);
    ylabel(axMain, '$s$',   'Interpreter','latex','FontSize',16);
    title(axMain, sprintf('$D_{env}(\\rho,s)$ at $T/T_{cross}=%d$', 5*snap_idx), ...
        'Interpreter','latex','FontSize',16);

    grid(axMain,'on');
    box(axMain,'on');

    fname = fullfile(outDir, ...
        sprintf('phase_Denv_rho_speed_snap_%05d.png', snap_idx));
    print(gcf, fname, '-dpng', '-r600');
    close(gcf);

    %% Figure 2: Phase diagram for Dmid
    figure('Units','pixels','Position',[200 200 800 800]);
    axMain = axes('Units','pixels','Position',[100 100 600 600]);
    hold(axMain,'on');

    imagesc(axMain, DMID);
    set(axMain,'YDir','normal');
    axis(axMain,'tight');
    colormap(axMain, parula(256));

    cb = colorbar(axMain);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;

    axMain.FontSize = 16;
    ax.XTick = 1:numel(rhos);
    ax.XTickLabel = compose('%g', rhos);
    ax.YTick = 1:numel(speeds);
    ax.YTickLabel = compose('%.2f', speeds);
    axMain.TickLabelInterpreter = 'latex';
    
    xlabel(axMain, '$\rho$', 'Interpreter','latex','FontSize',16);
    ylabel(axMain, '$s$',   'Interpreter','latex','FontSize',16);
    title(axMain, sprintf('$D_{mid}(\\rho,s)$ at $T/T_{cross}=%d$', 5*snap_idx), ...
        'Interpreter','latex','FontSize',16);

    grid(axMain,'on');
    box(axMain,'on');

    fname = fullfile(outDir, ...
        sprintf('phase_Dmid_rho_speed_snap_%05d.png', snap_idx));
    print(gcf, fname, '-dpng', '-r600');
    close(gcf);

    fprintf('Done snap_idx = %d (T/Tcross = %d)\n', snap_idx, 5*snap_idx);
end

disp('CSR phase diagrams and matrices generated for all tvals.');
