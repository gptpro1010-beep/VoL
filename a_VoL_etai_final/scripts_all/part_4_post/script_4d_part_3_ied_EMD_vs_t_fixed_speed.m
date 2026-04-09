clearvars; close all; clc;

noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = 1:20:1000;

rho_styles = {'-', '--', ':'};

dataDir = '../part_3d_properties_data/part_3d_ied';
envFile = fullfile('data', 'data_4d_ied_csr','csr_envelope_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4d_part_3_ied_EMD_vs_t_fixed_speed');
if ~isfolder(outDir)
    mkdir(outDir);
end


envDat  = readmatrix(envFile);
t_s_env = envDat(:,1);
H_low   = envDat(:,2);
H_high  = envDat(:,3);
H_mid   = 0.5*(H_low + H_high);

for iSpeed = 1:numel(speeds)

    s_fixed = speeds(iSpeed);
    fprintf('=== speed = %.2f ===\n', s_fixed);

    cc = parula(numel(rhos));

    Denv_all = nan(numel(kvals), numel(rhos));
    Dmid_all = nan(numel(kvals), numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fullSubDir = fullfile(dataDir, ...
            sprintf('ied_avg_rho_%.2f_speed_%.2f_noise_%.2f', rho, s_fixed, noise));

        for ik = 1:numel(kvals)
        
            snap_idx = kvals(ik);
            fpath = fullfile(fullSubDir, sprintf('k_%05d.txt', snap_idx));

            if ~isfile(fpath)
                warning('Missing file: %s', fpath);
                continue;
            end

            dat = readmatrix(fpath);
            if isempty(dat) || size(dat,2) < 2
                warning('Bad/empty file: %s', fpath);
                continue;
            end

            t_s_obs = dat(:,1);
            H_obs   = dat(:,2);

            % --- interpolate CSR onto observation grid ---
            H_low_i  = interp1(t_s_env, H_low,  t_s_obs, 'linear', 'extrap');
            H_high_i = interp1(t_s_env, H_high, t_s_obs, 'linear', 'extrap');
            H_mid_i  = 0.5*(H_low_i + H_high_i);
            
            % --- Dmid ---
            Dmid = trapz(t_s_obs, abs(H_obs - H_mid_i));
            
            % --- Denv ---
            d = zeros(size(H_obs));
            
            above = H_obs > H_high_i;
            below = H_obs < H_low_i;
            
            d(above) = H_obs(above) - H_high_i(above);
            d(below) = H_low_i(below) - H_obs(below);
            
            Denv = trapz(t_s_obs, d);

            Denv_all(ik, iRho) = Denv;
            Dmid_all(ik, iRho) = Dmid;
        end
    end

    % Figure 1: Denv vs t
    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(rhos),1);
    legStr = cell(numel(rhos),1);

    for iRho = 1:numel(rhos)
        rho = rhos(iRho);
        col = cc(iRho,:);
        ls  = rho_styles{mod(iRho-1, numel(rho_styles)) + 1};

        y = Denv_all(:, iRho);

        hList(iRho) = plot(ax, 5*kvals(:), y, ...
            'Color', col, ...
            'LineStyle', ls, ...
            'LineWidth', 1.8);

        legStr{iRho} = sprintf('$\\rho = %.2f$', rho);
    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{cross}$','Interpreter','latex','FontSize',16);
    ylabel('$D_{env}$', 'Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    legend(hList, legStr, ...
        'Location','eastoutside', ...
        'NumColumns', 1, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    text(ax, 0.05, 0.95, sprintf('$s = %.2f$', s_fixed), ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'BackgroundColor','w', ...
        'EdgeColor','k');

    filename = fullfile(outDir, sprintf('Denv_vs_Tcross_speed_%.2f.png', s_fixed));
    print(gcf, filename, '-dpng', '-r600');
    close(gcf);

    % Figure 2: Dmid vs t
    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(rhos),1);
    legStr = cell(numel(rhos),1);

    for iRho = 1:numel(rhos)
        rho = rhos(iRho);
        col = cc(iRho,:);
        ls  = rho_styles{mod(iRho-1, numel(rho_styles)) + 1};

        y = Dmid_all(:, iRho);

        hList(iRho) = plot(ax, 5*kvals(:), y, ...
            'Color', col, ...
            'LineStyle', ls, ...
            'LineWidth', 1.8);

        legStr{iRho} = sprintf('$\\rho = %.2f$', rho);
    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{cross}$','Interpreter','latex','FontSize',16);
    ylabel('$D_{mid}$', 'Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    legend(hList, legStr, ...
        'Location','eastoutside', ...
        'NumColumns', 1, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    text(ax, 0.05, 0.95, sprintf('$s = %.2f$', s_fixed), ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'BackgroundColor','w', ...
        'EdgeColor','k');

    filename = fullfile(outDir, sprintf('Dmid_vs_Tcross_speed_%.2f.png', s_fixed));
    print(gcf, filename, '-dpng', '-r600');
    close(gcf);
end

disp('Denv and Dmid plots generated for all fixed speeds.');
