clearvars; close all; clc;

noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = 1:20:1000;

speed_styles = {'-', '--', ':'};

dataDir = '../part_3d_properties_data/part_3d_ied';
envFile = fullfile('data', 'data_4d_ied_csr','csr_envelope_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4d_part_3_ied_EMD_vs_t_fixed_density');
if ~isfolder(outDir)
    mkdir(outDir);
end


envDat  = readmatrix(envFile);
t_s_env = envDat(:,1);
H_low   = envDat(:,2);
H_high  = envDat(:,3);
H_mid   = 0.5*(H_low + H_high);


for iRho = 1:numel(rhos)

    rho_fixed = rhos(iRho);

    fprintf('=== rho = %.2f ===\n', rho_fixed);

    cc = parula(numel(speeds));

    Denv_all = nan(numel(kvals), numel(speeds));
    Dmid_all = nan(numel(kvals), numel(speeds));


    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        fullSubDir = fullfile(dataDir, ...
            sprintf('ied_avg_rho_%.2f_speed_%.2f_noise_%.2f', rho_fixed, s, noise));

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

            % --- interpolate CSR onto simulation grid ---
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

            Denv_all(ik, iSpeed) = Denv;
            Dmid_all(ik, iSpeed) = Dmid;

        end
    end


    %% --------- Figure 1 : Denv vs T/Tcross ---------

    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(speeds),1);
    legStr = cell(numel(speeds),1);

    for iSpeed = 1:numel(speeds)

        s   = speeds(iSpeed);
        col = cc(iSpeed,:);
        ls  = speed_styles{mod(iSpeed-1,numel(speed_styles))+1};

        y = Denv_all(:,iSpeed);

        hList(iSpeed) = plot(ax, 5*kvals(:), y, ...
            'Color',col, ...
            'LineStyle',ls, ...
            'LineWidth',1.8);

        legStr{iSpeed} = sprintf('$s = %.2f$', s);

    end


    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{cross}$','Interpreter','latex','FontSize',16);
    ylabel('$D_{env}$','Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    legend(hList,legStr, ...
        'Location','eastoutside', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');


    text(ax,0.05,0.95,sprintf('$\\rho = %.2f$',rho_fixed), ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'BackgroundColor','w', ...
        'EdgeColor','k');


    filename = fullfile(outDir, sprintf('Denv_vs_Tcross_rho_%.2f.png', rho_fixed));
    print(gcf,filename,'-dpng','-r600');
    close(gcf);


    %% --------- Figure 2 : Dmid vs T/Tcross ---------

    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    for iSpeed = 1:numel(speeds)

        s   = speeds(iSpeed);
        col = cc(iSpeed,:);
        ls  = speed_styles{mod(iSpeed-1,numel(speed_styles))+1};

        y = Dmid_all(:,iSpeed);

        hList(iSpeed) = plot(ax, 5*kvals(:), y, ...
            'Color',col, ...
            'LineStyle',ls, ...
            'LineWidth',1.8);

    end


    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{cross}$','Interpreter','latex','FontSize',16);
    ylabel('$D_{mid}$','Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    legend(hList,legStr, ...
        'Location','eastoutside', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    text(ax,0.05,0.95,sprintf('$\\rho = %.2f$',rho_fixed), ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'BackgroundColor','w', ...
        'EdgeColor','k');


    filename = fullfile(outDir, sprintf('Dmid_vs_Tcross_rho_%.2f.png', rho_fixed));
    print(gcf,filename,'-dpng','-r600');
    close(gcf);

end


disp('CSR metric plots generated.');