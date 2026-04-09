clearvars; close all; clc;
%%
noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01];   % fix speed (like your other plots)

dataBase = '../part_3d_properties_data/part_3h_avgheading';

outDir = fullfile('figs','figs_4h_avgheading_time');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% ===== LOOP OVER SPEED =====
for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    figure('Units','pixels','Position',[200 200 800 500]);
    ax = axes('Units','pixels','Position',[100 100 600 300]);
    hold(ax,'on');

    colors = lines(numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf( ...
            'avgheading_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataBase, fname);

        if ~isfile(fpath)
            % continue;
        end

        D = readmatrix(fpath);

        % Columns:
        % 1 = snap_idx
        % 2 = theta_mean

        k   = D(:,1);
        th  = D(:,2);
        
        t_tc = 5*k;
        th = unwrap(th);
        % th_abs = abs(th);

        %% ===== plot =====
        plot(ax, t_tc, th, 'LineWidth',2, ...
            'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));
    end

    %% ===== AXIS =====
    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel(ax, '$t/T_{cross}$', 'Interpreter','latex','FontSize',16);
    ylabel(ax, '$|\theta_{\mathrm{mean}}|$', 'Interpreter','latex','FontSize',16);

    title(ax, sprintf('$|\\theta_{mean}(t)|$ vs time at $s=%.2f$', s), ...
        'Interpreter','latex','FontSize',16);

    legend(ax,'Location','northwest','Interpreter','latex');

    grid(ax,'on');
    box(ax,'on');

    %% ===== SAVE =====
    fname = sprintf('theta_abs_vs_time_speed_%.2f.png', s);
    fpath = fullfile(outDir, fname);

    print(gcf, fpath, '-dpng', '-r600');
end