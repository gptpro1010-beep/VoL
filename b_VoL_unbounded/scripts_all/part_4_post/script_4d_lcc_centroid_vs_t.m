clearvars; close all; clc;
%%
noise = 0.01;
save_every = 50;

rhos   = [1 2 4 8];
speeds = [0.01];

dataBase = '../part_3d_properties_data/part_3b_lcc_properties';

outDir = fullfile('figs','figs_4d_lcc_centroid_time');
if ~isfolder(outDir)
    mkdir(outDir);
end

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    fig = figure('Units','pixels','Position',[200 200 900 700]);

    ax1 = subplot(2,1,1);
    hold(ax1,'on');
    ax2 = subplot(2,1,2);
    hold(ax2,'on');

    colors = lines(numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf('lcc_properties_rho_%.2f_speed_%.2f_noise_%.2f.txt', rho, s, noise);
        fpath = fullfile(dataBase, fname);

        if ~isfile(fpath)
            continue;
        end

        D = readmatrix(fpath);

        k      = D(:,1);
        cx_m   = D(:,8);
        cx_s   = D(:,9);
        cy_m   = D(:,10);
        cy_s   = D(:,11);

        t = save_every*k;

        plot(ax1, t, cx_m, 'LineWidth',2, 'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));
        plot(ax2, t, cy_m, 'LineWidth',2, 'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));

        % errorbar(ax1, t, cx_m, cx_s, 'LineStyle','none', 'Color', colors(iRho,:), 'CapSize', 0);
        % errorbar(ax2, t, cy_m, cy_s, 'LineStyle','none', 'Color', colors(iRho,:), 'CapSize', 0);
    end

    ax1.FontSize = 16;
    ax1.TickLabelInterpreter = 'latex';
    ylabel(ax1, '$x_C$', 'Interpreter','latex','FontSize',16);
    title(ax1, sprintf('$x_C(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);
    legend(ax1,'Location','northwest','Interpreter','latex');
    grid(ax1,'on');
    box(ax1,'on');

    ax2.FontSize = 16;
    ax2.TickLabelInterpreter = 'latex';
    xlabel(ax2, '$t$', 'Interpreter','latex','FontSize',16);
    ylabel(ax2, '$y_C$', 'Interpreter','latex','FontSize',16);
    title(ax2, sprintf('$y_C(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);
    legend(ax2,'Location','northwest','Interpreter','latex');
    grid(ax2,'on');
    box(ax2,'on');

    fname = sprintf('lcc_centroid_vs_time_speed_%.2f.png', s);
    fpath = fullfile(outDir, fname);
    print(fig, fpath, '-dpng', '-r600');
end
