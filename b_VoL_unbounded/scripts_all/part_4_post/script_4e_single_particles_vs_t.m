clearvars; close all; clc;
%%
noise = 0.01;
save_every = 200;

rhos   = [1 2 4 8];
speeds = [0.01];

dataBase = '../part_3d_properties_data/part_3c_single_particles';

outDir = fullfile('figs','figs_4e_single_particles_time');
if ~isfolder(outDir)
    mkdir(outDir);
end

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    figure('Units','pixels','Position',[200 200 800 500]);
    ax = axes('Units','pixels','Position',[100 100 600 300]);
    hold(ax,'on');

    colors = lines(numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf('single_particles_rho_%.2f_speed_%.2f_noise_%.2f.txt', rho, s, noise);
        fpath = fullfile(dataBase, fname);

        if ~isfile(fpath)
            continue;
        end

        D = readmatrix(fpath);

        k      = D(:,1);
        n_mean = D(:,2);
        n_std  = D(:,3);

        t = save_every*k;

        plot(ax, t, n_mean, 'LineWidth',2, ...
            'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));

        % errorbar(ax, t, n_mean, n_std, 'LineStyle','none', 'Color', colors(iRho,:), 'CapSize', 0);
    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel(ax, '$t$', 'Interpreter','latex','FontSize',16);
    ylabel(ax, '$N_{single}$', 'Interpreter','latex','FontSize',16);

    title(ax, sprintf('$N_{single}(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);

    legend(ax,'Location','northwest','Interpreter','latex');

    grid(ax,'on');
    box(ax,'on');

    fname = sprintf('single_particles_vs_time_speed_%.2f.png', s);
    fpath = fullfile(outDir, fname);
    print(gcf, fpath, '-dpng', '-r600');
end
