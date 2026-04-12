clearvars; close all; clc;
%%
noise = 0.01;
save_every = 200;

rhos   = [1 2 4 8];
speeds = [0.01];

dataBase = '../part_3d_properties_data/part_3b_lcc_properties';

outDir = fullfile('figs','figs_4c_lcc_polarization_time');
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

        fname = sprintf('lcc_properties_rho_%.2f_speed_%.2f_noise_%.2f.txt', rho, s, noise);
        fpath = fullfile(dataBase, fname);

        if ~isfile(fpath)
            continue;
        end

        D = readmatrix(fpath);

        k       = D(:,1);
        p_mean  = D(:,6);
        p_std   = D(:,7);

        t = save_every*k;

        plot(ax, t, p_mean, 'LineWidth',2, ...
            'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));

        % errorbar(ax, t, p_mean, p_std, 'LineStyle','none', 'Color', colors(iRho,:), 'CapSize', 0);
    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel(ax, '$t$', 'Interpreter','latex','FontSize',16);
    ylabel(ax, '$\\phi_C$', 'Interpreter','latex','FontSize',16);

    title(ax, sprintf('$\\phi_C(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);

    legend(ax,'Location','northwest','Interpreter','latex');

    grid(ax,'on');
    box(ax,'on');

    fname = sprintf('lcc_polarization_vs_time_speed_%.2f.png', s);
    fpath = fullfile(outDir, fname);
    print(gcf, fpath, '-dpng', '-r600');
end
