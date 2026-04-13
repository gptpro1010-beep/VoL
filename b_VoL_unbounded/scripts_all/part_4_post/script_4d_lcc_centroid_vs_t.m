clearvars; close all; clc;
%%
noise = 0.01;
N = 100;

rhos   = [1 2 4 8];
speeds = [0.01];

dataBase = '../part_3d_properties_data/part_3b_lcc_properties';

outDir = fullfile('figs','figs_4d_lcc_centroid_time');
if ~isfolder(outDir)
    mkdir(outDir);
end

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    figX = figure('Units','pixels','Position',[200 200 900 700]);
    axX = axes(figX);
    hold(axX,'on');

    figY = figure('Units','pixels','Position',[220 220 900 700]);
    axY = axes(figY);
    hold(axY,'on');

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
        cy_m   = D(:,10);

        L_init = sqrt(N / rho);
        T = max(1, round(L_init / s));
        save_every = 1 * T;
        t = save_every*k;

        plot(axX, t, cx_m, 'LineWidth',2, 'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));
        plot(axY, t, cy_m, 'LineWidth',2, 'Color', colors(iRho,:), ...
            'DisplayName', sprintf('$\\rho = %g$', rho));
    end

    axX.FontSize = 16;
    axX.TickLabelInterpreter = 'latex';
    xlabel(axX, '$t$', 'Interpreter','latex','FontSize',16);
    ylabel(axX, '$x_C$', 'Interpreter','latex','FontSize',16);
    title(axX, sprintf('$x_C(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);
    legend(axX,'Location','northwest','Interpreter','latex');
    grid(axX,'on');
    box(axX,'on');

    axY.FontSize = 16;
    axY.TickLabelInterpreter = 'latex';
    xlabel(axY, '$t$', 'Interpreter','latex','FontSize',16);
    ylabel(axY, '$y_C$', 'Interpreter','latex','FontSize',16);
    title(axY, sprintf('$y_C(t)$ at $s=%.2f$', s), 'Interpreter','latex','FontSize',16);
    legend(axY,'Location','northwest','Interpreter','latex');
    grid(axY,'on');
    box(axY,'on');

    fnameX = sprintf('lcc_centroid_x_vs_time_speed_%.2f.png', s);
    fpathX = fullfile(outDir, fnameX);
    print(figX, fpathX, '-dpng', '-r600');

    fnameY = sprintf('lcc_centroid_y_vs_time_speed_%.2f.png', s);
    fpathY = fullfile(outDir, fnameY);
    print(figY, fpathY, '-dpng', '-r600');
end
