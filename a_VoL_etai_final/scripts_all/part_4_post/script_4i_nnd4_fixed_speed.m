clearvars; close all; clc;

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

% Line styles
rho_styles = {'-', '--', ':'};

% Data location
dataDir = '../part_3d_properties_data/part_3i_nnd4';

outDir = fullfile('figs','figs_4i_nnd4');
if ~isfolder(outDir)
    mkdir(outDir);
end

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);
    fprintf('Plotting NND4 for speed = %.2f\n', s);

    cc = parula(numel(rhos));

    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(rhos),1);
    legStr = cell(numel(rhos),1);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);
        col = cc(iRho,:);
        ls  = rho_styles{mod(iRho-1, numel(rho_styles)) + 1};

        fname = sprintf( ...
            'nnd4_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            warning('Missing file: %s', fname);
            continue;
        end

        dat      = readmatrix(fpath);

        T        = dat(:,1);
        meanNND4 = dat(:,2);

        % Convert to dimensionless time
        Tnorm = 20*T;

        hList(iRho) = plot(ax, Tnorm, meanNND4, ...
            'Color', col, ...
            'LineStyle', ls, ...
            'LineWidth', 1.8);

        legStr{iRho} = sprintf('$\\rho = %.2f$', rho);

    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{\mathrm{cross}}$', 'Interpreter','latex','FontSize',16);
    ylabel('$NND4$', 'Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    valid = isgraphics(hList);

    legend(hList(valid), legStr(valid), ...
        'Location','eastoutside', ...
        'NumColumns',2, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    title(sprintf('$s = %.2f$', s), ...
        'Interpreter','latex','FontSize',16);

    filename = fullfile(outDir, ...
        sprintf('nnd4_vs_Tcross_speed_%.2f.png', s));

    print(gcf, filename, '-dpng', '-r600');
    close(gcf);

end

disp('All NND4 fixed-speed plots generated.');
