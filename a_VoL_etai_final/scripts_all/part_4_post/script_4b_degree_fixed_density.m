clearvars; close all; clc;

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

% Line styles
speed_styles = {'-', '--', ':'};

% Data location
dataDir = '../part_3d_properties_data/part_3b_degree';

outDir = fullfile('figs','figs_4b_degree');
if ~isfolder(outDir)
    mkdir(outDir);
end


for iRho = 1:numel(rhos)

    rho = rhos(iRho);
    fprintf('Plotting rho = %.2f\n', rho);

    cc = parula(numel(speeds));

    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(speeds),1);
    legStr = cell(numel(speeds),1);

    for iSpeed = 1:numel(speeds)

        s   = speeds(iSpeed);
        col = cc(iSpeed,:);
        ls  = speed_styles{mod(iSpeed-1, numel(speed_styles)) + 1};

        fname = sprintf( ...
            'avgdeg_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            warning('Missing file: %s', fname);
            continue;
        end

        dat     = readmatrix(fpath);

        snap_idx = dat(:,1);
        meanK   = dat(:,2);
        stdSim  = dat(:,3);

        % Convert to dimensionless time
        Tnorm    = 5 * snap_idx;

        % Plot
        hList(iSpeed) = plot(ax, Tnorm, meanK, ...
            'Color', col, ...
            'LineStyle', ls, ...
            'LineWidth', 1.8);

        legStr{iSpeed} = sprintf('$s = %.2f$', s);

    end


    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$T/T_{\mathrm{cross}}$', 'Interpreter','latex','FontSize',16);
    ylabel('$\langle k \rangle$', 'Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    % xlim([0 4000]);
    % xticks(0:500:4000);

    % Remove empty legend entries if files missing
    valid = isgraphics(hList);

    legend(hList(valid), legStr(valid), ...
        'Location','eastoutside', ...
        'NumColumns',2, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    title(sprintf('$\\rho = %.2f$', rho), ...
        'Interpreter','latex','FontSize',16);

    filename = fullfile(outDir, ...
        sprintf('avgdeg_vs_Tcross_rho_%.2f.png', rho));

    print(gcf, filename, '-dpng', '-r600');
    close(gcf);

end


disp('All density plots generated.');