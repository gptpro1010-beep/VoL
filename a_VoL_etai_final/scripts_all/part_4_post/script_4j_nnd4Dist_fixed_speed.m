clearvars; close all; clc;

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

% Data location from part_3j
dataDir = '../part_3d_properties_data/part_3j_nnd4_dist';

outDir = fullfile('figs','figs_4j_nnd4Dist');
if ~isfolder(outDir)
    mkdir(outDir);
end

% Fractions of simulation time (in terms of available snapshot index range)
timeFrac = [0.00 0.10 0.50 1.00];
timeCols = lines(numel(timeFrac));

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);
    fprintf('Plotting NND4 distribution for speed = %.2f\n', s);

    figure('Units','pixels','Position',[120 120 1200 750]);
    tlo = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

    hLegend = gobjects(numel(timeFrac),1);
    legText = cell(numel(timeFrac),1);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf('nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);
        fpath = fullfile(dataDir, fname);

        nexttile;
        hold on;

        if ~isfile(fpath)
            title(sprintf('\\rho = %.2f (missing)', rho), 'Interpreter','tex');
            axis off;
            continue;
        end

        dat = readmatrix(fpath, 'FileType','text', 'CommentStyle','#');

        if isempty(dat)
            title(sprintf('\\rho = %.2f (empty)', rho), 'Interpreter','tex');
            axis off;
            continue;
        end

        snap = dat(:,1);
        rval = dat(:,2);
        Fbar = dat(:,3);

        snapU = unique(snap);
        rU    = unique(rval);

        nS = numel(snapU);
        nR = numel(rU);

        % reshape to matrix [snap x r]
        Fmat = reshape(Fbar, [nR, nS]).';

        % pick target snapshots
        targetIdx = unique(max(1, min(nS, 1 + round((nS - 1) * timeFrac))));

        for it = 1:numel(targetIdx)
            idx = targetIdx(it);

            h = plot(rU, Fmat(idx,:), ...
                'LineWidth', 1.8, ...
                'Color', timeCols(it,:), ...
                'LineStyle', '-');

            if iRho == 1
                hLegend(it) = h;
                legText{it} = sprintf('k = %d (T/T_{cross} = %.2f)', ...
                    snapU(idx), 20*snapU(idx));
            end
        end

        grid on;
        box on;

        xlabel('r / L', 'Interpreter','tex');
        ylabel('F_{NND4}(r)', 'Interpreter','tex');
        title(sprintf('\\rho = %.2f', rho), 'Interpreter','tex');

        xlim([min(rU), max(rU)]);
        ylim([0, 1]);

    end

    title(tlo, sprintf('NND4 CDF at fixed speed s = %.2f', s), ...
        'Interpreter','tex', 'FontWeight','normal');

    valid = isgraphics(hLegend);
    if any(valid)
        lg = legend(hLegend(valid), legText(valid), ...
            'Location','southoutside', ...
            'NumColumns',2, ...
            'Interpreter','tex', ...
            'Box','off');
        lg.Layout.Tile = 'south';
    end

    outName = fullfile(outDir, sprintf('nnd4Dist_fixed_speed_%.2f.png', s));
    print(gcf, outName, '-dpng', '-r400');
    close(gcf);
end

disp('All NND4 distribution fixed-speed plots generated.');
