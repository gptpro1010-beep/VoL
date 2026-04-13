clearvars;close all;clc;

%% PARAMETERS

noise = 0.01;

rhos = [1 2 4 8 16];
speeds  = [0.01 0.02 0.04 0.08 0.16];

dataDir = '../../part_3d_properties_data/part_3b_degree_nosingle';

outDir = fullfile('figs','figs_5b_degree');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% STORAGE

finalK = nan(numel(rhos), numel(speeds));

%% READ DATA

for iRho = 1:numel(rhos)

    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        fname = sprintf( ...
            'avgdeg_nosingle_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            warning('Missing file: %s', fname);
            continue
        end

        dat = readmatrix(fpath);

        meanK = dat(:,2);

        if numel(meanK) < 10
            warning('Too few time points in %s', fname);
            continue
        end

        % Final steady-state average (last 10 points)
        finalK(iRho,iSpeed) = mean(meanK(end-9:end));

    end
end

%% PLOT

figure('Units','pixels','Position',[200 200 820 420]);
ax = axes('Units','pixels','Position',[90 80 680 300]);
hold(ax,'on');

cc = parula(numel(rhos));

hList  = gobjects(numel(rhos),1);
legStr = cell(numel(rhos),1);

for iRho = 1:numel(rhos)

    hList(iRho) = plot(ax, speeds, finalK(iRho,:), ...
        'LineWidth',2, ...
        'Color',cc(iRho,:));

    legStr{iRho} = sprintf('$\\rho = %.2f$', rhos(iRho));

end

%% AXIS FORMAT

ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$s$', 'Interpreter','latex','FontSize',18);
ylabel('$\langle k \rangle_{\mathrm{final}}$', ...
    'Interpreter','latex','FontSize',18);

grid(ax,'on');
box(ax,'on');

xlim([0.01 0.16])

valid = isgraphics(hList);

legend(hList(valid), legStr(valid), ...
    'Location','eastoutside', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'Box','off');

%% SAVE
filename = fullfile(outDir,'part_1b_final_degree_vs_speed_multi_density_nosingle.png');
print(gcf, filename, '-dpng', '-r600');

disp('Plot saved.');