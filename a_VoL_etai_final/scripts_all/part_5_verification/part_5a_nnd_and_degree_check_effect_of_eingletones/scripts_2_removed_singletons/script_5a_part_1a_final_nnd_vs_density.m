clearvars; close all; clc;

%% PARAMETERS

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

dataDir = '../../../part_3d_properties_data/part_3a_nnd_nosingle';

outDir = fullfile('figs','figs_5a_nnd');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% STORAGE

finalNND = nan(numel(speeds), numel(rhos));

%% READ DATA

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf( ...
            'nnd_nosingle_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            warning('Missing file: %s', fname);
            continue
        end

        dat = readmatrix(fpath);

        meanNND = dat(:,2);

        if numel(meanNND) < 10
            warning('Too few points: %s', fname);
            continue
        end

        finalNND(iSpeed,iRho) = mean(meanNND(end-9:end));

    end
end

%% PLOT

figure('Units','pixels','Position',[200 200 820 420]);
ax = axes('Units','pixels','Position',[90 80 680 300]);
hold(ax,'on');

cc = parula(numel(speeds));

hList  = gobjects(numel(speeds),1);
legStr = cell(numel(speeds),1);

for iSpeed = 1:numel(speeds)

    hList(iSpeed) = plot(ax, rhos, finalNND(iSpeed,:), ...
        'LineWidth',2, ...
        'Color',cc(iSpeed,:));

    legStr{iSpeed} = sprintf('$s = %.2f$', speeds(iSpeed));

end

ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$\\rho$', 'Interpreter','latex','FontSize',18);
ylabel('$NND_{\mathrm{final}}$', 'Interpreter','latex','FontSize',18);

grid(ax,'on'); box(ax,'on');
xlim([1 16])

valid = isgraphics(hList);

legend(hList(valid), legStr(valid), ...
    'Location','eastoutside', ...
    'Interpreter','latex','FontSize',16,'Box','off');

filename = fullfile(outDir,'part_1a_final_nnd_vs_density_multi_speed_nosingle.png');
print(gcf, filename, '-dpng', '-r600');

disp('Plot saved.');