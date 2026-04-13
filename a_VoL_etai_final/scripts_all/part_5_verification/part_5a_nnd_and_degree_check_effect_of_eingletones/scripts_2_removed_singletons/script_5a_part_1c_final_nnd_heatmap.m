clearvars; close all; clc;

%% PARAMETERS

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

dataDir = '../../part_3d_properties_data/part_3a_nnd_nosingle';

outDir = fullfile('figs','figs_5a_nnd');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% STORAGE

finalNND = nan(numel(rhos), numel(speeds));

%% READ DATA

for iRho = 1:numel(rhos)

    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

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

        finalNND(iRho,iSpeed) = mean(meanNND(end-9:end));

    end
end

%% HEATMAP

figure('Units','pixels','Position',[200 200 700 500]);

imagesc(speeds, rhos, finalNND);
set(gca,'YDir','normal')

colormap(parula)
colorbar

ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$s$', 'Interpreter','latex','FontSize',18)
ylabel('$\\rho$', 'Interpreter','latex','FontSize',18)

grid off
box on

filename = fullfile(outDir,'part_1c_final_nnd_heatmap_nosingle.png');
print(gcf, filename, '-dpng', '-r600');

disp('Heatmap saved.');