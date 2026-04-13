clearvars;
close all;
clc;

%% PARAMETERS

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

dataDir = '../../part_3d_properties_data/part_3b_degree_nosingle';

outDir = fullfile('figs','figs_5b_degree');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% MATRIX TO STORE FINAL VALUES

finalK = nan(numel(rhos), numel(speeds));

%% LOOP THROUGH PARAMETER SPACE

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

        %% TAKE MEAN OF LAST 'few' VALUES

        n = length(meanK);

        if n < 10
            warning('File too short: %s', fname);
            continue
        end

        finalK(iRho,iSpeed) = mean(meanK(end-9:end));

    end
end

%% HEATMAP

figure('Units','pixels','Position',[200 200 700 500]);

imagesc(speeds, rhos, finalK);

set(gca,'YDir','normal')

colormap(parula)
colorbar

ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$s$', 'Interpreter','latex','FontSize',18)
ylabel('$\rho$', 'Interpreter','latex','FontSize',18)

title('Final $\langle k \rangle$', ...
    'Interpreter','latex','FontSize',18)

grid off
box on

%% SAVE FIGURE
filename = fullfile(outDir,'part_1c_final_degree_heatmap_nosingle.png');
print(gcf, filename, '-dpng', '-r600');
disp('Heatmap saved.');