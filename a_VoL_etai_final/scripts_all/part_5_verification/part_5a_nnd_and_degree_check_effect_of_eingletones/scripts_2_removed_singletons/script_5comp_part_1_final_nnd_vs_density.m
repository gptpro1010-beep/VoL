clearvars; close all; clc;

%% PARAMETERS

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

dataDir_with = '../../../part_3d_properties_data/part_3a_nnd';
dataDir_no   = '../../../part_3d_properties_data/part_3a_nnd_nosingle';

outDir = fullfile('figs','figs_5a_nnd');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% STORAGE

finalNND_with = nan(numel(speeds), numel(rhos));
finalNND_no   = nan(numel(speeds), numel(rhos));

%% READ DATA

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        % ===== WITH SINGLETONS =====
        fname1 = sprintf( ...
            'nnd_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);
        
        fpath1 = fullfile(dataDir_with, fname1);
        
        if isfile(fpath1)
            dat1 = readmatrix(fpath1);
            m1 = dat1(:,2);
        
            if numel(m1) >= 10
                finalNND_with(iSpeed,iRho) = mean(m1(end-9:end));
            end
        end
        
        % ===== WITHOUT SINGLETONS =====
        fname2 = sprintf( ...
            'nnd_nosingle_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);
        
        fpath2 = fullfile(dataDir_no, fname2);
        
        if isfile(fpath2)
            dat2 = readmatrix(fpath2);
            m2 = dat2(:,2);
        
            if numel(m2) >= 10
                finalNND_no(iSpeed,iRho) = mean(m2(end-9:end));
            end
        end
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

    % ----- WITH SINGLETONS (SOLID) -----
    plot(ax, rhos, finalNND_with(iSpeed,:), ...
        'LineWidth',2, ...
        'Color',cc(iSpeed,:), ...
        'LineStyle','-');

    % ----- WITHOUT SINGLETONS (DASHED) -----
    hList(iSpeed) = plot(ax, rhos, finalNND_no(iSpeed,:), ...
        'LineWidth',2, ...
        'Color',cc(iSpeed,:), ...
        'LineStyle','--');

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

annotation('textbox', [0.40 0.75 0.3 0.1], ...
    'String', {'Solid: all particles', 'Dashed: no singletons'}, ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'EdgeColor','none');
    
filename = fullfile(outDir,'part_1a_final_nnd_vs_density_comparison.png');
print(gcf, filename, '-dpng', '-r600');

disp('Plot saved.');