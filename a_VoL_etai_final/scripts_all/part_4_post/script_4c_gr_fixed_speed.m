clearvars; close all; clc;

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

% snapshot indices
kvals = [1 5 10 20 50 100 200 500 1000];

rho_styles = {'-', '-', '-'};

dataDir = '../part_3d_properties_data/part_3c_gr';

outDir = fullfile('figs','figs_4c_gr_fixed_speed');
if ~isfolder(outDir)
    mkdir(outDir);
end


for ik = 1:numel(kvals)

    snap_idx = kvals(ik);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        fprintf('Plotting speed = %.2f, k = %d\n', s, snap_idx);

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

            fullSubDir = fullfile(dataDir, ...
                sprintf('gr_avg_rho_%.2f_speed_%.2f_noise_%.2f', rho, s, noise));

            fpath = fullfile(fullSubDir, sprintf('k_%05d.txt', snap_idx));

            if ~isfile(fpath)
                warning('Missing file: %s', fpath);
                break
            end

            dat = readmatrix(fpath);

            if isempty(dat) || size(dat,2) < 5
                warning('Bad/empty file: %s', fpath);
                break
            end

            r_mid    = dat(:,1);
            meanDens = dat(:,4);
            stdDens  = dat(:,5);

            gVals  = meanDens / rho;
            stdSim = stdDens  / rho;

            hList(iRho) = plot(ax, r_mid, gVals, ...
                'Color', col, ...
                'LineStyle', ls, ...
                'LineWidth', 1.8);

            legStr{iRho} = sprintf('$\\rho = %.2f$', rho);

        end


        ax.FontSize = 16;
        ax.TickLabelInterpreter = 'latex';

        xlabel('$r$',   'Interpreter','latex','FontSize',16);
        ylabel('$g(r)$','Interpreter','latex','FontSize',16);

        grid(ax,'on');
        box(ax,'on');

        %ylim([0 3]);

        valid = isgraphics(hList);

        legend(hList(valid), legStr(valid), ...
            'Location','eastoutside', ...
            'NumColumns',1, ...
            'Interpreter','latex', ...
            'FontSize',16, ...
            'Box','off');

        title(sprintf('$s = %.2f,\\; T/T_{cross} = %d$', s, 5*snap_idx), ...
            'Interpreter','latex','FontSize',16);

        filename = fullfile(outDir, ...
            sprintf('gr_vs_r_speed_%.2f_snap_%05d.png', s, snap_idx));

        print(gcf, filename, '-dpng', '-r600');
        close(gcf);

    end
end

disp('All fixed-speed g(r) plots generated.');