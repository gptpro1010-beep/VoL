clearvars; close all; clc;
%%
noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

simMin = 1;
simMax = 20;

% Desired times
snap_idxs = [2 4 10 20 40 60 100 150 200 500 1000];

min_op = 0.9995;
max_op = 0.9999;

dataBase = '../../sim_data';

outDir = fullfile('figs','figs_4_1_op');
if ~isfolder(outDir)
    mkdir(outDir);
end

for it = 1:numel(snap_idxs)
    k_snap = snap_idxs(it);

    OP = nan(numel(speeds), numel(rhos));

    for iRho = 1:numel(rhos)
        rho = rhos(iRho);

        for iSpeed = 1:numel(speeds)
            s = speeds(iSpeed);

            caseDir = fullfile(dataBase, ...
                sprintf('vicsek_density_%.2f_speed_%.2f_noise_%.2f', rho, s, noise));

            vals = nan(simMax - simMin + 1, 1);

            k = 0;
            for sim = simMin:simMax

                fpath = fullfile(caseDir, sprintf('sim_%02d', sim), ...
                    sprintf('order_param_noise_%.2f_sim_%02d_timeseries.txt', noise, sim));

                P = readmatrix(fpath);
                
                row = k_snap + 1;
                
                k = k + 1;
                vals(k) = P(row);
            end

            OP(iSpeed, iRho) = mean(vals(1:k));
        end
    end

    figure('Units','pixels','Position',[200 200 800 800]);
    axMain = axes('Units','pixels','Position',[100 100 600 600]);
    hold(axMain,'on');

    imagesc(axMain, OP);
    set(axMain,'YDir','normal');
    axis(axMain,'tight');

    colormap(axMain, parula(256));
    clim(axMain, [min_op max_op]);

    cb = colorbar(axMain);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;
    % set(axMain,'ColorScale','log');

    axMain.FontSize = 16;
    axMain.XTick = 1:numel(rhos);
    axMain.XTickLabel = compose('%g', rhos);
    axMain.YTick = 1:numel(speeds);
    axMain.YTickLabel = compose('%.2f', speeds);
    axMain.TickLabelInterpreter = 'latex';

    xlabel(axMain, '$\rho$', 'Interpreter','latex','FontSize',16);
    ylabel(axMain, '$s$',   'Interpreter','latex','FontSize',16);

    title(axMain, sprintf('$Order\\ Parameter(\\rho,s)$ at $k=%d$', k_snap), ...
        'Interpreter','latex','FontSize',16);

    grid(axMain,'on');
    box(axMain,'on');

    fname = sprintf('phase_op_rho_speed_snap_%04d.png', k_snap);
    fpath = fullfile(outDir, fname);
    print(gcf, fpath, '-dpng', '-r600');
end