clearvars; close all; clc;

noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3e_quadrat';

outDir = fullfile('figs','figs_4e_part_1_quadrat_I');
if ~isfolder(outDir)
    mkdir(outDir);
end

[RHO, S] = meshgrid(rhos, speeds);

cmin = 0.25;
cmax = 8.00;


for ik = 1:numel(kvals)

    snap_idx = kvals(ik);

    IVAL = nan(numel(speeds), numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        for iSpeed = 1:numel(speeds)

            spd = speeds(iSpeed);

            fname = sprintf( ...
                'quadrat_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                rho, spd, noise);

            fpath = fullfile(dataDir, fname);

            if ~isfile(fpath)
                % continue;
            end


            M = readmatrix(fpath);

            if isempty(M) || size(M,2) < 2
                continue;
            end


            kk = M(:,1);

            idx = find(kk == snap_idx, 1);

            if isempty(idx)
                continue;
            end


            val = M(idx,2);   % mean_I


            if ~isfinite(val)
                continue;
            end


            if val < cmin
                val = cmin;
            end

            if val > cmax
                val = cmax;
            end


            IVAL(iSpeed, iRho) = val;

        end
    end


    figure('Units','pixels','Position',[200 200 800 800]);

    axMain = axes('Units','pixels','Position',[100 100 600 600]);
    hold(axMain,'on');

    imagesc(axMain, IVAL);
    set(axMain,'YDir','normal');
    axis(axMain,'tight');


    colormap(axMain, parula(256));

    set(axMain,'ColorScale','log');
    clim(axMain,[cmin cmax]);


    cb = colorbar(axMain);

    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;
    cb.Ticks = [0.25 0.5 1 2 4 8];

    axMain.FontSize = 16;
    axMain.XTick = 1:numel(rhos);
    axMain.XTickLabel = compose('%g', rhos);
    axMain.YTick = 1:numel(speeds);
    axMain.YTickLabel = compose('%.2f', speeds);
    axMain.TickLabelInterpreter = 'latex';


    xlabel(axMain,'$\rho$','Interpreter','latex','FontSize',16);
    ylabel(axMain,'$s$','Interpreter','latex','FontSize',16);

    title(axMain, ...
        sprintf('$Quadrat\\ I(\\rho,s)$ at $T/T_{cross}=%d$',5*snap_idx), ...
        'Interpreter','latex','FontSize',16);


    grid(axMain,'on');
    box(axMain,'on');


    fname = fullfile(outDir, ...
        sprintf('phase_quadrat_I_rho_speed_snap_%05d.png',snap_idx));

    print(gcf, fname, '-dpng', '-r600');

end


disp('Quadrat I phase diagrams generated.');