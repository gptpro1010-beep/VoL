clearvars; close all; clc;
dddd
noise = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3e_quadrat';

outDir = fullfile('figs','figs_4e_part_3_quadrat_p_fisher');
if ~isfolder(outDir)
    mkdir(outDir);
end


[RHO, S] = meshgrid(rhos, speeds);

cmin = 1e-4;
cmax = 1.0;


for ik = 1:numel(kvals)

    k_sel = kvals(ik);

    PFISH = nan(numel(speeds), numel(rhos));

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        for iSpeed = 1:numel(speeds)

            spd = speeds(iSpeed);

            fname = sprintf( ...
                'quadrat_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                rho, spd, noise);

            fpath = fullfile(dataDir, fname);

            if ~isfile(fpath)
                continue;
            end


            M = readmatrix(fpath);

            if isempty(M) || size(M,2) < 4
                continue;
            end


            kk = M(:,1);

            idx = find(kk == k_sel, 1);

            if isempty(idx)
                continue;
            end


            val = M(idx,4);   % p_fisher


            if ~isfinite(val)
                continue;
            end


            if val < cmin
                val = cmin;
            end

            if val > cmax
                val = cmax;
            end


            PFISH(iSpeed, iRho) = val;

        end
    end


    figure('Units','pixels','Position',[200 200 800 800]);

    ax = axes('Units','pixels','Position',[100 100 600 600]);
    hold(ax,'on');

    imagesc(PFISH);
    set(ax,'YDir','normal');
    axis(ax, 'tight')

    colormap(ax, parula(256));

    set(ax,'ColorScale','log');
    clim(ax,[cmin cmax]);


    cb = colorbar(ax);

    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 16;

    cb.Ticks = [1e-4 1e-3 1e-2 1e-1 1];
    cb.TickLabels = {'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$1$'};


    ax.FontSize = 16;
    ax.XTick = 1:numel(rhos);
    ax.XTickLabel = compose('%g', rhos);
    ax.YTick = 1:numel(speeds);
    ax.YTickLabel = compose('%.2f', speeds);
    ax.TickLabelInterpreter = 'latex';

    xlabel(ax,'$\rho$','Interpreter','latex','FontSize',16);
    ylabel(ax,'$s$','Interpreter','latex','FontSize',16);

    title(ax, ...
        sprintf('$Quadrat\\ CSR\\ p$-value (Fisher) at $T/T_{cross}=%d$',5*k_sel), ...
        'Interpreter','latex','FontSize',16);


    grid(ax,'on');
    box(ax,'on');


    fname = fullfile(outDir, ...
        sprintf('phase_quadrat_p_fisher_rho_speed_Tcross_%05d.png',5*k_sel));

    print(gcf, fname, '-dpng', '-r600');

end


disp('Quadrat p_fisher phase diagrams generated.');