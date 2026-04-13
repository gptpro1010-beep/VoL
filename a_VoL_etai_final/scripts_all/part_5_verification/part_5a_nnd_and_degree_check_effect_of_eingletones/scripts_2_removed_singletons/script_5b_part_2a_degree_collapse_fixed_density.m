clearvars; close all; clc;

noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.02 0.04 0.08 0.16];

N = 100;
n_side = 10;

speed_styles = {'-', '--', ':'};

dataDir = '../../part_3d_properties_data/part_3b_degree_nosingle';

outDir  = fullfile('figs','figs_5b_degree');
if ~isfolder(outDir)
    mkdir(outDir);
end

for iRho = 1:numel(rhos)

    rho = rhos(iRho);
    fprintf('Plotting rho = %.2f\n', rho);

    % box length
    L = sqrt(N / rho);

    % lattice spacing
    d0 = L / n_side;

    cc = parula(numel(speeds));

    figure('Units','pixels','Position',[200 200 820 420]);
    ax = axes('Units','pixels','Position',[90 80 680 300]);
    hold(ax,'on');

    hList  = gobjects(numel(speeds),1);
    legStr = cell(numel(speeds),1);

    for iSpeed = 1:numel(speeds)

        s   = speeds(iSpeed);
        col = cc(iSpeed,:);
        ls  = speed_styles{mod(iSpeed-1, numel(speed_styles)) + 1};

        fname = sprintf( ...
            'avgdeg_nosingle_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            warning('Missing file: %s', fname);
            continue;
        end

        dat      = readmatrix(fpath);
        snap_idx = dat(:,1);
        meanK    = dat(:,2);

        %% ---- scaled time from theory ----
        tau = (5 * s * noise^2 * L / (d0^2)) * snap_idx;

        %% ---- plot ----
        hList(iSpeed) = plot(ax, tau, meanK, ...
            'Color', col, ...
            'LineStyle', ls, ...
            'LineWidth', 1.8);

        legStr{iSpeed} = sprintf('$s = %.2f$', s);

    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';

    xlabel('$\tau \ (\mathrm{scaled\ time})$', ...
        'Interpreter','latex','FontSize',16);
    ylabel('$\langle k \rangle$', ...
        'Interpreter','latex','FontSize',16);

    grid(ax,'on');
    box(ax,'on');

    % robust legend
    valid = isgraphics(hList);

    legend(hList(valid), legStr(valid), ...
        'Location','eastoutside', ...
        'NumColumns', 2, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Box','off');

    title(sprintf('$\\rho = %.2f$', rho), ...
        'Interpreter','latex','FontSize',16);

    filename = fullfile(outDir, ...
        sprintf('figs_5b_part_2a_degree_collapse_fixed_density_rho_%.2f_nosingle.png', rho));

    print(gcf, filename, '-dpng', '-r600');
    close(gcf);

end

disp('All scaled plots generated.');