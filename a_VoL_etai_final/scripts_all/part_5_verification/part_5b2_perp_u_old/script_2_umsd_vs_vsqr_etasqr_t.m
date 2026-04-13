clearvars; close all; clc;

%% PARAMETERS
N = 100;
noise  = 0.01;

rhos   = [1 2 4 8 16];
speeds = [0.01 0.16];

dataDir = '../../part_3d_properties_data/part_3g_uperp';

outDir = fullfile('figs','figs_5g_uperp');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% EFFECTIVE NOISE
eta_eff = 2*pi*noise;

%% PLOT SETUP
figure('Units','pixels','Position',[200 200 820 420]);
ax = axes('Units','pixels','Position',[90 80 680 300]);
hold(ax,'on');

cc = parula(numel(speeds));
markers = {'o','s','d','^','v'};

%% LOOP
for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf( ...
            'uperp_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            continue
        end

        dat = readmatrix(fpath);

        k      = dat(:,1);
        u_msd  = dat(:,2);

        %% TIME
        L = sqrt(N / rho);
        T = round(L / s);
        dt_out = 5 * T;

        t = k * dt_out;

        %% REMOVE EARLY TIME
        valid = t > 0;
        t_plot = t(valid);
        u_plot = u_msd(valid);

        %% LATTICE SPACING
        n_side = sqrt(N);
        d0 = L / n_side;

        %% SCALED TIME (FULL THEORY)
        tau = (s^2 * eta_eff^2 .* t_plot) / (d0^2);
        u_plot = u_plot / d0^2;

        %% LOG-LOG SLOPE CHECK
        p = polyfit(log10(tau), log10(u_plot), 1);
        slope = p(1);

        fprintf('s=%.2f rho=%.2f → loglog slope = %.3f\n', s, rho, slope);

        %% PLOT
        loglog(ax, tau, u_plot, ...
            'LineStyle','none', ...
            'Color', cc(iSpeed,:), ...
            'Marker', markers{iRho}, ...
            'MarkerSize',6);

    end
end

%% AXIS
set(ax, 'XScale', 'log', 'YScale', 'log');

ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$\tau = \frac{v^2 \eta_{\mathrm{eff}}^2 t}{d_0^2}$', ...
    'Interpreter','latex','FontSize',18);

ylabel('$\langle u_\perp^2 \rangle$', ...
    'Interpreter','latex','FontSize',18);

grid(ax,'on');
box(ax,'on');

%% LEGEND (speed = color)
hSpeed = gobjects(numel(speeds),1);
legSpeed = cell(numel(speeds),1);

for iSpeed = 1:numel(speeds)
    hSpeed(iSpeed) = plot(nan, nan, ...
        'Color', cc(iSpeed,:), ...
        'LineWidth', 2);
    legSpeed{iSpeed} = sprintf('$s = %.2f$', speeds(iSpeed));
end

leg1 = legend(hSpeed, legSpeed, ...
    'Position',[0.15 0.55 0.15 0.3], ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'Box','off');

set(leg1, 'AutoUpdate', 'off');

%% LEGEND (rho = marker)
ax2 = axes('Position', ax.Position, ...
           'Color','none', ...
           'XTick',[], 'YTick',[], ...
           'Box','off');

hold(ax2,'on');

hRho = gobjects(numel(rhos),1);
legRho = cell(numel(rhos),1);

for iRho = 1:numel(rhos)
    hRho(iRho) = plot(nan, nan, ...
        'LineStyle','none', ...
        'Color','k', ...
        'Marker', markers{iRho}, ...
        'MarkerSize',6);
    legRho{iRho} = sprintf('$\\rho = %.0f$', rhos(iRho));
end

legend(ax2, hRho, legRho, ...
    'Position',[0.75 0.20 0.15 0.3], ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'Box','off');

%% SAVE
filename = fullfile(outDir,'figs_5g_uperp_tau_collapse.png');
print(gcf, filename, '-dpng', '-r600');

disp('τ-collapse plot saved.');