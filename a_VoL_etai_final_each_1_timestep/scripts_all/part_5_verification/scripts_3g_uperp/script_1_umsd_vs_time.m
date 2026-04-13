clearvars; close all; clc;

%% PARAMETERS
N = 100;
noise  = 0.01;

rhos   = [0.25 8];
% speeds = [0.01 0.02 0.04 0.08 0.16];
speeds = [0.01 0.04 0.16];

dataDir = '../../part_3d_properties_data/part_3g_uperp';

outDir = fullfile('figs','figs_5g_uperp');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% PLOT

figure('Units','pixels','Position',[200 200 820 420]);
ax = axes('Units','pixels','Position',[90 80 680 300]);
hold(ax,'on');

cc = parula(numel(speeds)+2);
markers = {'o','s','d','^','v'};

for iSpeed = 1:numel(speeds)

    s = speeds(iSpeed);

    for iRho = 1:numel(rhos)

        rho = rhos(iRho);

        fname = sprintf( ...
            'uperp_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
            rho, s, noise);

        fpath = fullfile(dataDir, fname);

        if ~isfile(fpath)
            % continue
        end

        dat = readmatrix(fpath);

        k    = dat(:,1);
        u_msd = dat(:,2);
        u_std = dat(:,3);

        
        %% remove t = 0
        valid = k > 0;

        k    = k(valid);
        u_msd = u_msd(valid);
        u_std = u_std(valid);

        %% time
        T = round(sqrt(N/rho) / s);   % same as simulation
        dt_out = 1;
        
        t = k * dt_out;

        %% plot
        loglog(ax, t, u_msd, ...
            'LineStyle','none', ...
            'Color', cc(iSpeed,:), ...
            'Marker', markers{iRho}, ...
            'MarkerSize', 6);
            
        %% shaded std band (mean ± std)
        upper = u_msd + u_std;
        lower = u_msd - u_std;

        % avoid negative values (important for log scale)
        lower(lower <= 0) = min(u_msd(u_msd>0)) * 0.1;

        fill([t; flipud(t)], ...
             [upper; flipud(lower)], ...
             cc(iSpeed,:), ...
             'FaceAlpha', 0.15, ...
             'EdgeColor', 'none');
        
        %% scaling
        p = polyfit(log10(t), log10(u_msd), 1);
        slope = p(1);
        ycut  = p(2);
        
        % predicted values
        y_fit = polyval(p, log10(t));
        
        % actual values
        y_actual = log10(u_msd);
        
        % R^2 calculation
        SS_res = sum((y_actual - y_fit).^2);
        SS_tot = sum((y_actual - mean(y_actual)).^2);
        R2 = 1 - SS_res / SS_tot;
        
        fprintf('s = %.2f, rho = %.2f: slope = %.3f, R^2 = %.4f\n', s, rho, slope, R2);
        fprintf('s = %.2f, rho = %.2f: ycut = %.3f\n', s, rho, ycut);

    end
end

%% AXIS
set(ax, 'XScale', 'log', 'YScale', 'log');

ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$\langle u_\perp^2 \rangle $', 'Interpreter','latex','FontSize',18);
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
    legRho{iRho} = sprintf('$\\rho = %.2f$', rhos(iRho));
end

legend(ax2, hRho, legRho, ...
    'Position',[0.75 0.20 0.15 0.3], ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'Box','off');

%% SAVE

filename = fullfile(outDir,'figs_5g_uperp_loglog.png');
print(gcf, filename, '-dpng', '-r600');

disp('Log-log plot saved.');