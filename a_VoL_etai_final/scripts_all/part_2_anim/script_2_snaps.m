clearvars; close all; clc;

%% Parameters
N     = 100;
noise = 0.01;

rhos   = [1 4 16];
speeds = [0.01 0.04 0.16];

% rhos   = [7];
% speeds = [0.09];

% Snapshot indices (snap_idx)
snap_idxs = [1 2 5 10 15 25 50 100 200 500 1000];

sim_idx = 1;

%% Output directory
outRoot = fullfile('figs','frames_quiver');
if ~isfolder(outRoot), mkdir(outRoot); end

%% Data path format
dataRootFmt = fullfile('..','..','sim_data', ...
    'vicsek_density_%.2f_speed_%.2f_noise_%.2f', ...
    sprintf('sim_%02d',sim_idx), ...
    'sim_data');

%% Figure
fig = figure('Units','pixels','Position',[100 100 650 650],'Color','w');
ax  = axes(fig);
set(ax,'FontSize',18,'TickLabelInterpreter','latex');
set(fig,'Renderer','painters');

%% Loop
for iR = 1:numel(rhos)

    rho = rhos(iR);
    L = sqrt(N / rho);

    for iS = 1:numel(speeds)

        spd = speeds(iS);

        sim_folder = sprintf(dataRootFmt,rho,spd,noise);

        outDir = fullfile(outRoot, ...
            sprintf('rho_%.2f',rho), ...
            sprintf('speed_%.2f',spd), ...
            sprintf('sim_%02d',sim_idx));

        if ~isfolder(outDir), mkdir(outDir); end

        fprintf('rho=%.2f speed=%.2f\n',rho,spd);

        for it = 1:numel(snap_idxs)

            k = snap_idxs(it);

            fname = fullfile(sim_folder, ...
                sprintf('sim_data_noise_%.2f_sim_%02d_snap_%05d.txt', ...
                noise,sim_idx,k));

            data = readmatrix(fname);

            x  = data(1:N,1);
            y  = data(1:N,2);
            th = data(1:N,3);

            cla(ax)

            u = cos(th);
            v = sin(th);

            quiver(ax,x,y,u,v,0.2,'k')

            xlim(ax,[0 L])
            ylim(ax,[0 L])
            axis(ax,'square')

            title(ax, ...
                sprintf('$\\rho=%.2f,\\ v=%.2f,\\ \\eta=%.2f,\\ k=%d$', ...
                rho,spd,noise,k), ...
                'Interpreter','latex','FontSize',18)

            set(ax,'XTick',[],'YTick',[])
            box(ax,'on')

            outPng = fullfile(outDir, ...
                sprintf('frame_snap_%05d.png',k));

            exportgraphics(fig,outPng,'Resolution',600)

            fprintf('  saved snap=%d\n',k)

        end
    end
end

disp('DONE')