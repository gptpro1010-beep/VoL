clearvars; close all; clc;

%% PARAMETERS
N = 100;
noise = 0.01;

rho = 1.00;
speed = 0.04;
sim_nums = 1:5;

dataDir = '../../part_3d_properties_data/part_3g_uperp_each_particles_1to5';
outDir = fullfile('figs','figs_5g_uperp_each_particles_1to5');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% PLOT ONE FIGURE PER SIM (one line per particle)
cc = parula(N);

for iSim = 1:numel(sim_nums)

    sim = sim_nums(iSim);

    fname = sprintf('uperp2_particles_rho_%.2f_speed_%.2f_noise_%.2f_sim_%02d.txt', ...
        rho, speed, noise, sim);
    fpath = fullfile(dataDir, fname);

    if ~isfile(fpath)
        fprintf('Missing file: %s\n', fpath);
        continue;
    end

    dat = readmatrix(fpath);

    k = dat(:,1);
    u2_particles = dat(:,2:end);

    % optional: remove t=0
    valid = k >= 0;
    k = k(valid);
    u2_particles = u2_particles(valid,:);

    dt_out = 1;
    t = k * dt_out;

    fig = figure('Units','pixels','Position',[120 120 900 520]);
    ax = axes('Units','pixels','Position',[90 80 760 380]);
    hold(ax,'on');

    for iParticle = 1:N
        plot(ax, t, u2_particles(:,iParticle), ...
            'Color', cc(iParticle,:), ...
            'LineWidth', 0.9);
    end

    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';
    xlabel('$t$','Interpreter','latex','FontSize',18);
    ylabel('$u_{\perp,i}^2(t)$','Interpreter','latex','FontSize',18);
    title(sprintf('$\\rho=%.2f,\\ s=%.2f,\\ \eta=%.2f,\\ sim=%02d$',rho,speed,noise,sim), ...
        'Interpreter','latex','FontSize',16);
    grid(ax,'on');
    box(ax,'on');

    colormap(ax, cc);
    cb = colorbar(ax);
    cb.Label.String = 'particle index';
    cb.Label.Interpreter = 'latex';
    cb.FontSize = 12;

    outPng = fullfile(outDir, sprintf('uperp2_particles_lines_rho_%.2f_speed_%.2f_sim_%02d.png',rho,speed,sim));
    print(fig, outPng, '-dpng', '-r500');

    close(fig);
    fprintf('Saved %s\n', outPng);
end

fprintf('DONE plotting particle-wise u_perp^2 lines for sim 1..5.\n');
