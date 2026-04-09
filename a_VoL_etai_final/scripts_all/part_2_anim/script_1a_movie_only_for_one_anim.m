clear; clc; close all;
N = 100;
rho   = 1.0;
speed = 0.01;
L = sqrt(N / rho);

sim_idx = 1;
sim_folder = sprintf('../../sim_data/vicsek_density_%.2f_speed_%.2f_noise_0.01/sim_%02d/sim_data', rho, speed, sim_idx);

files = dir(fullfile(sim_folder, 'sim_data_noise_0.01_sim_*_snap_*.txt'));

snap_idx = zeros(1, numel(files));
for i = 1:numel(files)
    snap_idx(i) = sscanf(files(i).name, 'sim_data_noise_%*f_sim_%*d_snap_%d.txt');
end

desiredSteps = 1:50;

[snap_idx, ia] = intersect(snap_idx, desiredSteps);
files = files(ia);

[snap_idx, idx] = sort(snap_idx);
simDataFiles = files(idx);

numFrames = numel(simDataFiles);

video_name = sprintf('AniMotion_N%d_rho_%.2f_speed_%.2f_sim_%02d', ...
                     N, rho, speed, sim_idx);

video_writer = VideoWriter(video_name, 'Motion JPEG AVI');
video_writer.FrameRate = 5;
open(video_writer);

figure_handle = figure('Position', [100, 100, 600, 600]);
set(figure_handle, 'Toolbar', 'none');

% --- Generate Video ---
for i = 1:numFrames
    fname = fullfile(sim_folder, simDataFiles(i).name);
    data = readmatrix(fname);

    x = data(:, 1);
    y = data(:, 2);
    theta = data(:, 3);

    plt_motion(x, y, theta, snap_idx(i), L, figure_handle, rho, speed);
    frame = getframe(gcf);
    writeVideo(video_writer, frame);
end

close(video_writer);
fprintf('Saved: %s\n', video_name);


function plt_motion(xx, yy, theta, t, L, figure_handle, rho, speed)
    figure(figure_handle);
    clf;
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

    u = cos(theta);
    v = sin(theta);

    % scatter(xx, yy, 60, 'k', 'filled'); hold on;
    quiver(xx, yy, u, v, 0.2, 'k');

    xlim([0 L]);
    ylim([0 L]);
    axis square;

    title(sprintf('$\\rho = %.2f,\\ v = %.2f,\\ t = %d$', rho, speed, t), ...
          'Interpreter', 'latex', 'FontSize', 20);

    % set(gca,'XTick',[],'YTick',[]);
    box on;
end
