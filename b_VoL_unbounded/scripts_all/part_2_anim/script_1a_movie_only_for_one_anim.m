clear; clc; close all;
N = 100;
rho   = 1.0;
speed = 0.01;
noise = 0.01;
save_every = 200;

sim_idx = 1;
L_init = sqrt(N / rho);
view_box_length = 3.0 * L_init;

sim_folder = sprintf('../../sim_data/vicsek_unbounded_density_%.2f_speed_%.2f_noise_%.2f/sim_%02d/sim_data', ...
                     rho, speed, noise, sim_idx);

lcc_idx_file = sprintf(['../part_3d_properties_data/part_3a_lcc_indices/' ...
    'lcc_indices_rho_%.2f_speed_%.2f_noise_%.2f_sim_%02d.txt'], ...
    rho, speed, noise, sim_idx);

lcc_nodes_by_snap = load_lcc_nodes(lcc_idx_file);

files = dir(fullfile(sim_folder, sprintf('sim_data_noise_%.2f_sim_*_snap_*.txt', noise)));
if isempty(files)
    error('No simulation snapshot files found in: %s', sim_folder);
end

snap_idx = nan(1, numel(files));
for i = 1:numel(files)
    parsed_idx = sscanf(files(i).name, 'sim_data_noise_%*f_sim_%*d_snap_%d.txt');
    if ~isempty(parsed_idx)
        snap_idx(i) = parsed_idx;
    end
end

valid = ~isnan(snap_idx);
files = files(valid);
snap_idx = snap_idx(valid);
if isempty(files)
    error('Snapshot filenames did not match expected pattern in: %s', sim_folder);
end

[snap_idx, idx] = sort(snap_idx);
simDataFiles = files(idx);
numFrames = numel(simDataFiles);

video_name = sprintf('AniMotion_unbounded_N%d_rho_%.2f_speed_%.2f_sim_%02d', ...
                     N, rho, speed, sim_idx);

video_writer = VideoWriter(video_name, 'Motion JPEG AVI');
video_writer.FrameRate = 5;
open(video_writer);

figure_handle = figure('Position', [100, 100, 600, 600]);
set(figure_handle, 'Toolbar', 'none');

for i = 1:numFrames
    fname = fullfile(sim_folder, simDataFiles(i).name);
    data = readmatrix(fname);

    x = data(:, 1);
    y = data(:, 2);
    theta = data(:, 3);

    k = snap_idx(i);
    t = k * save_every;

    if (k + 1) <= numel(lcc_nodes_by_snap) && ~isempty(lcc_nodes_by_snap{k + 1})
        lcc_nodes = lcc_nodes_by_snap{k + 1};
        lcc_nodes = lcc_nodes(lcc_nodes >= 1 & lcc_nodes <= numel(x));
    else
        lcc_nodes = [];
    end

    if isempty(lcc_nodes)
        cx = mean(x);
        cy = mean(y);
    else
        cx = mean(x(lcc_nodes));
        cy = mean(y(lcc_nodes));
    end

    plt_motion(x, y, theta, t, figure_handle, rho, speed, cx, cy, view_box_length);
    frame = getframe(figure_handle);
    writeVideo(video_writer, frame);
end

close(video_writer);
fprintf('Saved: %s\n', video_name);


function lcc_nodes_by_snap = load_lcc_nodes(filename)
    lcc_nodes_by_snap = {};

    if ~isfile(filename)
        warning('LCC indices file not found: %s', filename);
        return;
    end

    fid = fopen(filename, 'r');
    if fid < 0
        warning('Could not open LCC indices file: %s', filename);
        return;
    end

    cleaner = onCleanup(@() fclose(fid));

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line, '#')
            continue;
        end

        vals = sscanf(line, '%d');
        if numel(vals) < 2
            continue;
        end

        k = vals(1);
        m = vals(2);

        if m > 0 && numel(vals) >= (2 + m)
            node_ids_1based = vals(3:2 + m) + 1;
        else
            node_ids_1based = [];
        end

        if numel(lcc_nodes_by_snap) < (k + 1)
            lcc_nodes_by_snap{k + 1} = [];
        end
        lcc_nodes_by_snap{k + 1} = node_ids_1based;
    end
end


function plt_motion(xx, yy, theta, t, figure_handle, rho, speed, cx, cy, box_length)
    figure(figure_handle);
    clf;
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

    u = cos(theta);
    v = sin(theta);

    quiver(xx, yy, u, v, 0.2, 'k');

    half_span = 0.5 * box_length;
    xlim([cx - half_span, cx + half_span]);
    ylim([cy - half_span, cy + half_span]);
    axis square;

    title(sprintf('$\\rho = %.2f,\\ v = %.2f,\\ t = %d$', rho, speed, t), ...
          'Interpreter', 'latex', 'FontSize', 20);

    box on;
end
