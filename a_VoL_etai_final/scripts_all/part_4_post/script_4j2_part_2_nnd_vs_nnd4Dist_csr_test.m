clearvars;
close all;
clc;

%% Parameters
noise = 0.01;

rhos   = [1.00 4.00 16.00];
speeds = [0.01 0.04 0.16];

kvals = [1 5 10 20 50 100 200 500 1000];

nndDir  = '../part_3d_properties_data/part_3f_nnd_dist';
nnd4Dir = '../part_3d_properties_data/part_3j_nnd4_dist';

envNndFile = fullfile('data', 'data_4f_nndDist_csr', ...
    'csr_envelope_nnd_dist_periodic_N_100_B_5000.txt');

envNnd4File = fullfile('data', 'data_4j_nnd4Dist_csr', ...
    'csr_envelope_nnd4_dist_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4j2_part_2_nnd_vs_nnd4Dist_csr_test');
if ~isfolder(outDir), mkdir(outDir); end

%% Load CSR envelopes
envNnd  = readmatrix(envNndFile);
envNnd4 = readmatrix(envNnd4File);

t_s_env_nnd   = envNnd(:,1);
F_nnd_low     = envNnd(:,2);
F_nnd_high    = envNnd(:,3);

t_s_env_nnd4  = envNnd4(:,1);
F_nnd4_low    = envNnd4(:,2);
F_nnd4_high   = envNnd4(:,3);

%% Colors
cc = parula(6);
col_nnd      = cc(2,:);
col_nnd4     = cc(5,:);
col_env_nnd  = cc(1,:);
col_env_nnd4 = cc(4,:);

%% Loop over parameters
for iRho = 1:numel(rhos)

    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        %% Load averaged NND file
        fNnd = fullfile(nndDir, ...
            sprintf('nnd_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                    rho, s, noise));

        %% Load averaged NND4 file
        fNnd4 = fullfile(nnd4Dir, ...
            sprintf('nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt', ...
                    rho, s, noise));

        if ~isfile(fNnd)
            warning('Missing NND file: %s', fNnd);
            continue;
        end

        if ~isfile(fNnd4)
            warning('Missing NND4 file: %s', fNnd4);
            continue;
        end

        datNnd  = readmatrix(fNnd);
        datNnd4 = readmatrix(fNnd4);

        if isempty(datNnd) || size(datNnd,2) < 4
            warning('Bad NND format: %s', fNnd);
            continue;
        end

        if isempty(datNnd4) || size(datNnd4,2) < 4
            warning('Bad NND4 format: %s', fNnd4);
            continue;
        end

        snap_nnd = datNnd(:,1);
        t_nnd    = datNnd(:,2);
        F_nnd    = datNnd(:,3);
        S_nnd    = datNnd(:,4);

        snap_nnd4 = datNnd4(:,1);
        t_nnd4    = datNnd4(:,2);
        F_nnd4    = datNnd4(:,3);
        S_nnd4    = datNnd4(:,4);

        %% Loop over snapshots
        for ik = 1:numel(kvals)

            snap_idx = kvals(ik);

            maskNnd  = (snap_nnd  == snap_idx);
            maskNnd4 = (snap_nnd4 == snap_idx);

            if ~any(maskNnd) || ~any(maskNnd4)
                warning('Missing snapshot=%d for rho=%.2f, s=%.2f', ...
                    snap_idx, rho, s);
                continue;
            end

            %% Extract and sort
            tObsNnd  = t_nnd(maskNnd);
            FObsNnd  = F_nnd(maskNnd);
            SObsNnd  = S_nnd(maskNnd);
            [tObsNnd, idx1] = sort(tObsNnd);
            FObsNnd = FObsNnd(idx1);
            SObsNnd = SObsNnd(idx1);

            tObsNnd4 = t_nnd4(maskNnd4);
            FObsNnd4 = F_nnd4(maskNnd4);
            SObsNnd4 = S_nnd4(maskNnd4);
            [tObsNnd4, idx2] = sort(tObsNnd4);
            FObsNnd4 = FObsNnd4(idx2);
            SObsNnd4 = SObsNnd4(idx2);

            %% Plot
            figure('Units','pixels','Position',[200 200 620 440]);
            ax = axes('Units','pixels','Position',[90 70 470 320]);
            hold(ax,'on');

            % CSR envelopes
            plot(ax, t_s_env_nnd,  F_nnd_low,  '--', 'Color', col_env_nnd,  'LineWidth', 1.3);
            plot(ax, t_s_env_nnd,  F_nnd_high, '--', 'Color', col_env_nnd,  'LineWidth', 1.3);
            plot(ax, t_s_env_nnd4, F_nnd4_low, '--', 'Color', col_env_nnd4, 'LineWidth', 1.3);
            plot(ax, t_s_env_nnd4, F_nnd4_high,'--', 'Color', col_env_nnd4, 'LineWidth', 1.3);

            % Main curves
            hNnd  = plot(ax, tObsNnd,  FObsNnd,  '-', 'Color', col_nnd,  'LineWidth', 2.0);
            hNnd4 = plot(ax, tObsNnd4, FObsNnd4, '-', 'Color', col_nnd4, 'LineWidth', 2.0);

            % Standard deviation shading (keep commented)
            % fill(ax, [tObsNnd; flipud(tObsNnd)], [FObsNnd-SObsNnd; flipud(FObsNnd+SObsNnd)], col_nnd,  'FaceAlpha',0.20, 'EdgeColor','none');
            % fill(ax, [tObsNnd4; flipud(tObsNnd4)], [FObsNnd4-SObsNnd4; flipud(FObsNnd4+SObsNnd4)], col_nnd4, 'FaceAlpha',0.20, 'EdgeColor','none');

            box(ax,'on');
            grid(ax,'on');
            axis(ax,'square');

            xlabel('$t_s$', 'Interpreter','latex','FontSize',16);
            ylabel('$H(t_s)$ and $H_4(t_s)$', 'Interpreter','latex','FontSize',16);

            ax.FontSize = 15;
            ax.TickLabelInterpreter = 'latex';
            xlim(ax,[0 max([tObsNnd(:); tObsNnd4(:)])]);
            ylim(ax,[0 1]);

            lg = legend(ax, [hNnd hNnd4], {'NND', 'NND4'}, ...
                'Location','southeast', ...
                'Interpreter','latex', ...
                'FontSize',14, ...
                'Box','on');
            lg.Color = 'w';

            txt = sprintf('$\\rho = %.2f$\\n$s = %.2f$\\n$T/T_{cross} = %d$', ...
                rho, s, 5*snap_idx);
            text(ax,0.04,0.95,txt, ...
                'Units','normalized', ...
                'VerticalAlignment','top', ...
                'HorizontalAlignment','left', ...
                'Interpreter','latex', ...
                'FontSize',14, ...
                'BackgroundColor','w', ...
                'EdgeColor','k');

            outname = sprintf('csr_test_nnd_vs_nnd4_rho_%.2f_speed_%.2f_snap_%05d.png', ...
                rho, s, snap_idx);
            print(gcf, fullfile(outDir, outname), '-dpng', '-r600');
            close(gcf);

        end
    end
end

disp('NND vs NND4 CSR comparison plots generated.');
