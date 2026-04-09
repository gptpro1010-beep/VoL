clearvars;close all;clc;

noise = 0.01;

rhos   = [1 4 16];
speeds = [0.01 0.04 0.16];

kvals = [1 2 5 10 20 50 100 200 500 1000];

dataDir = '../part_3d_properties_data/part_3d_ied';

envFile = fullfile('data', 'data_4d_ied_csr','csr_envelope_periodic_N_100_B_5000.txt');

outDir = fullfile('figs','figs_4d_part_2_ied_csr_test');
if ~isfolder(outDir)
    mkdir(outDir);
end


envDat = readmatrix(envFile);

t_s_env  = envDat(:,1);
H_low    = envDat(:,2);
H_high   = envDat(:,3);
H_mid    = 0.5*(H_low + H_high);


cc = parula(6);
cc1 = cc(1,:);
cc2 = cc(4,:);

line_colors = {cc1,cc2};


for iRho = 1:numel(rhos)

    rho = rhos(iRho);

    for iSpeed = 1:numel(speeds)

        s = speeds(iSpeed);

        for ik = 1:numel(kvals)

            snap_idx = kvals(ik);

            fullSubDir = fullfile(dataDir, ...
                sprintf('ied_avg_rho_%.2f_speed_%.2f_noise_%.2f', rho, s, noise));

            fpath = fullfile(fullSubDir, sprintf('k_%05d.txt', snap_idx));

            if ~isfile(fpath)
                warning('Missing file: %s', fpath);
                continue
            end

            dat = readmatrix(fpath);

            t_s_obs = dat(:,1);
            H_obs   = dat(:,2);


            % --- interpolate CSR envelope onto observation grid ---
            H_low_i  = interp1(t_s_env, H_low,  t_s_obs, 'linear', 'extrap');
            H_high_i = interp1(t_s_env, H_high, t_s_obs, 'linear', 'extrap');
            H_mid_i  = 0.5*(H_low_i + H_high_i);
            
            % --- Dmid (distance from CSR mean band center) ---
            Dmid = trapz(t_s_obs, abs(H_obs - H_mid_i));
            
            % --- deviation outside envelope ---
            d = zeros(size(H_obs));
            
            above = H_obs > H_high_i;
            below = H_obs < H_low_i;
            
            d(above) = H_obs(above) - H_high_i(above);
            d(below) = H_low_i(below) - H_obs(below);
            
            % --- Denv (distance outside envelope) ---
            Denv = trapz(t_s_obs, d);


            figure('Units','pixels','Position',[200 200 520 420]);
            ax = axes('Units','pixels','Position',[80 70 380 300]);
            hold(ax,'on');

            plot(ax, t_s_env, H_low,  '--', 'Color', line_colors{1}, 'LineWidth',1.5);
            plot(ax, t_s_env, H_high, '--', 'Color', line_colors{1}, 'LineWidth',1.5);
            plot(ax, t_s_env, H_obs,  '-',  'Color', line_colors{2}, 'LineWidth',2.0);

            box(ax,'on');
            grid(ax,'on');
            axis(ax,'square');

            xlabel(ax,'$t_s$','Interpreter','latex','FontSize',16);
            ylabel(ax,'$H(t_s)$','Interpreter','latex','FontSize',16);

            ax.FontSize = 16;
            ax.TickLabelInterpreter = 'latex';

            xlim(ax,[t_s_env(1) t_s_env(end)]);
            ylim(ax,[0 1]);


            txt = sprintf('$D_{mid}=%.4f$\n$D_{env}=%.4f$', Dmid, Denv);

            text(ax,0.05,0.95,txt, ...
                'Units','normalized', ...
                'VerticalAlignment','top', ...
                'HorizontalAlignment','left', ...
                'Interpreter','latex', ...
                'FontSize',16, ...
                'BackgroundColor','w', ...
                'EdgeColor','k');


                txt2 = sprintf('$\\rho=%.2f$\n$s=%.2f$\n$T/T_{cross}=%d$', ...
                    rho, s, 5*snap_idx);

            text(ax,0.95,0.05,txt2, ...
                'Units','normalized', ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right', ...
                'Interpreter','latex', ...
                'FontSize',16, ...
                'BackgroundColor','w', ...
                'EdgeColor','k');


                outname = sprintf('csr_test_rho_%.2f_speed_%.2f_snap_%05d.png', ...
                    rho, s, snap_idx);

            print(gcf, fullfile(outDir,outname), '-dpng','-r600');

            close(gcf);

        end

    end

end


disp('CSR test plots generated.');