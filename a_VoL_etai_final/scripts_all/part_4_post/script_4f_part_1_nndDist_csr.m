clearvars; close all; clc;

N = 100;

ell = 1 / sqrt(N);        % nominal spacing
t_s = linspace(0, 2*ell, 141);

% Monte Carlo samples
B = 5000;

% storage
H_sim = zeros(B, numel(t_s));

%% Monte Carlo CSR (periodic unit box)
for b = 1:B
    % CSR points
    xs = rand(N,1);
    ys = rand(N,1);

    % nearest-neighbour distances
    dnn = zeros(N,1);

    for i = 1:N
        dx = xs - xs(i);
        dy = ys - ys(i);

        % periodic minimum image
        dx = dx - round(dx);
        dy = dy - round(dy);

        d2 = dx.^2 + dy.^2;
        d2(i) = inf;            % exclude self

        dnn(i) = sqrt(min(d2));
    end

    % empirical CDF of NND
    H_sim(b,:) = arrayfun(@(tt) mean(dnn <= tt), t_s);

    if mod(b,50) == 0
        fprintf('Monte Carlo %d / %d\n', b, B);
    end
end

%% CSR envelope
H_low  = quantile(H_sim, 0.025, 1);
H_high = quantile(H_sim, 0.975, 1);

%% Output directory
outDir = fullfile('data', 'data_4f_nndDist_csr');
if ~isfolder(outDir)
    mkdir(outDir);
end

outFile = fullfile(outDir, ...
    sprintf('csr_envelope_nnd_dist_periodic_N_%d_B_%d.txt', N, B));

%% Write file
fid = fopen(outFile, 'w');
fprintf(fid, '# CSR envelope for NND distribution (CDF) (periodic unit box)\n');
fprintf(fid, '# N=%d  B=%d  grid=141  tmax=2/sqrt(N)\n', N, B);
fprintf(fid, '# columns: t_s  H_low  H_high\n');

for k = 1:numel(t_s)
    fprintf(fid, '%.10f %.10f %.10f\n', ...
        t_s(k), H_low(k), H_high(k));
end

fclose(fid);

disp(['Saved: ' outFile]);
