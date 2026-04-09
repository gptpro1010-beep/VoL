clearvars;
close all;
clc;

N = 100;

t_s = linspace(0, sqrt(2)/2, 141);

B = 5000;
H_sim = zeros(B, numel(t_s));

M = N*(N-1)/2;

for b = 1:B
    xs = rand(N,1);
    ys = rand(N,1);

    Dsim = zeros(M,1);
    idx = 1;

    for i = 1:N-1
        for j = i+1:N
            dx = xs(j) - xs(i);
            dy = ys(j) - ys(i);

            dx = dx - round(dx);
            dy = dy - round(dy);

            Dsim(idx) = sqrt(dx^2 + dy^2);
            idx = idx + 1;
        end
    end

    H_sim(b,:) = arrayfun(@(tt) mean(Dsim <= tt), t_s);

    if mod(b,50) == 0
        fprintf('Monte Carlo %d / %d\n', b, B);
    end
end

H_low  = quantile(H_sim, 0.025, 1);
H_high = quantile(H_sim, 0.975, 1);

outDir = fullfile('data', 'data_4d_ied_csr');
if ~isfolder(outDir)
    mkdir(outDir);
end

outFile = fullfile(outDir, sprintf('csr_envelope_periodic_N_%d_B_%d.txt', N, B));

fid = fopen(outFile, 'w');
fprintf(fid, '# CSR envelope (periodic unit box)\n');
fprintf(fid, '# N=%d  B=%d  grid=141  tmax=sqrt(2)/2\n', N, B);
fprintf(fid, '# columns: t_s  H_low  H_high\n');

for k = 1:numel(t_s)
    fprintf(fid, '%.10f %.10f %.10f\n', t_s(k), H_low(k), H_high(k));
end

fclose(fid);

disp(['Saved: ' outFile]);
