#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <queue>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

void mean_std(const std::vector<double>& v, double& mean, double& stddev)
{
    int n = (int)v.size();
    if (n == 0) { mean = stddev = 0.0; return; }

    double s = 0.0;
    for (double x : v) s += x;
    mean = s / n;

    double var = 0.0;
    for (double x : v) {
        double d = x - mean;
        var += d*d;
    }
    stddev = (n > 1) ? std::sqrt(var/(n-1)) : 0.0;
}

static int count_single_particles(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  double r)
{
    const int N = (int)x.size();
    const double r2 = r * r;

    std::vector<std::vector<int>> adj(N);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double d2 = dx*dx + dy*dy;
            if (d2 <= r2) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    std::vector<int> vis(N, 0);
    int singles = 0;

    for (int i = 0; i < N; ++i) {
        if (vis[i]) continue;

        std::queue<int> q;
        q.push(i);
        vis[i] = 1;

        int comp_size = 0;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            comp_size++;

            for (int v : adj[u]) {
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }

        if (comp_size == 1) singles++;
    }

    return singles;
}

int main()
{
    const int N = 100;
    const double noise = 0.01;
    const int SIM_MIN = 1;
    const int SIM_MAX = 20;
    const int total_time = 50000;
    const int save_every = 200;
    const int k_max = total_time / save_every;
    const double r = 1.0;

    std::vector<double> rho_vals = {1.00, 2.00, 4.00, 8.00};
    std::vector<double> speed_vals = {0.01, 0.02, 0.04, 0.08};

    fs::create_directories("../../part_3d_properties_data/part_3c_single_particles");

    for (double rho : rho_vals)
    for (double spd : speed_vals)
    {
        std::map<int, std::vector<double>> acc_singles;

        #pragma omp parallel
        {
            std::map<int, std::vector<double>> loc_singles;

            #pragma omp for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim) {
                for (int k = 0; k <= k_max; ++k) {
                    char fpath[1024];
                    std::snprintf(fpath, sizeof(fpath),
                        "../../../sim_data/vicsek_unbounded_density_%.2f_speed_%.2f_noise_%.2f/"
                        "sim_%02d/sim_data/"
                        "sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                        rho, spd, noise,
                        sim,
                        noise, sim, k);

                    std::ifstream in(fpath);
                    if (!in) continue;

                    std::vector<double> x(N), y(N), th(N);
                    int cnt = 0;
                    while (cnt < N && (in >> x[cnt] >> y[cnt] >> th[cnt])) cnt++;
                    if (cnt != N) continue;

                    int singles = count_single_particles(x, y, r);
                    loc_singles[k].push_back((double)singles);
                }
            }

            #pragma omp critical
            {
                for (auto& kv : loc_singles)
                    acc_singles[kv.first].insert(acc_singles[kv.first].end(), kv.second.begin(), kv.second.end());
            }
        }

        char outname[512];
        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3c_single_particles/"
            "single_particles_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);
        out << "# snap_idx singles_mean singles_std\n";

        for (const auto& kv : acc_singles) {
            int k = kv.first;
            double m, s;
            mean_std(acc_singles[k], m, s);
            out << k << " " << std::setprecision(10) << m << " " << s << "\n";
        }

        std::cout << "wrote " << outname << std::endl;
    }

    std::cout << "DONE.\n";
    return 0;
}
