#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <queue>
#include <limits>

namespace fs = std::filesystem;

struct Component {
    std::vector<int> nodes;
    double cx = 0.0;
    double cy = 0.0;
};

static Component choose_lcc(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double r,
                            bool has_prev,
                            double prev_cx,
                            double prev_cy)
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
    std::vector<Component> comps;

    for (int i = 0; i < N; ++i) {
        if (vis[i]) continue;

        std::queue<int> q;
        q.push(i);
        vis[i] = 1;

        Component c;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            c.nodes.push_back(u);
            c.cx += x[u];
            c.cy += y[u];

            for (int v : adj[u]) {
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }

        c.cx /= (double)c.nodes.size();
        c.cy /= (double)c.nodes.size();
        comps.push_back(c);
    }

    int best = -1;
    for (int i = 0; i < (int)comps.size(); ++i) {
        if (best == -1 || comps[i].nodes.size() > comps[best].nodes.size()) {
            best = i;
        } else if (comps[i].nodes.size() == comps[best].nodes.size()) {
            if (has_prev) {
                double di = (comps[i].cx - prev_cx) * (comps[i].cx - prev_cx)
                          + (comps[i].cy - prev_cy) * (comps[i].cy - prev_cy);
                double db = (comps[best].cx - prev_cx) * (comps[best].cx - prev_cx)
                          + (comps[best].cy - prev_cy) * (comps[best].cy - prev_cy);
                if (di < db) best = i;
            }
        }
    }

    return comps[best];
}

int main()
{
    const int N = 100;
    const double noise = 0.01;
    const int SIM_MIN = 1;
    const int SIM_MAX = 20;
    const int k_max = 20000;
    const double r = 1.0;

    std::vector<double> rho_vals = {1.00, 2.00, 4.00, 8.00};
    std::vector<double> speed_vals = {0.01, 0.02, 0.04, 0.08};

    fs::create_directories("../../part_3d_properties_data/part_3a_lcc_indices");

    for (int irho = 0; irho < (int)rho_vals.size(); ++irho)
    for (int ispd = 0; ispd < (int)speed_vals.size(); ++ispd) {
            const double rho = rho_vals[irho];
            const double spd = speed_vals[ispd];
            #pragma omp parallel for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim) {
            char outname[512];
            std::snprintf(outname, sizeof(outname),
                "../../part_3d_properties_data/part_3a_lcc_indices/"
                "lcc_indices_rho_%.2f_speed_%.2f_noise_%.2f_sim_%02d.txt",
                rho, spd, noise, sim);

            std::ofstream out(outname);
            out << "# snap_idx size indices_0based\n";

            bool has_prev = false;
            double prev_cx = 0.0, prev_cy = 0.0;

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

                Component best = choose_lcc(x, y, r, has_prev, prev_cx, prev_cy);
                has_prev = true;
                prev_cx = best.cx;
                prev_cy = best.cy;

                if ((int)best.nodes.size() < 5) continue;

                out << k << " " << best.nodes.size();
                for (int id : best.nodes) out << " " << id;
                out << "\n";
                break;
            }

            #pragma omp critical
            std::cout << "wrote " << outname << std::endl;
            }
    }

    std::cout << "DONE.\n";
    return 0;
}
