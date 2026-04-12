#ifndef M_PI
#define M_PI 3.1415926
#endif

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

struct Component {
    std::vector<int> nodes;
    double cx = 0.0;
    double cy = 0.0;
};

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
    const int total_time = 50000;
    const int save_every = 200;
    const int k_max = total_time / save_every;
    const double r = 1.0;

    std::vector<double> rho_vals = {1.00, 2.00, 4.00, 8.00};
    std::vector<double> speed_vals = {0.01, 0.02, 0.04, 0.08};

    fs::create_directories("../../part_3d_properties_data/part_3b_lcc_properties");

    for (double rho : rho_vals)
    for (double spd : speed_vals)
    {
        std::map<int, std::vector<double>> acc_size;
        std::map<int, std::vector<double>> acc_heading_cos;
        std::map<int, std::vector<double>> acc_heading_sin;
        std::map<int, std::vector<double>> acc_pol;
        std::map<int, std::vector<double>> acc_cx;
        std::map<int, std::vector<double>> acc_cy;

        #pragma omp parallel
        {
            std::map<int, std::vector<double>> loc_size;
            std::map<int, std::vector<double>> loc_heading_cos;
            std::map<int, std::vector<double>> loc_heading_sin;
            std::map<int, std::vector<double>> loc_pol;
            std::map<int, std::vector<double>> loc_cx;
            std::map<int, std::vector<double>> loc_cy;

            #pragma omp for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim) {
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

                    double sx = 0.0, sy = 0.0;
                    for (int id : best.nodes) {
                        sx += std::cos(th[id]);
                        sy += std::sin(th[id]);
                    }

                    int m = (int)best.nodes.size();
                    double heading = std::atan2(sy, sx);
                    double pol = std::sqrt(sx*sx + sy*sy) / (double)m;

                    loc_size[k].push_back((double)m);
                    loc_heading_cos[k].push_back(std::cos(heading));
                    loc_heading_sin[k].push_back(std::sin(heading));
                    loc_pol[k].push_back(pol);
                    loc_cx[k].push_back(best.cx);
                    loc_cy[k].push_back(best.cy);
                }
            }

            #pragma omp critical
            {
                for (auto& kv : loc_size)
                    acc_size[kv.first].insert(acc_size[kv.first].end(), kv.second.begin(), kv.second.end());
                for (auto& kv : loc_heading_cos)
                    acc_heading_cos[kv.first].insert(acc_heading_cos[kv.first].end(), kv.second.begin(), kv.second.end());
                for (auto& kv : loc_heading_sin)
                    acc_heading_sin[kv.first].insert(acc_heading_sin[kv.first].end(), kv.second.begin(), kv.second.end());
                for (auto& kv : loc_pol)
                    acc_pol[kv.first].insert(acc_pol[kv.first].end(), kv.second.begin(), kv.second.end());
                for (auto& kv : loc_cx)
                    acc_cx[kv.first].insert(acc_cx[kv.first].end(), kv.second.begin(), kv.second.end());
                for (auto& kv : loc_cy)
                    acc_cy[kv.first].insert(acc_cy[kv.first].end(), kv.second.begin(), kv.second.end());
            }
        }

        char outname[512];
        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3b_lcc_properties/"
            "lcc_properties_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);
        out << "# snap_idx mC_mean mC_std heading_mean heading_std pol_mean pol_std cx_mean cx_std cy_mean cy_std\n";

        for (const auto& kv : acc_size) {
            int k = kv.first;

            double m_mean, m_std;
            double c_mean, c_std;
            double s_mean, s_std;
            double p_mean, p_std;
            double cx_mean, cx_std;
            double cy_mean, cy_std;

            mean_std(acc_size[k], m_mean, m_std);
            mean_std(acc_heading_cos[k], c_mean, c_std);
            mean_std(acc_heading_sin[k], s_mean, s_std);
            mean_std(acc_pol[k], p_mean, p_std);
            mean_std(acc_cx[k], cx_mean, cx_std);
            mean_std(acc_cy[k], cy_mean, cy_std);

            double heading_mean = std::atan2(s_mean, c_mean);
            double R = std::sqrt(c_mean*c_mean + s_mean*s_mean);
            double heading_std = std::sqrt(std::max(0.0, -2.0 * std::log(std::max(1e-12, R))));

            out << k << " " << std::setprecision(10)
                << m_mean << " " << m_std << " "
                << heading_mean << " " << heading_std << " "
                << p_mean << " " << p_std << " "
                << cx_mean << " " << cx_std << " "
                << cy_mean << " " << cy_std << "\n";
        }

        std::cout << "wrote " << outname << std::endl;
    }

    std::cout << "DONE.\n";
    return 0;
}
