#ifndef M_PI
#define M_PI 3.1415926
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <map>
#include <limits>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;


/* ============================================================
   Compute NND4 for one snapshot
   (for each particle: average distance to its 4 nearest neighbors)
   ============================================================ */
void compute_nnd4(const std::vector<float>& xx,
                  const std::vector<float>& yy,
                  float L,
                  std::vector<float>& nnd4)
{
    const int N = (int)xx.size();
    nnd4.assign(N, 0.0f);

    if (N <= 1) return;

    const int k_use = std::min(4, N - 1);

    std::vector<float> d2;
    d2.reserve(std::max(0, N - 1));

    for (int i = 0; i < N; ++i)
    {
        d2.clear();

        const float xi = xx[i];
        const float yi = yy[i];

        for (int j = 0; j < N; ++j)
        {
            if (i == j) continue;

            float dx = xx[j] - xi;
            float dy = yy[j] - yi;

            dx -= L * std::round(dx / L);
            dy -= L * std::round(dy / L);

            d2.push_back(dx * dx + dy * dy);
        }

        std::partial_sort(d2.begin(), d2.begin() + k_use, d2.end());

        double sum4 = 0.0;

        for (int k = 0; k < k_use; ++k)
            sum4 += std::sqrt(d2[k]);

        nnd4[i] = (float)(sum4 / k_use);
    }
}


/* ============================================================
   Mean + std across simulations
   ============================================================ */
void mean_std(const std::vector<double>& v,
              double& mean,
              double& stddev)
{
    const int n = (int)v.size();

    if (n == 0)
    {
        mean = stddev = 0.0;
        return;
    }

    double s = 0.0;

    for (double x : v)
        s += x;

    mean = s / n;

    double var = 0.0;

    for (double x : v)
    {
        const double d = x - mean;
        var += d * d;
    }

    stddev = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;
}


/* ============================================================
   MAIN
   ============================================================ */
int main()
{
    const int   N       = 100;
    const float noise   = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;

    fs::create_directories("../../part_3d_properties_data/part_3i_nnd4");

    std::vector<float> rho_vals   = {1.00f, 2.00f, 4.00f, 8.00f, 16.00f};
    std::vector<float> speed_vals = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};

    for (float rho : rho_vals)
    for (float spd : speed_vals)
    {
        const float L = std::sqrt((float)N / rho);

        const int T      = (int)std::round(L / spd);
        const int T_max  = 5000 * T;
        const int dt_out = 5 * T;

        const int k_max = T_max / dt_out;

        std::map<int, std::vector<double>> accum;

        std::cout << "rho=" << std::fixed << std::setprecision(2)
                  << rho << "  speed=" << spd << std::endl;


        /* ============================================================
           Parallel over simulations
           ============================================================ */

        #pragma omp parallel
        {
            std::vector<float> xx(N), yy(N), nnd4;
            std::map<int, std::vector<double>> local_accum;

            #pragma omp for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim)
            {
                for (int k = 0; k <= k_max; ++k)
                {
                    char fpath[1024];

                    std::snprintf(fpath, sizeof(fpath),
                        "../../../sim_data/vicsek_density_%.2f_speed_%.2f_noise_%.2f/"
                        "sim_%02d/sim_data/"
                        "sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                        rho, spd, noise,
                        sim,
                        noise, sim, k);

                    std::ifstream in(fpath);

                    if (!in) continue;

                    float x, y, th;
                    int cnt = 0;

                    while (cnt < N && (in >> x >> y >> th))
                    {
                        xx[cnt] = x;
                        yy[cnt] = y;
                        cnt++;
                    }

                    if (cnt != N) continue;

                    compute_nnd4(xx, yy, L, nnd4);

                    double avg = 0.0;

                    for (float v : nnd4)
                        avg += v;

                    avg /= N;

                    local_accum[k].push_back(avg);
                }
            }

            /* ============================================================
               Merge thread-local data
               ============================================================ */

            #pragma omp critical
            {
                for (auto& kv : local_accum)
                {
                    auto& dst = accum[kv.first];

                    dst.insert(dst.end(),
                               kv.second.begin(),
                               kv.second.end());
                }
            }

        } // end omp


        if (accum.empty()) continue;


        /* ============================================================
           Write output
           ============================================================ */

        char outname[256];

        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3i_nnd4/"
            "nnd4_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);

        out << "# avg NND4(k): mean +- std across simulations\n";
        out << "# rho=" << std::fixed << std::setprecision(2)
            << rho << " speed=" << spd
            << " noise=" << noise << "\n";

        out << "# snap_idx  mean  std\n";

        for (const auto& kv : accum)
        {
            double m, sd;
            mean_std(kv.second, m, sd);

            out << kv.first << " "
                << std::setprecision(10) << m << " "
                << sd << "\n";
        }

        std::cout << "  wrote " << outname << std::endl;
    }

    std::cout << "DONE.\n";

    return 0;
}
