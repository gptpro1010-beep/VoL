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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;


/* ============================================================
   Compute mean heading and order parameter R
   ============================================================ */
void compute_mean_heading_and_R(const std::vector<float>& theta,
                                float& theta_mean,
                                float& R)
{
    int N = (int)theta.size();

    float sx = 0.0f, sy = 0.0f;

    for (float th : theta) {
        sx += std::cos(th);
        sy += std::sin(th);
    }

    theta_mean = std::atan2(sy, sx);
    R = std::sqrt(sx*sx + sy*sy) / N;
}


/* ============================================================
   Mean + std
   ============================================================ */
void mean_std(const std::vector<double>& v,
              double& mean,
              double& stddev)
{
    int n = (int)v.size();

    if (n == 0) {
        mean = stddev = 0.0;
        return;
    }

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


/* ============================================================
   MAIN
   ============================================================ */
int main()
{
    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;

    fs::create_directories("../../part_3d_properties_data/part_3h_avgheading");

    std::vector<float> rho_vals   = {1.00f, 2.00f, 4.00f, 8.00f, 16.00f};
    std::vector<float> speed_vals = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};


    for (float rho : rho_vals)
    for (float spd : speed_vals)
    {
        float L = std::sqrt((float)N / rho);

        int T      = (int)std::round(L / spd);
        int T_max  = 5000 * T;
        int dt_out = 5 * T;

        int k_max = T_max / dt_out;

        std::map<int, std::vector<double>> accum_cos;
        std::map<int, std::vector<double>> accum_sin;
        std::map<int, std::vector<double>> accum_R;

        std::cout << "rho=" << std::fixed << std::setprecision(2)
                  << rho << "  speed=" << spd << std::endl;


        /* ============================================================
           Parallel over simulations
           ============================================================ */
        #pragma omp parallel
        {
            std::vector<float> theta(N);

            std::map<int, std::vector<double>> local_cos;
            std::map<int, std::vector<double>> local_sin;
            std::map<int, std::vector<double>> local_R;


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

                    while (cnt < N && (in >> x >> y >> th)) {
                        theta[cnt] = th;
                        cnt++;
                    }

                    if (cnt != N) continue;

                    float theta_mean, R;
                    compute_mean_heading_and_R(theta, theta_mean, R);

                    local_cos[k].push_back(std::cos(theta_mean));
                    local_sin[k].push_back(std::sin(theta_mean));
                    local_R[k].push_back(R);
                }
            }


            /* ============================================================
               Merge thread-local data
               ============================================================ */
            #pragma omp critical
            {
                for (auto& kv : local_cos)
                    accum_cos[kv.first].insert(
                        accum_cos[kv.first].end(),
                        kv.second.begin(), kv.second.end());

                for (auto& kv : local_sin)
                    accum_sin[kv.first].insert(
                        accum_sin[kv.first].end(),
                        kv.second.begin(), kv.second.end());

                for (auto& kv : local_R)
                    accum_R[kv.first].insert(
                        accum_R[kv.first].end(),
                        kv.second.begin(), kv.second.end());
            }

        } // end omp


        if (accum_cos.empty()) continue;


        /* ============================================================
           Write output
           ============================================================ */
        char outname[256];

        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3h_avgheading/"
            "avgheading_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);

        out << "# mean heading and R across simulations\n";
        out << "# rho=" << std::fixed << std::setprecision(2)
            << rho << " speed=" << spd
            << " noise=" << noise << "\n";

        out << "# snap_idx  theta_mean  theta_std  R_mean  R_std\n";


        for (const auto& kv : accum_cos)
        {
            int k = kv.first;

            double mean_cos, std_cos;
            double mean_sin, std_sin;
            double mean_R, std_R;

            mean_std(accum_cos[k], mean_cos, std_cos);
            mean_std(accum_sin[k], mean_sin, std_sin);
            mean_std(accum_R[k], mean_R, std_R);

            double theta_mean = std::atan2(mean_sin, mean_cos);

            // optional: angular std approximation
            double theta_std = std::sqrt(-2.0 * std::log(std::max(1e-12, mean_R)));

            out << k << " "
                << std::setprecision(10)
                << theta_mean << " "
                << theta_std << " "
                << mean_R << " "
                << std_R << "\n";
        }

        std::cout << "  wrote " << outname << std::endl;
    }


    std::cout << "DONE.\n";
    return 0;
}