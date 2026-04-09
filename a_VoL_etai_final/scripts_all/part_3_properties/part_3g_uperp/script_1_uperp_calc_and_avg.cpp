#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <map>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

/* ============================================================
   Minimum image
   ============================================================ */
static inline float min_image(float dx, float L)
{
    dx -= L * std::round(dx / L);
    return dx;
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
        mean = 0.0;
        stddev = 0.0;
        return;
    }

    double s = 0.0;
    for (double x : v) s += x;
    mean = s / n;

    double var = 0.0;
    for (double x : v) {
        double d = x - mean;
        var += d * d;
    }

    stddev = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;
}

/* ============================================================
   Read snapshot
   ============================================================ */
bool read_snapshot(const char* fpath,
                   int N,
                   std::vector<float>& x,
                   std::vector<float>& y,
                   std::vector<float>& theta)
{
    std::ifstream in(fpath);
    if (!in) return false;

    x.resize(N);
    y.resize(N);
    theta.resize(N);

    float xx, yy, th;
    int cnt = 0;

    while (cnt < N && (in >> xx >> yy >> th))
    {
        x[cnt] = xx;
        y[cnt] = yy;
        theta[cnt] = th;
        cnt++;
    }

    return (cnt == N);
}

/* ============================================================
   Compute mean direction (instantaneous)
   ============================================================ */
void compute_mean_direction(const std::vector<float>& theta,
                            float& ex,
                            float& ey)
{
    float sx = 0.0f;
    float sy = 0.0f;

    for (float th : theta) {
        sx += std::cos(th);
        sy += std::sin(th);
    }

    float norm = std::sqrt(sx * sx + sy * sy) + 1e-8f;
    ex = sx / norm;
    ey = sy / norm;
}

/* ============================================================
   Incremental u_perp^2 (rotating frame, correct drift removal)
   ============================================================ */
double compute_uperp2_incremental(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& theta,
    const std::vector<float>& x_prev,
    const std::vector<float>& y_prev,
    std::vector<float>& cum_u_perp,
    float L,
    float spd)
{
    const int N = (int)x.size();

    // --- mean direction at CURRENT time ---
    float ex, ey;
    compute_mean_direction(theta, ex, ey);

    // --- perpendicular direction ---
    float px = -ey;
    float py =  ex;

    // --- drift (mean velocity) ---
    float drift_x = spd * ex;
    float drift_y = spd * ey;

    double sum = 0.0;

    for (int i = 0; i < N; ++i)
    {
        float dx = min_image(x[i] - x_prev[i], L);
        float dy = min_image(y[i] - y_prev[i], L);

        // remove flock drift
        dx -= drift_x;
        dy -= drift_y;

        // project onto perpendicular direction
        float u_step = dx * px + dy * py;

        // accumulate
        cum_u_perp[i] += u_step;

        sum += (double)cum_u_perp[i] * (double)cum_u_perp[i];
    }

    return sum / N;
}

/* ============================================================
   MAIN
   ============================================================ */
int main()
{
    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 50;

    fs::create_directories("../../part_3d_properties_data/part_3g_uperp");

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

        std::map<int, std::vector<double>> accum;

        std::cout << "rho=" << rho
                  << " speed=" << spd << std::endl;

        #pragma omp parallel
        {
            std::vector<float> x, y, theta;
            std::vector<float> x_prev, y_prev, theta_prev;

            std::map<int, std::vector<double>> local_accum;

            #pragma omp for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim)
            {
                char fpath0[1024];

                std::snprintf(fpath0, sizeof(fpath0),
                    "../../../sim_data/vicsek_density_%.2f_speed_%.2f_noise_%.2f/"
                    "sim_%02d/sim_data/"
                    "sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                    rho, spd, noise,
                    sim,
                    noise, sim, 0);

                bool ok0 = read_snapshot(fpath0, N, x_prev, y_prev, theta_prev);
                if (!ok0) continue;

                std::vector<float> cum_u_perp(N, 0.0f);

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

                    bool ok = read_snapshot(fpath, N, x, y, theta);
                    if (!ok) continue;

                    if (k == 0)
                    {
                        local_accum[k].push_back(0.0);
                    }
                    else
                    {
                        double uperp2 = compute_uperp2_incremental(
                            x, y, theta,
                            x_prev, y_prev,
                            cum_u_perp,
                            L,
                            spd
                        );

                        local_accum[k].push_back(uperp2);
                    }

                    x_prev = x;
                    y_prev = y;
                    theta_prev = theta;
                }
            }

            #pragma omp critical
            {
                for (auto& kv : local_accum)
                {
                    auto& dst = accum[kv.first];
                    dst.insert(dst.end(), kv.second.begin(), kv.second.end());
                }
            }
        }

        if (accum.empty()) continue;

        char outname[256];

        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3g_uperp/"
            "uperp_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);

        out << "# k mean std\n";

        for (const auto& kv : accum)
        {
            double m, sd;
            mean_std(kv.second, m, sd);

            out << kv.first << " "
                << std::setprecision(10)
                << m << " "
                << sd << "\n";
        }

        std::cout << "wrote " << outname << std::endl;
    }

    std::cout << "DONE\n";
    return 0;
}