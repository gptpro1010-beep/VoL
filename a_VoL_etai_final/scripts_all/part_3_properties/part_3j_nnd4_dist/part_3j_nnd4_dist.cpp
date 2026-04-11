#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

static inline void min_image_unit(float& dx, float& dy)
{
    dx -= std::round(dx);
    dy -= std::round(dy);
}

static inline double sample_std_from_sums(double sum, double sumsq, int n)
{
    if (n <= 1) return 0.0;

    double var = (sumsq - (sum * sum) / (double)n) / (double)(n - 1);

    if (var < 0.0) var = 0.0;

    return std::sqrt(var);
}

int main()
{
    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;

    const int    K   = 141;
    const double r0  = 0.0;
    const double ell = 1.0 / std::sqrt((double)N);
    const double r1  = 2.0 * ell;
    const double dr  = (r1 - r0) / (double)(K - 1);

    std::vector<double> rgrid(K);

    for (int k = 0; k < K; k++)
        rgrid[k] = r0 + dr * (double)k;

    fs::create_directories("../../part_3d_properties_data/part_3j_nnd4_dist");

    std::vector<float> rhos   = {1, 2, 4, 8, 16};
    std::vector<float> speeds = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};

    for (float rho : rhos)
    for (float spd : speeds)
    {
        float L = std::sqrt((float)N / rho);

        int T      = (int)std::round(L / spd);
        int T_max  = 5000 * T;
        int dt_out = 5 * T;
        int k_max  = T_max / dt_out;

        char outname[256];

        std::snprintf(outname, sizeof(outname),
            "../../part_3d_properties_data/part_3j_nnd4_dist/"
            "nnd4_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho, spd, noise);

        std::ofstream out(outname);

        if (!out)
        {
            std::cerr << "ERROR writing " << outname << std::endl;
            return 1;
        }

        std::cout << "rho=" << std::fixed << std::setprecision(2) << rho
                  << " speed=" << std::fixed << std::setprecision(2) << spd
                  << std::endl;

        std::vector<double> sumF((k_max + 1) * K, 0.0);
        std::vector<double> sumF2((k_max + 1) * K, 0.0);
        std::vector<double> sumD(k_max + 1, 0.0);
        std::vector<double> sumD2(k_max + 1, 0.0);
        std::vector<int>    count(k_max + 1, 0);

        #pragma omp parallel
        {
            std::vector<float> x(N), y(N);
            std::vector<double> dnn4(N);
            std::vector<float> d2;
            std::vector<double> F(K);

            d2.reserve(N - 1);

            #pragma omp for schedule(dynamic)
            for (int sim = SIM_MIN; sim <= SIM_MAX; sim++)
            {
                for (int k = 0; k <= k_max; k++)
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

                    float xx, yy, th;
                    int cnt = 0;

                    while (cnt < N && (in >> xx >> yy >> th))
                    {
                        x[cnt] = xx / L;
                        y[cnt] = yy / L;
                        cnt++;
                    }

                    if (cnt != N) continue;

                    const int k_use = std::min(4, N - 1);

                    for (int i = 0; i < N; i++)
                    {
                        d2.clear();

                        for (int j = 0; j < N; j++)
                        {
                            if (i == j) continue;

                            float dx = x[j] - x[i];
                            float dy = y[j] - y[i];

                            min_image_unit(dx, dy);

                            d2.push_back(dx * dx + dy * dy);
                        }

                        std::partial_sort(d2.begin(), d2.begin() + k_use, d2.end());

                        double sum4 = 0.0;

                        for (int q = 0; q < k_use; q++)
                            sum4 += std::sqrt((double)d2[q]);

                        dnn4[i] = sum4 / (double)k_use;
                    }

                    std::sort(dnn4.begin(), dnn4.end());

                    for (int kk = 0; kk < K; kk++)
                    {
                        auto it = std::upper_bound(dnn4.begin(), dnn4.end(), rgrid[kk]);

                        F[kk] = (double)std::distance(dnn4.begin(), it) / (double)N;
                    }

                    double dbar = 0.0;

                    for (int i = 0; i < N; i++)
                        dbar += dnn4[i];

                    dbar /= (double)N;

                    #pragma omp critical
                    {
                        int base = k * K;

                        for (int kk = 0; kk < K; kk++)
                        {
                            sumF[base + kk]  += F[kk];
                            sumF2[base + kk] += F[kk] * F[kk];
                        }

                        sumD[k]  += dbar;
                        sumD2[k] += dbar * dbar;
                        count[k] += 1;
                    }
                }
            }
        }

        out << "# NND4 distribution (CDF): for each particle, average distance to its 4 nearest neighbors\n";
        out << "# rho=" << rho << " speed=" << spd << " noise=" << noise << "\n";
        out << "# snap_idx  r  mean_F  std_F  mean_d4  std_d4  nsamples\n";

        for (int k = 0; k <= k_max; k++)
        {
            int base = k * K;
            int n    = count[k];

            if (n == 0) continue;

            double meanD = sumD[k] / (double)n;
            double sdD   = sample_std_from_sums(sumD[k], sumD2[k], n);

            for (int kk = 0; kk < K; kk++)
            {
                double meanF = sumF[base + kk] / (double)n;
                double sdF   = sample_std_from_sums(sumF[base + kk], sumF2[base + kk], n);

                out << k << " "
                    << std::setprecision(10)
                    << rgrid[kk] << " "
                    << meanF << " "
                    << sdF << " "
                    << meanD << " "
                    << sdD << " "
                    << n << "\n";
            }
        }

        std::cout << " wrote " << outname << std::endl;
    }

    std::cout << "DONE.\n";

    return 0;
}
