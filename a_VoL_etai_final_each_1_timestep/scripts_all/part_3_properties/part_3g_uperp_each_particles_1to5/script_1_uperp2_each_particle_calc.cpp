#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

static inline float min_image(float dx, float L)
{
    dx -= L * std::round(dx / L);
    return dx;
}

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

void compute_mean_direction(const std::vector<float>& theta_prev,
                            float& ex,
                            float& ey)
{
    float sx = 0.0f;
    float sy = 0.0f;

    for (float th : theta_prev) {
        sx += std::cos(th);
        sy += std::sin(th);
    }

    float norm = std::sqrt(sx * sx + sy * sy) + 1e-8f;
    ex = sx / norm;
    ey = sy / norm;
}

void accumulate_uperp_increment(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& x_prev,
    const std::vector<float>& y_prev,
    const std::vector<float>& theta_prev,
    std::vector<float>& cum_u_perp,
    float L,
    float spd,
    int dt_out)
{
    const int N = (int)x.size();

    float ex, ey;
    compute_mean_direction(theta_prev, ex, ey);

    float px = -ey;
    float py =  ex;

    float drift_x = spd * dt_out * ex;
    float drift_y = spd * dt_out * ey;

    for (int i = 0; i < N; ++i)
    {
        float dx = min_image(x[i] - x_prev[i], L);
        float dy = min_image(y[i] - y_prev[i], L);

        dx -= drift_x;
        dy -= drift_y;

        float u_step = dx * px + dy * py;
        cum_u_perp[i] += u_step;
    }
}

int main()
{
    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 5;

    fs::create_directories("../../part_3d_properties_data/part_3g_uperp_each_particles_1to5");

    std::vector<float> rho_vals   = {1.00f, 2.00f, 4.00f, 8.00f, 0.16f};
    std::vector<float> speed_vals = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};

    for (float rho : rho_vals)
    for (float spd : speed_vals)
    {
        float L = std::sqrt((float)N / rho);

        int T      = (int)std::round(L / spd);
        int T_max  = 5 * T;
        int dt_out = 1;

        int k_max = T_max / dt_out;

        std::cout << "rho=" << rho
                  << " speed=" << spd << std::endl;

        #pragma omp parallel
        {
            std::vector<float> x, y, theta;
            std::vector<float> x_prev, y_prev, theta_prev;

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
                if (!ok0) {
                    std::cerr << "ERROR: missing initial snapshot:\n"
                              << fpath0 << std::endl;
                    std::exit(EXIT_FAILURE);
                }

                std::vector<float> cum_u_perp(N, 0.0f);

                char outname[1024];
                std::snprintf(outname, sizeof(outname),
                    "../../part_3d_properties_data/part_3g_uperp_each_particles_1to5/"
                    "uperp2_particles_rho_%.2f_speed_%.2f_noise_%.2f_sim_%02d.txt",
                    rho, spd, noise, sim);

                std::ofstream out(outname);
                out << "# k uperp2_particle_001 ... uperp2_particle_" << std::setw(3)
                    << std::setfill('0') << N << "\n";
                out << std::setfill(' ');

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
                    if (!ok) {
                        std::cerr << "ERROR: missing snapshot file:\n"
                                  << fpath << std::endl;
                        std::exit(EXIT_FAILURE);
                    }

                    if (k > 0) {
                        accumulate_uperp_increment(
                            x, y,
                            x_prev, y_prev, theta_prev,
                            cum_u_perp,
                            L,
                            spd,
                            dt_out
                        );
                    }

                    out << k;
                    for (int i = 0; i < N; ++i)
                    {
                        double u2 = (double)cum_u_perp[i] * (double)cum_u_perp[i];
                        out << " " << std::setprecision(10) << std::fixed << u2;
                    }
                    out << "\n";

                    x_prev = x;
                    y_prev = y;
                    theta_prev = theta;
                }

                out.close();

                #pragma omp critical
                {
                    std::cout << "wrote " << outname << std::endl;
                }
            }
        }
    }

    std::cout << "DONE\n";
    return 0;
}
