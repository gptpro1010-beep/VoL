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

static inline void min_image(float& dx, float& dy, float L) {
    dx -= L * std::round(dx / L);
    dy -= L * std::round(dy / L);
}

static inline double sample_std_from_sums(double sum, double sumsq, int n) {
    if (n <= 1) return 0.0;
    double var = (sumsq - (sum * sum) / (double)n) / (double)(n - 1);
    if (var < 0.0) var = 0.0;
    return std::sqrt(var);
}

int main() {

    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;
    const int NSIMS   = SIM_MAX - SIM_MIN + 1;

    fs::create_directories("../../part_3d_properties_data/part_3c_gr");

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

        const float rmax = L / std::sqrt(2.0f);

        const int   n_side  = (int)std::sqrt((float)N);
        const float spacing = L / (float)n_side;
        const float dr      = spacing / 10.0f;

        const int nbins = (int)std::floor(rmax / dr);

        std::vector<float> r_in(nbins), r_out(nbins), r_mid(nbins), area(nbins);

        for (int b = 0; b < nbins; ++b) {
            r_in[b]  = b * dr;
            r_out[b] = (b + 1) * dr;
            r_mid[b] = 0.5f * (r_in[b] + r_out[b]);
            area[b]  = (float)M_PI * (r_out[b]*r_out[b] - r_in[b]*r_in[b]);
        }

        char buf[256];
        std::snprintf(buf,sizeof(buf),
            "../../part_3d_properties_data/part_3c_gr/"
            "gr_avg_rho_%.2f_speed_%.2f_noise_%.2f",
            rho,spd,noise);

        fs::path out_dir = fs::path(buf);
        fs::create_directories(out_dir);

        std::cout << "rho=" << rho << " speed=" << spd
                  << " k=[0," << k_max << "]\n";

        std::vector<double> sum ((k_max+1)*nbins,0.0);
        std::vector<double> sum2((k_max+1)*nbins,0.0);

        #pragma omp parallel
        {

            std::vector<float> xx(N),yy(N);
            std::vector<int> total_counts(nbins,0);

            #pragma omp for schedule(dynamic)

            for (int sim = SIM_MIN; sim <= SIM_MAX; ++sim)
            {

                for (int k = 0; k <= k_max; ++k)
                {

                    char fpath[1024];

                    std::snprintf(fpath,sizeof(fpath),
                        "../../../sim_data/vicsek_density_%.2f_speed_%.2f_noise_%.2f/"
                        "sim_%02d/sim_data/"
                        "sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                        rho,spd,noise,sim,noise,sim,k);

                    std::ifstream in(fpath);
                    if (!in) continue;

                    float x,y,th;
                    int cnt=0;

                    while (cnt<N && (in>>x>>y>>th)) {
                        xx[cnt]=x;
                        yy[cnt]=y;
                        cnt++;
                    }

                    if (cnt!=N) continue;

                    std::fill(total_counts.begin(),total_counts.end(),0);

                    for (int i=0;i<N;i++) {

                        float xi=xx[i], yi=yy[i];

                        for (int j=0;j<N;j++) {

                            if (j==i) continue;

                            float dx=xx[j]-xi;
                            float dy=yy[j]-yi;

                            min_image(dx,dy,L);

                            float r = std::sqrt(dx*dx+dy*dy);

                            if (r>=rmax) continue;

                            int b = (int)(r/dr);

                            if (b>=0 && b<nbins)
                                total_counts[b]++;
                        }
                    }

                    int base = k*nbins;

                    #pragma omp critical
                    {

                        for (int b=0;b<nbins;b++) {

                            double dens = (double)total_counts[b] /
                                          ((double)N * (double)area[b]);

                            sum [base+b] += dens;
                            sum2[base+b] += dens*dens;
                        }
                    }

                }

            }

        }

        for (int k=0;k<=k_max;k++)
        {

            char fname[64];
            std::snprintf(fname,sizeof(fname),"k_%05d.txt",k);

            fs::path outpath = out_dir / fname;

            std::ofstream out(outpath);

            out << "# radial density vs r\n";
            out << "# rho=" << rho
                << " speed=" << spd
                << " noise=" << noise
                << " snap_idx=" << k << "\n";

            out << "# r_mid r_in r_out mean_dens std_dens nSims\n";

            int base = k*nbins;

            for (int b=0;b<nbins;b++) {

                double mean = sum[base+b]/(double)NSIMS;

                double sd = sample_std_from_sums(
                    sum[base+b],sum2[base+b],NSIMS);

                out << std::setprecision(10)
                    << r_mid[b] << " "
                    << r_in[b]  << " "
                    << r_out[b] << " "
                    << mean << " "
                    << sd   << " "
                    << NSIMS << "\n";
            }

        }

        std::cout << "  wrote " << out_dir.string() << "/k_*.txt\n";

    }

    std::cout << "DONE.\n";
    return 0;
}