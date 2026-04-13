#ifndef M_PI
#define M_PI 3.1415926
#endif

#include <bits/stdc++.h>
#include <omp.h>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

vector<float> generateSequenceRandom(int count, mt19937& rng) {
    vector<float> seq(count);
    uniform_real_distribution<float> dist(0.0f, 1.0f);
    for (int i = 0; i < count; i++) seq[i] = dist(rng);
    return seq;
}

struct Params {
    float noise;
    float r;
    float v;
    int   N;
    float L_init;
};

float computeOrder(const vector<float>& theta, const Params& pm) {
    float sx = 0.0f, sy = 0.0f;
    for (int i = 0; i < pm.N; i++) {
        sx += cosf(theta[i]);
        sy += sinf(theta[i]);
    }
    return sqrtf(sx*sx + sy*sy) / pm.N;
}

void update(vector<float>& x, vector<float>& y, vector<float>& theta, const Params& pm, mt19937& rng) {
    int N = pm.N;
    float r2 = pm.r * pm.r;

    vector<float> nx(N), ny(N), nt(N);
    vector<float> c(N), s(N);
    for (int i = 0; i < N; i++) {
        c[i] = cosf(theta[i]);
        s[i] = sinf(theta[i]);
    }

    auto noise_seq = generateSequenceRandom(N, rng);

    for (int i = 0; i < N; i++) {
        float xi = x[i], yi = y[i];
        float sc = 0.0f, ss = 0.0f;

        for (int j = 0; j < N; j++) {
            float dxp = x[j] - xi;
            float dyp = y[j] - yi;

            if (dxp*dxp + dyp*dyp < r2) {
                sc += c[j];
                ss += s[j];
            }
        }

        float th = atan2f(ss, sc)
                 + pm.noise * 2.0f * M_PI * (noise_seq[i] - 0.5f);

        nt[i] = th;
        nx[i] = xi + pm.v*cosf(th);
        ny[i] = yi + pm.v*sinf(th);
    }

    x = move(nx);
    y = move(ny);
    theta = move(nt);
}

int main() {

    const float rho   = {RHO_VALUE};
    const float speed = {SPEED_VALUE};
    const float noise = {NOISE_VALUE};

    const int N = 100;
    const int n_side = 10;
    const int n_sims = 20;

    Params pm;
    pm.N = N;
    pm.v = speed;
    pm.r = 1.0f;
    pm.noise = noise;
    pm.L_init = sqrt(N / rho);

    int T = static_cast<int>(round(pm.L_init / pm.v));
    int T_max = 20000 * T;
    int save_every = 1 * T;
    int dt_out = save_every;

    #pragma omp parallel for schedule(static)
    for (int sim = 1; sim <= n_sims; sim++) {

        mt19937 rng(12345u + sim);

        char sim_dir[64];
        sprintf(sim_dir, "sim_%02d", sim);

        fs::create_directories(string(sim_dir) + "/sim_data");

        char order_file[256];
        sprintf(order_file,
                "%s/order_param_noise_%.2f_sim_%02d_timeseries.txt",
                sim_dir, noise, sim);

        ofstream fout_order(order_file);

        vector<float> x(N), y(N), theta(N);

        float spacing = pm.L_init / n_side;
        int idx = 0;
        float theta0 = 0.0f;
        for (int i = 0; i < n_side; i++) {
            for (int j = 0; j < n_side; j++) {
                x[idx] = i * spacing;
                y[idx] = j * spacing;
                theta[idx] = theta0;
                idx++;
            }
        }

        for (int t = 0; t <= T_max; t++) {

            fout_order << computeOrder(theta, pm) << "\n";

            if (t % dt_out == 0) {
                int k = t / dt_out;
                char snap[256];
                sprintf(snap,
                        "%s/sim_data/sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                        sim_dir, noise, sim, k);

                ofstream fout_snap(snap);
                for (int i = 0; i < N; i++)
                    fout_snap << x[i] << " " << y[i] << " " << theta[i] << "\n";
            }

            update(x, y, theta, pm, rng);
        }

        fout_order.close();
    }

    return 0;
}
