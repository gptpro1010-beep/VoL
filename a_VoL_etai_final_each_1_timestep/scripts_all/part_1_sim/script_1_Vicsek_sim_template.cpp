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
    float L;
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
    float L = pm.L;
    float Lh = 0.5f * L;
    float r2 = pm.r * pm.r;

    int nc = floor(L / pm.r);
    float cs = L / nc;

    vector<int> head(nc*nc, -1), next(N, -1);

    for (int i = 0; i < N; i++) {
        int cx = min(int(x[i] / cs), nc-1);
        int cy = min(int(y[i] / cs), nc-1);
        int c  = cy*nc + cx;
        next[i] = head[c];
        head[c] = i;
    }

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

        int cx = xi / cs;
        int cy = yi / cs;

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx_c = (cx + dx + nc) % nc;
                int ny_c = (cy + dy + nc) % nc;
                int cell = ny_c*nc + nx_c;

                for (int j = head[cell]; j != -1; j = next[j]) {
                    float dxp = x[j] - xi;
                    float dyp = y[j] - yi;
                    if (dxp >  Lh) dxp -= L;
                    if (dxp < -Lh) dxp += L;
                    if (dyp >  Lh) dyp -= L;
                    if (dyp < -Lh) dyp += L;

                    if (dxp*dxp + dyp*dyp < r2) {
                        sc += c[j];
                        ss += s[j];
                    }
                }
            }
        }

        float th = atan2f(ss, sc)
                 + pm.noise * 2.0f * M_PI * (noise_seq[i] - 0.5f);

        nt[i] = th;
        nx[i] = fmodf(xi + pm.v*cosf(th) + L, L);
        ny[i] = fmodf(yi + pm.v*sinf(th) + L, L);
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
    pm.L = sqrt(N / rho);
    int T = (int)round(pm.L / pm.v);
    int T_max = 5 * T;
    int dt_out = 1;

    #pragma omp parallel for schedule(static)
    for (int sim = 1; sim <= n_sims; sim++) {
    
        mt19937 rng(12345u + sim);   // independent RNG per realization
        uniform_real_distribution<float> angle_dist(0.0f, 2.0f * M_PI);

        char sim_dir[64];
        sprintf(sim_dir, "sim_%02d", sim);

        fs::create_directories(string(sim_dir) + "/sim_data");

        char order_file[256];
        sprintf(order_file,
                "%s/order_param_noise_%.2f_sim_%02d_timeseries.txt",
                sim_dir, noise, sim);

        ofstream fout_order(order_file);

        vector<float> x(N), y(N), theta(N);

        float spacing = pm.L / n_side;
        int idx = 0;
        float theta0 = angle_dist(rng);
        for (int i = 0; i < n_side; i++) {
            for (int j = 0; j < n_side; j++) {
                x[idx] = i * spacing;
                y[idx] = j * spacing;
                theta[idx] = theta0;
                idx++;
            }
        }

        for (int t = 0; t <= T_max; t++) {

            if (t % T == 0)
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
