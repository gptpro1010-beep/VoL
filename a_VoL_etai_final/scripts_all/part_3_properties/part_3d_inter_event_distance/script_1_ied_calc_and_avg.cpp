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
#include <map>
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

    double var = (sumsq - (sum*sum)/n) / (n-1);

    if (var < 0.0) var = 0.0;

    return std::sqrt(var);
}



int main()
{
    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;
    const int NSIMS   = SIM_MAX - SIM_MIN + 1;

    const int   K   = 141;
    const double T0 = 0.0;
    const double T1 = std::sqrt(2.0)/2.0;
    const double DT = (T1 - T0)/(K-1);

    std::vector<double> tau(K);

    for (int k=0;k<K;k++)
        tau[k] = T0 + DT*k;

    fs::create_directories("../../part_3d_properties_data/part_3d_ied");


    std::vector<float> rho_vals   = {1.00f, 2.00f, 4.00f, 8.00f, 16.00f};
    std::vector<float> speed_vals = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};
    
    for (float rho : rho_vals)
    for (float spd : speed_vals)
    {

        float L = std::sqrt((float)N / rho);

        int T      = (int)std::round(L / spd);
        int T_max  = 5000*T;
        int dt_out = 5*T;

        int k_max = T_max / dt_out;

        std::map<int,std::vector<double>> sumH;
        std::map<int,std::vector<double>> sumH2;

        std::cout<<"rho="<<std::fixed<<std::setprecision(2)<<rho
                 <<" speed="<<spd<<std::endl;


#pragma omp parallel
        {

            std::vector<float> x(N),y(N);
            std::vector<float> dists;

            dists.reserve((size_t)N*(N-1)/2);

            std::map<int,std::vector<double>> local_sum;
            std::map<int,std::vector<double>> local_sum2;


#pragma omp for schedule(dynamic)

            for (int sim=SIM_MIN; sim<=SIM_MAX; sim++)
            {

                for (int k=0;k<=k_max;k++)
                {

                    char fpath[1024];

                    std::snprintf(fpath,sizeof(fpath),
                        "../../../sim_data/vicsek_density_%.2f_speed_%.2f_noise_%.2f/"
                        "sim_%02d/sim_data/"
                        "sim_data_noise_%.2f_sim_%02d_snap_%05d.txt",
                        rho,spd,noise,
                        sim,
                        noise,sim,k);

                    std::ifstream in(fpath);

                    if(!in) continue;


                    float xx,yy,th;

                    int cnt=0;

                    while(cnt<N && (in>>xx>>yy>>th))
                    {
                        x[cnt]=xx/L;
                        y[cnt]=yy/L;
                        cnt++;
                    }

                    if(cnt!=N) continue;


                    dists.clear();

                    for(int i=0;i<N;i++)
                    {
                        float xi=x[i];
                        float yi=y[i];

                        for(int j=i+1;j<N;j++)
                        {
                            float dx=x[j]-xi;
                            float dy=y[j]-yi;

                            min_image_unit(dx,dy);

                            dists.push_back(std::sqrt(dx*dx+dy*dy));
                        }
                    }

                    std::sort(dists.begin(),dists.end());

                    double M = (double)dists.size();

                    if(local_sum[k].empty())
                    {
                        local_sum[k].assign(K,0.0);
                        local_sum2[k].assign(K,0.0);
                    }

                    for(int kk=0;kk<K;kk++)
                    {
                        auto it = std::upper_bound(dists.begin(),dists.end(),(float)tau[kk]);

                        double c = std::distance(dists.begin(),it);

                        double H = c/M;

                        local_sum[k][kk]+=H;
                        local_sum2[k][kk]+=H*H;
                    }

                }
            }


#pragma omp critical
            {
                for(auto &kv : local_sum)
                {
                    int k = kv.first;

                    if(sumH[k].empty())
                    {
                        sumH[k].assign(K,0.0);
                        sumH2[k].assign(K,0.0);
                    }

                    for(int kk=0;kk<K;kk++)
                    {
                        sumH[k][kk]  += kv.second[kk];
                        sumH2[k][kk] += local_sum2[k][kk];
                    }
                }
            }

        }



        fs::path out_dir;

        char buf[256];

        std::snprintf(buf,sizeof(buf),
            "../../part_3d_properties_data/part_3d_ied/"
            "ied_avg_rho_%.2f_speed_%.2f_noise_%.2f",
            rho,spd,noise);

        out_dir = fs::path(buf);

        fs::create_directories(out_dir);



        for(auto &kv : sumH)
        {

            int k = kv.first;

            char fname[256];

            std::snprintf(fname,sizeof(fname),"k_%05d.txt",k);

            fs::path outpath = out_dir / fname;

            std::ofstream out(outpath);

            out<<"# IED EDF\n";
            out<<"# rho="<<rho
               <<" speed="<<spd
               <<" noise="<<noise
               <<" snap_idx="<<k<<"\n";
            out<<"# tau  mean_H  std_H\n";

            for(int kk=0;kk<K;kk++)
            {

                double mean = sumH[k][kk]/NSIMS;

                double sd = sample_std_from_sums(sumH[k][kk],sumH2[k][kk],NSIMS);

                out<<std::setprecision(10)
                   <<tau[kk]<<" "
                   <<mean<<" "
                   <<sd<<"\n";

            }

        }

        std::cout<<" wrote "<<out_dir.string()<<"/k_*.txt"<<std::endl;

    }


    std::cout<<"DONE\n";

    return 0;
}