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

static inline void min_image_unit(float& dx, float& dy){
    dx -= std::round(dx);
    dy -= std::round(dy);
}

static inline double sample_std_from_sums(double sum,double sumsq,int n){
    if(n<=1) return 0.0;
    
    double var=(sumsq-(sum*sum)/(double)n)/(double)(n-1);

    if(var<0.0) var=0.0;

    return std::sqrt(var);
}



int main(){

    const int   N     = 100;
    const float noise = 0.01f;

    const int SIM_MIN = 1;
    const int SIM_MAX = 20;
    const int NSIMS   = SIM_MAX-SIM_MIN+1;

    const int    K = 141;

    const double r0 = 0.0;
    const double ell = 1.0/std::sqrt((double)N);
    const double r1 = 2.0*ell;
    const double dr = (r1-r0)/(double)(K-1);

    std::vector<double> rgrid(K);

    for(int k=0;k<K;k++)
        rgrid[k]=r0+dr*(double)k;

    fs::create_directories("../../part_3d_properties_data/part_3f_nnd_dist");

    std::vector<float> rhos   = {1, 2, 4, 8, 16};
    std::vector<float> speeds = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};
    
    for(float rho : rhos)
    for(float spd : speeds)
    {
        float L=std::sqrt((float)N/rho);

        int T=(int)std::round(L/spd);

        int T_max=5000*T;

        int dt_out=5*T;

        int k_max=T_max/dt_out;
        
        char outname[256];

        std::snprintf(outname,sizeof(outname),
            "../../part_3d_properties_data/part_3f_nnd_dist/"
            "nnd_dist_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho,spd,noise);
        
        std::ofstream out(outname);
        
        if(!out){
            std::cerr<<"ERROR writing "<<outname<<std::endl;
            return 1;
        }

        std::cout<<"rho="<<std::fixed<<std::setprecision(2)<<rho
                 <<" speed="<<std::fixed<<std::setprecision(2)<<spd
                 <<" sims="<<NSIMS<<std::endl;

        std::vector<double> sumF((k_max+1)*K,0.0);
        std::vector<double> sumF2((k_max+1)*K,0.0);
        std::vector<double> sumD(k_max+1,0.0);
        std::vector<double> sumD2(k_max+1,0.0);
    
#pragma omp parallel
{
            std::vector<float> x(N),y(N);
            std::vector<double> dnn(N);
            std::vector<double> F(K);


#pragma omp for schedule(dynamic)
            for(int sim=SIM_MIN;sim<=SIM_MAX;sim++){
                for(int k=0;k<=k_max;k++){
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

                    while(cnt<N&&(in>>xx>>yy>>th)){
                        x[cnt]=xx/L;
                        y[cnt]=yy/L;
                        cnt++;
                    }

                    if(cnt!=N) continue;

                    for(int i=0;i<N;i++){

                        double dmin=1e9;

                        for(int j=0;j<N;j++){

                            if(i==j) continue;

                            float dx=x[j]-x[i];
                            float dy=y[j]-y[i];

                            min_image_unit(dx,dy);

                            double d=std::sqrt(dx*dx+dy*dy);

                            if(d<dmin) dmin=d;
                        }
                        dnn[i]=dmin;
                    }

                    std::sort(dnn.begin(),dnn.end());

                    for(int kk=0;kk<K;kk++){
                        auto it=std::upper_bound(dnn.begin(),dnn.end(),rgrid[kk]);

                        F[kk]=(double)std::distance(dnn.begin(),it)/(double)N;
                    }



                    double dbar=0.0;

                    for(int i=0;i<N;i++)
                        dbar+=dnn[i];

                    dbar/=(double)N;



#pragma omp critical
{
                        int base=k*K;

                        for(int kk=0;kk<K;kk++){
                            sumF[base+kk]+=F[kk];
                            sumF2[base+kk]+=F[kk]*F[kk];
                        }

                        sumD[k]+=dbar;
                        sumD2[k]+=dbar*dbar;
                    }

                }

            }

        }


        out << "# NND distribution (CDF)\n";
        out << "# rho=" << rho << " speed=" << spd << " noise=" << noise << "\n";
        out << "# snap_idx  r  mean_F  std_F  mean_d\n";
        
        for(int k=0;k<=k_max;k++){
            int base = k*K;
        
            double meanD = sumD[k]/NSIMS;
            double sdD   = sample_std_from_sums(sumD[k],sumD2[k],NSIMS);
        
            for(int kk=0;kk<K;kk++){
                double meanF = sumF[base+kk]/NSIMS;
        
                double sdF = sample_std_from_sums(
                                sumF[base+kk],
                                sumF2[base+kk],
                                NSIMS);
        
                out << k << " "
                    << std::setprecision(10)
                    << rgrid[kk] << " "
                    << meanF << " "
                    << sdF << " "
                    << meanD << "\n";
            }
        }
        std::cout<<" wrote "<<outname<<std::endl;
    }

    std::cout<<"DONE\n";

    return 0;

}