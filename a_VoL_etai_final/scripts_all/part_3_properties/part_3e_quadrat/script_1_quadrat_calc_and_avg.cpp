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
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

static inline double sample_std_from_sums(double sum,double sumsq,int n)
{
    if(n<=1) return 0.0;

    double var=(sumsq-(sum*sum)/n)/(n-1);

    if(var<0.0) var=0.0;

    return std::sqrt(var);
}



static inline double gammln(double xx)
{
    return std::lgamma(xx);
}

static double gser(double a,double x)
{
    const int ITMAX=200;
    const double EPS=3e-14;

    if(x<=0.0) return 0.0;

    double sum=1.0/a;
    double del=sum;
    double ap=a;

    for(int n=1;n<=ITMAX;n++)
    {
        ap+=1.0;
        del*=x/ap;
        sum+=del;

        if(std::fabs(del)<std::fabs(sum)*EPS) break;
    }

    return sum*std::exp(-x+a*std::log(x)-gammln(a));
}

static double gcf(double a,double x)
{
    const int ITMAX=200;
    const double EPS=3e-14;
    const double FPMIN=1e-300;

    double b=x+1.0-a;
    double c=1.0/FPMIN;
    double d=1.0/b;
    double h=d;

    for(int i=1;i<=ITMAX;i++)
    {
        double an=-(double)i*((double)i-a);
        b+=2.0;

        d=an*d+b;
        if(std::fabs(d)<FPMIN) d=FPMIN;

        c=b+an/c;
        if(std::fabs(c)<FPMIN) c=FPMIN;

        d=1.0/d;

        double del=d*c;

        h*=del;

        if(std::fabs(del-1.0)<EPS) break;
    }

    return std::exp(-x+a*std::log(x)-gammln(a))*h;
}

static double gammp(double a,double x)
{
    if(x<0.0||a<=0.0) return 0.0;

    if(x<a+1.0) return gser(a,x);

    return 1.0-gcf(a,x);
}

static inline double chi2_cdf(double x,double df)
{
    if(x<=0.0) return 0.0;

    return gammp(0.5*df,0.5*x);
}

static inline double clamp_p(double p)
{
    const double pmin=1e-300;

    if(p<pmin) return pmin;
    if(p>1.0)  return 1.0;

    return p;
}



int main()
{

    const int   N=100;
    const float noise=0.01f;

    const int SIM_MIN=1;
    const int SIM_MAX=20;
    const int NSIMS=SIM_MAX-SIM_MIN+1;

    const int nq=8;
    const int m=nq*nq;
    const double dx=1.0/nq;
    const double df=(double)(m-1);

    fs::create_directories("../../part_3d_properties_data/part_3e_quadrat");


    std::vector<float> rho_vals   = {1.00f, 2.00f, 4.00f, 8.00f, 16.00f};
    std::vector<float> speed_vals = {0.01f, 0.02f, 0.04f, 0.08f, 0.16f};
    
    for(float rho : rho_vals)
    for(float spd : speed_vals)
    {

        float L=std::sqrt((float)N/rho);

        int T=(int)std::round(L/spd);

        int T_max=5000*T;

        int dt_out=5*T;

        int k_max=T_max/dt_out;


        std::vector<double> sumI(k_max+1,0.0);
        std::vector<double> sumI2(k_max+1,0.0);
        std::vector<double> sumLogP(k_max+1,0.0);
        std::vector<int> sumReject05(k_max+1,0);
        std::vector<int> sumReject10(k_max+1,0);


#pragma omp parallel
        {

            std::vector<float> x(N),y(N);
            std::vector<int> counts(m,0);

            std::vector<double> local_sumI(k_max+1,0.0);
            std::vector<double> local_sumI2(k_max+1,0.0);
            std::vector<double> local_sumLogP(k_max+1,0.0);
            std::vector<int> local_sumReject05(k_max+1,0);
            std::vector<int> local_sumReject10(k_max+1,0);


#pragma omp for schedule(dynamic)

            for(int sim=SIM_MIN;sim<=SIM_MAX;sim++)
            {

                for(int k=0;k<=k_max;k++)
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

                    while(cnt<N&&(in>>xx>>yy>>th))
                    {
                        double xn=(double)xx/L;
                        double yn=(double)yy/L;

                        xn-=std::floor(xn);
                        yn-=std::floor(yn);

                        x[cnt]=(float)xn;
                        y[cnt]=(float)yn;

                        cnt++;
                    }

                    if(cnt!=N) continue;


                    std::fill(counts.begin(),counts.end(),0);


                    for(int i=0;i<N;i++)
                    {
                        int ix=(int)std::floor((double)x[i]/dx);
                        int iy=(int)std::floor((double)y[i]/dx);

                        if(ix<0) ix=0;
                        if(iy<0) iy=0;

                        if(ix>=nq) ix=nq-1;
                        if(iy>=nq) iy=nq-1;

                        counts[ix*nq+iy]+=1;
                    }


                    double Nbar=(double)N/m;

                    double s2=0.0;

                    for(int q=0;q<m;q++)
                    {
                        double d=(double)counts[q]-Nbar;
                        s2+=d*d;
                    }

                    s2/=m;

                    double I=(Nbar>0.0)?(s2/Nbar):0.0;

                    double chi2=(m-1)*I;

                    double cdf=chi2_cdf(chi2,df);

                    double p=1.0-cdf;

                    p=clamp_p(p);


                    local_sumI[k]+=I;
                    local_sumI2[k]+=I*I;
                    local_sumLogP[k]+=std::log(p);

                    if(p<0.05) local_sumReject05[k]++;
                    if(p<0.10) local_sumReject10[k]++;

                }

            }


#pragma omp critical
            {

                for(int k=0;k<=k_max;k++)
                {
                    sumI[k]+=local_sumI[k];
                    sumI2[k]+=local_sumI2[k];
                    sumLogP[k]+=local_sumLogP[k];
                    sumReject05[k]+=local_sumReject05[k];
                    sumReject10[k]+=local_sumReject10[k];
                }

            }

        }


        char outname[256];

        std::snprintf(outname,sizeof(outname),
            "../../part_3d_properties_data/part_3e_quadrat/"
            "quadrat_avg_rho_%.2f_speed_%.2f_noise_%.2f.txt",
            rho,spd,noise);

        std::ofstream out(outname);


        out<<"# Quadrat CSR test\n";
        out<<"# rho="<<rho<<" speed="<<spd<<" noise="<<noise<<"\n";
        out<<"# columns: k  mean_I  std_I  p_fisher  frac_reject_05  frac_reject_10\n";


        for(int k=0;k<=k_max;k++)
        {
            double meanI=sumI[k]/NSIMS;

            double sdI=sample_std_from_sums(sumI[k],sumI2[k],NSIMS);

            double X=-2.0*sumLogP[k];

            double p_fisher=1.0-chi2_cdf(X,2.0*NSIMS);

            double frac05=(double)sumReject05[k]/NSIMS;

            double frac10=(double)sumReject10[k]/NSIMS;


            out<<k<<" "
               <<std::setprecision(10)
               <<meanI<<" "
               <<sdI<<" "
               <<p_fisher<<" "
               <<frac05<<" "
               <<frac10<<"\n";

        }

        std::cout<<" wrote "<<outname<<std::endl;

    }


    std::cout<<"DONE\n";

    return 0;
}