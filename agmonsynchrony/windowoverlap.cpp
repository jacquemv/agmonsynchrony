#include <math.h>

int time_windows_c(double* x, int n, double tau, int bilateral, double* interv)
{
    int i;
    if (n==0) return 0;
    double tau1, tau2;
    if (bilateral) {
        tau1 = tau2 = fabs(tau);
    } else {
        if (tau < 0) {
            tau1 = -tau;
            tau2 = 0;
        } else {
            tau1 = 0;
            tau2 = tau;
        }
    }
    interv[2*0] = x[0]-tau1;
    interv[2*0+1] = x[0]+tau2;
    int na = 1;
    for (i=1;i<n;i++) {
        if (interv[2*(na-1)+1]<x[i]-tau1) {
            interv[2*na] = x[i]-tau1;
            interv[2*na+1] = x[i]+tau2;
            na++;
        } else {
            interv[2*(na-1)+1] = x[i]+tau2;
        }
    }
    return na;
}

//-----------------------------------------------------------------------------
void coincidence_distribution_c(double* p, int n, double* P)
{
    if (n == 0) {
        P[0] = 1;
        return;
    }    
    P[0] = 1-p[0];
    P[1] = p[0];
    for (int m=1;m<n;m++) {
        for (int M=m+1;M>0;M--) {
            P[M] = p[m]*P[M-1]+(1-p[m])*P[M];
        }
        P[0] *= 1-p[m];
    }
}

//-----------------------------------------------------------------------------
void overlap_c(double* x, int na, double* y, int n, double tauj, int& Nc, 
               double* p)
{
    int i,j,j0 = 0;
    double c,d;
    Nc = 0;
    if (na==0)
        return;
    for (i=0;i<n;i++) {
        while (j0<na-1 && x[2*(j0+1)]<y[i])
            j0++;
        Nc += x[2*j0]<=y[i] && y[i]<=x[2*j0+1];
        p[i] = 0;
        for (j=j0;j<na;j++) {
            double aj = x[2*j], bj = x[2*j+1];
            if (aj>=y[i]+tauj) break;
            c = y[i]-tauj;
            if (aj>c) c = aj;
            d = y[i]+tauj;
            if (bj<d) d = bj;
            if (c<d) p[i] += d-c;
        }
        for (j=j0-1;j>=0;j--) {
            double aj = x[2*j], bj = x[2*j+1];
            if (bj<=y[i]-tauj) break;
            c = y[i]-tauj;
            if (aj>c) c = aj;
            d = y[i]+tauj;
            if (bj<d) d = bj;
            if (c<d) p[i] += d-c;
        }
        p[i] /= 2*tauj;
    }
}
