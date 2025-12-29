/********************************************************************************************************
--> processObservables.c

    : Processes (e.g. averaging and filtering) radar observables

    # Functions
    1. DP averaging, filtering, and unfolding
    2. Kdp derivation
    3. lookup table for azimuth angles in different elevation angles
    4. Alpha estimation for specific attenuation (default: quadratic function)

***********************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

#define ELEM_SWAP(a,b) {float t = (a); (a) = (b); (b) = t;}

void getvaziLookup(struct volume_info header,
                   int nangle,
                   float **azimuth,
                   short **lookup)
{
    int ang, nazi, bazi, azi;
    float base, dmin, adiff;

    for (ang=1; ang<nangle; ang++) {
         //if (ang<=nsuper-1)
         if (header.typeangle[ang]==SUPER)
             nazi=MAXIMUM_AZIMUTH_SUPER;
         else
             nazi=MAXIMUM_AZIMUTH;

         for (bazi=0; bazi<MAXIMUM_AZIMUTH_SUPER; bazi++) {
              base=azimuth[bazi][0];
              dmin=1000.;
              for (azi=0; azi<nazi; azi++) {
                   adiff=fabs(base-azimuth[azi][ang]);
                   if (adiff<=dmin) {
                       dmin=adiff;
                       lookup[bazi][ang]=azi;
                       if (dmin<0.1)
                           break;
                   }
              }
         }
    }
}

void unfoldPhidp(struct volume_info header,
                 int ***Phidp_3d,
                 int ***Rho_3d,
                 int angle,
                 int scale)
{
    float *Phidp_ray, *Rho_ray, Rho_thresh=0.85;
    int ang, azi, rng, nazi, flag;

    /* arrays */
    Phidp_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    Rho_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         flag=0;
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (Phidp_3d[azi][rng][ang]!=NODATA)
                  Phidp_ray[rng]=(float)Phidp_3d[azi][rng][ang]/scale;
              else
                  Phidp_ray[rng]=NODATA;

              if (Rho_3d[azi][rng][ang]!=NODATA)
                  Rho_ray[rng]=(float)Rho_3d[azi][rng][ang]/scale;
              else
                  Rho_ray[rng]=NODATA;
         }
         flag=unfoldRay(Phidp_ray,Rho_ray,Rho_thresh);
         if (flag>0) {
             for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                  if (Phidp_3d[azi][rng][ang]!=NODATA)
                      Phidp_3d[azi][rng][ang]=(int)(Phidp_ray[rng]*scale+0.5);
             }
         }
    }
    free(Phidp_ray);
    free(Rho_ray);
}

int unfoldRay(float *Phidp,
              float *Rho,
              float Rho_thresh)
{

    int i, j, count, v_size;
    int max_hist_median_size, hist_median_thresh, min_valid_data, min_bin_unfold_start, valid_data,
        unfold_flag, unfold, unfold_yes=1, unfold_no=0;
    float *unfolded_Phidp, *historical_median_vector;
    float stddev;
    float fold_in_degrees, max_stddev, half_fold_in_degrees, historical_median,
          phi_diff, phi_single_fold, phi_double_fold, single_fold_diff, double_fold_diff;

    /* init constants */
    unfold=unfold_no;
    fold_in_degrees=360.;
    max_hist_median_size=30;
    hist_median_thresh=25;
    max_stddev=fold_in_degrees/3.0;
    min_valid_data=15;
    half_fold_in_degrees=fold_in_degrees/2.0;
    min_bin_unfold_start=240; /* 60 km assuming 0.25 km bin size */

    /* array */
    unfolded_Phidp=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    historical_median_vector=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    /* init variables */
    valid_data=0;
    count=0;
    historical_median=0;

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         unfold_flag=0;
         if (Phidp[i]==NODATA) {
             unfolded_Phidp[i]=Phidp[i];
             continue;
         }

         if (Rho[i]!=NODATA && Rho[i]>=Rho_thresh)
             valid_data++;

         if (i>=max_hist_median_size) {
             count=0;
             v_size=0;
             for (j=1; j<=max_hist_median_size; j++) {
                  if (Rho[i-j]>=Rho_thresh && unfolded_Phidp[i-j]!=NODATA && unfolded_Phidp[i-j]!=INITIAL_VALUE) {
                      historical_median_vector[v_size]=unfolded_Phidp[i-j];
                      v_size++;
                      count++;
                  }
             }

             if (count>hist_median_thresh) {
                 stddev=Standard_deviation(v_size,historical_median_vector);
                 if (stddev<max_stddev)
                     historical_median=medianfilter(historical_median_vector,v_size); // median of the previous 30 gates
             }
         }

         phi_diff=Phidp[i]-historical_median;
         if (phi_diff<0)
             phi_diff=-phi_diff;

         if (phi_diff>=half_fold_in_degrees && valid_data>min_valid_data && i>min_bin_unfold_start) {
             phi_single_fold=Phidp[i]+fold_in_degrees;
             phi_double_fold=Phidp[i]+2*fold_in_degrees;

             single_fold_diff=historical_median-phi_single_fold;
             if (single_fold_diff<0)
                 single_fold_diff=-single_fold_diff;

             double_fold_diff=historical_median-phi_double_fold;
             if (double_fold_diff<0)
                 double_fold_diff=-double_fold_diff;

             if (phi_diff>single_fold_diff)
                 unfold_flag=1;

             if (single_fold_diff>double_fold_diff)
                 unfold_flag=2;
	     }

	     if (unfold_flag==1) {
	         unfolded_Phidp[i]=phi_single_fold;
	         Rho[i]=Rho_thresh;
	         unfold=unfold_yes;
	     }
	     else if (unfold_flag==2) {
	         unfolded_Phidp[i]=phi_double_fold;
	         Rho[i]=Rho_thresh;
	         unfold=unfold_yes;
	     }
	     else
	         unfolded_Phidp[i]=Phidp[i];
    }

    /* copy array */
    if (unfold==unfold_yes) {
        for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
             Phidp[i]=unfolded_Phidp[i];
        }
    }

    return unfold;
}

float Standard_deviation(int num,
                         float *data)
{
    float sum, sq, d;
    int i;

    if (num<=1)
	    return 0;

    sum=0;
    sq=0;

    for (i=0; i<num; i++) {
	     d=data[i];
	     sum += d;
	     sq +=d*d;
    }
    d=sum/num;
    d=(sq-d*d*num)/(num - 1);		/* Non-biased STD */

    if (d<0)
	    d=0;

    return (sqrt(d));
}

void averageZ(struct volume_info header,
              short ***z,
              int w,
              int angle)

{
    int ang, azi, rng, nazi;
    float *z_ray;

    /* arrays */
    z_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE)
                  z_ray[rng]=int_dBZ2float_dBZ(z[azi][rng][ang]);
              else if (z[azi][rng][ang]==INITIAL_VALUE)
                  z_ray[rng]=INITIAL_VALUE;
              else if (z[azi][rng][ang]==NODATA)
                  z_ray[rng]=NODATA;
         }
         averageZ_ray(z_ray,w);

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z_ray[rng]!=NODATA && z_ray[rng]!=INITIAL_VALUE)
                  z[azi][rng][ang]=float_dBZ2int_dBZ(z_ray[rng]);
              else if (z_ray[rng]==INITIAL_VALUE)
                  z[azi][rng][ang]=INITIAL_VALUE;
              else if (z_ray[rng]==NODATA)
                  z[azi][rng][ang]=NODATA;
         }
    }
    free(z_ray);
}

void averageZ_ray(float *ray,
                  int w)
{
    int i, hw, cnt, j, k, p1, p2, end, cnt_nodata, cnt_initial, range_start=8;
    float one_over[30], *smd, sum;

    /* array */
    smd=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    for (i=1; i<=w; i++) {
         one_over[i]=1/(float)i;
    }
    hw=w/2;
    p1=-hw-1;
    p2=hw;
    end=MAXIMUM_RANGE_SUPER-1;
    cnt=0;
    sum=0;

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if ((i%128)==0) {	/* reinitialize sum and cnt for precision */
             sum=0;
             cnt=0;
             cnt_nodata=0;
             cnt_initial=0;
             for (j=-hw-1; j<hw; j++) {
                  k=i+j;
                  if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]!=NODATA && ray[k]!=INITIAL_VALUE) {
                      sum+=pow(10.,ray[k]/10.);
                      cnt++;
                  }
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]==NODATA)
                      cnt_nodata++;
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]==INITIAL_VALUE)
                      cnt_initial++;
             }
         }

         if (p1>=0 && ray[p1]!=NODATA && ray[p1]!=INITIAL_VALUE) {
             sum-=pow(10.,ray[p1]/10.);
             cnt--;
         }
         else if (p1>0 && ray[p1]==NODATA)
             cnt_nodata--;
         else if (p1>0 && ray[p1]==INITIAL_VALUE)
             cnt_initial--;

         if (p2<=end && ray[p2]!=NODATA && ray[p2]!=INITIAL_VALUE) {
             sum+=pow(10.,ray[p2]/10.);
             cnt++;
         }
         else if (p2<=end && ray[p2]==NODATA)
             cnt_nodata++;
         else if (p2<=end && ray[p2]==INITIAL_VALUE)
             cnt_initial++;

         if (cnt>0)
             smd[i]=sum*one_over[cnt];
         else {
             if (cnt_nodata>cnt_initial)
                 smd[i]=NODATA;
             else
                 smd[i]=INITIAL_VALUE;
         }
         p1++;
         p2++;
    }

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (i<range_start)
             ray[i]=NODATA;
         else
             ray[i]=10*log10(smd[i]);
    }
    free(smd);
}

void averageDPvariables(struct volume_info header,
                        int ***DP,
                        int w,
                        int angle,
                        int scale)

{
    int ang, azi, rng, nazi;
    float *DP_ray;

    /* arrays */
    DP_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (DP[azi][rng][ang]!=NODATA && DP[azi][rng][ang]!=INITIAL_VALUE)
                  DP_ray[rng]=(float)DP[azi][rng][ang]/scale;
              else if (DP[azi][rng][ang]==INITIAL_VALUE)
                  DP_ray[rng]=INITIAL_VALUE;
              else if (DP[azi][rng][ang]==NODATA)
                  DP_ray[rng]=NODATA;
         }
         average_ray(DP_ray,w);

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (DP_ray[rng]!=NODATA)
                  DP[azi][rng][ang]=(int)(DP_ray[rng]*scale+0.5);
              else
                  DP[azi][rng][ang]=NODATA;
         }
    }
    free(DP_ray);
}

void average_ray(float *ray,
                 int w)
{
    int i, hw, cnt, j, k, p1, p2, end, cnt_nodata, cnt_initial, range_start=8;
    float one_over[30], *smd, sum;

    /* array */
    smd=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    for (i=1; i<=w; i++) {
         one_over[i]=1/(float)i;
    }
    hw=w/2;
    p1=-hw-1;
    p2=hw;
    end=MAXIMUM_RANGE_SUPER-1;
    cnt=0;
    sum=0;

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if ((i%128)==0) {	/* reinitialize sum and cnt for precision */
             sum=0;
             cnt=0;
             cnt_nodata=0;
             cnt_initial=0;
             for (j=-hw-1; j<hw; j++) {
                  k=i+j;
                  if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]!=NODATA &&ray[k]!=INITIAL_VALUE) {
                      sum+=ray[k];
                      cnt++;
                  }
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]==NODATA)
                      cnt_nodata++;
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && ray[k]==INITIAL_VALUE)
                      cnt_initial++;
             }
         }

         if (p1>=0 && ray[p1]!=NODATA && ray[p1]!=INITIAL_VALUE) {
             sum-=ray[p1];
             cnt--;
         }
         else if (p1>0 && ray[p1]==NODATA)
             cnt_nodata--;
         else if (p1>0 && ray[p1]==INITIAL_VALUE)
             cnt_initial--;

         if (p2<=end && ray[p2]!=NODATA && ray[p2]!=INITIAL_VALUE) {
             sum+=ray[p2];
             cnt++;
         }
         else if (p2<=end && ray[p2]==NODATA)
             cnt_nodata++;
         else if (p2<=end && ray[p2]==INITIAL_VALUE)
             cnt_initial++;

         if (cnt>0)
             smd[i]=sum*one_over[cnt];
         else {
             if (cnt_nodata>cnt_initial)
                 smd[i]=NODATA;
             else
                 smd[i]=INITIAL_VALUE;
         }
         p1++;
         p2++;
    }

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (i<range_start)
             ray[i]=NODATA;
         else
             ray[i]=smd[i];
    }
    free(smd);
}

void medianDPvariables(struct volume_info header,
                       int ***DP,
                       int w,
                       int angle,
                       int scale)
{
    int ang, azi, rng, nazi;
    float *DP_ray;

    /* arrays */
    DP_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (DP[azi][rng][ang]!=NODATA && DP[azi][rng][ang]!=INITIAL_VALUE)
                  DP_ray[rng]=(float)DP[azi][rng][ang]/scale;
              else if (DP[azi][rng][ang]==INITIAL_VALUE)
                  DP_ray[rng]=INITIAL_VALUE;
              else if (DP[azi][rng][ang]==NODATA)
                  DP_ray[rng]=NODATA;
         }
         medianfilter_ray(DP_ray,w);

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (DP_ray[rng]!=NODATA)
                  DP[azi][rng][ang]=(int)(DP_ray[rng]*scale+0.5);
              else
                  DP[azi][rng][ang]=NODATA;
         }
    }
    free(DP_ray);
}

void medianfilter_ray(float *ray,
                      int w)
{
    int i, j, k, hw, b1, b2, cnt, cnt_nodata, cnt_initial, range_start=8;
    float *med, buf[30], v;

    /* array */
    med=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    hw=w/2;
    b1=hw;
    b2=MAXIMUM_RANGE_SUPER-hw;

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         cnt=0;
         cnt_nodata=0;
         cnt_initial=0;
         if (i>b1 && i<b2) {		/* efficient in the middle of array */
             for (j=-hw; j<=hw; j++) {
                  v=ray[i+j];
                  if (v!=NODATA && v!=INITIAL_VALUE) {
                      buf[cnt]=v;
                      cnt++;
                  }
                  else if (v==NODATA)
                      cnt_nodata++;
                  else if (v==INITIAL_VALUE)
                      cnt_initial++;
             }
         }
         else {				/* near the boundaries of array */
             for (j=-hw; j<=hw; j++) {
                  k=i+j;
                  if (k>=0 && k<MAXIMUM_RANGE_SUPER) {
                      if (ray[k]!=NODATA && ray[k]!=INITIAL_VALUE) {
                          buf[cnt]=ray[k];
                          cnt++;
                      }
                      else if (ray[k]==NODATA)
                          cnt_nodata++;
                      else if (ray[k]==INITIAL_VALUE)
                          cnt_initial++;
                  }
             }
         }

         if (cnt>0)
             med[i]=medianfilter(buf,cnt);
         else {
             if (cnt_nodata>cnt_initial)
                 med[i]=NODATA;
             else
                 med[i]=INITIAL_VALUE;
         }
    }

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (i<range_start)
             ray[i]=NODATA;
         else
             ray[i]=med[i];
    }
    free(med);
}

float medianfilter(float *arr,
                   int n)
{

/**************************************************************************

    Performs median filtering on array "arr" of size "n". The medial
    value is returned. If n <= 0, 0 is returned. This routine is adapted
    from scit_filter_quickselect.c (Yukuan Song) and based on the
    algorithm described in "numerical recipies in C", Second edition,
    Cambridge University Press, 1992, section 8.5, ISBN 0-521-43108-5.
    ((low + high + 1) / 2) selects the upper element in case of n is
    even while ((low + high) / 2) selects the lower.

***************************************************************************/

    int low, high, median;
    int middle, ll, hh;

    if (n<=0)
        return (0);

    low=0;
    high=n-1;
    median=(low+high+1)/2;

    for (;;) {
        if (high <= low)		/* One element only */
            return arr[median] ;

        if (high==low+1) {		/* Two elements only */
            if (arr[low]>arr[high])
                ELEM_SWAP(arr[low],arr[high]);
            return arr[median];
        }

	    /* Find median of low, middle and high items; swap into position low */
	    middle=(low+high+1)/2;
	    if (arr[middle]>arr[high])
	        ELEM_SWAP(arr[middle],arr[high]);
	    if (arr[low]>arr[high])
	        ELEM_SWAP(arr[low],arr[high]);
	    if (arr[middle]>arr[low])
	        ELEM_SWAP(arr[middle],arr[low]);

	    /* Swap low item (now in position middle) into position (low + 1) */
	    ELEM_SWAP(arr[middle],arr[low+1]);

	    /* Nibble from each end towards middle, swapping items when stuck */
	    ll=low + 1;
	    hh=high;
	    for (;;) {
            do ll++;
            while (arr[low]>arr[ll]);
            do hh--;
            while (arr[hh]>arr[low]);

            if (hh<ll)
            break;

            ELEM_SWAP(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low],arr[hh]);

        /* Re-set active partition */
        if (hh<= median)
            low=ll;
        if (hh>=median)
            high=hh - 1;
    }
}

void computeKdp(struct volume_info header,
                int ***phi,
                int ***rho,
                short ***z,
                int ***kdp,
                int w,
                int angle,
                int scale)
{
    int ang, azi, rng, nazi;
    float *phi_ray, *rho_ray, *kdp_ray, *z_ray;

    /* arrays */
    phi_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    rho_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    kdp_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    z_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (rho[azi][rng][ang]!=NODATA && rho[azi][rng][ang]!=INITIAL_VALUE)
                  rho_ray[rng]=(float)rho[azi][rng][ang]/scale;
              else if (rho[azi][rng][ang]==NODATA)
                  rho_ray[rng]=NODATA;
              else if (rho[azi][rng][ang]==INITIAL_VALUE)
                  rho_ray[rng]=INITIAL_VALUE;


              if (phi[azi][rng][ang]!=NODATA && phi[azi][rng][ang]!=INITIAL_VALUE)
                  phi_ray[rng]=(float)phi[azi][rng][ang]/scale;
              else if (phi[azi][rng][ang]==NODATA)
                  phi_ray[rng]=NODATA;
              else if (phi[azi][rng][ang]==INITIAL_VALUE)
                  phi_ray[rng]=INITIAL_VALUE;

              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE)
                  z_ray[rng]=int_dBZ2float_dBZ(z[azi][rng][ang]);
              else if (z[azi][rng][ang]==INITIAL_VALUE)
                  z_ray[rng]=INITIAL_VALUE;
              else if (z[azi][rng][ang]==NODATA)
                  z_ray[rng]=NODATA;
         }
         getkdp_ray(phi_ray,rho_ray,z_ray,w,kdp_ray);

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (kdp_ray[rng]!=NODATA)
                  kdp[azi][rng][ang]=(int)(kdp_ray[rng]*scale+0.5);
              else
                  kdp[azi][rng][ang]=NODATA;
         }
    }
    free(phi_ray);
    free(rho_ray);
    free(kdp_ray);
    free(z_ray);
}

void adjustKdp(struct volume_info header,
               int angle,
               short ***z,
               int ***shortkdp,
               int ***longkdp)
{
    int ang, azi, rng, nazi;
    float zval, zshort=40.;

    ang=angle-1;

    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  if (zval>=zshort && shortkdp[azi][rng][ang]>0)
                      longkdp[azi][rng][ang]=shortkdp[azi][rng][ang];
              }
         }
    }
}

void getkdp_ray(float *phi,
                float *rho,
                float *z,
                int w,
                float *kdp)
{
    int hw, b1, b2, i, j, k, cnt, cnt_nodata, cnt_initial;
    float buf[30], z_thresh, rho_thresh;

    z_thresh=40.;
    rho_thresh=0.85;

    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (rho[i]!=NODATA && rho[i]<rho_thresh) {
             kdp[i]=INITIAL_VALUE;
             continue;
         }
         else if (rho[i]==NODATA) {
             kdp[i]=NODATA;
             continue;
         }

         hw=w/2;
         b1=hw;
         b2=MAXIMUM_RANGE_SUPER-hw;

         cnt=0;
         cnt_nodata=0;
         cnt_initial=0;
         if (i>b1 && i<b2) {		/* efficient in the middle of array */
             for (j=-hw; j<=hw; j++) {
                  k=i+j;
                  if (phi[k]!=NODATA && phi[k]!=INITIAL_VALUE) {
                      buf[cnt]=phi[k];
                      cnt++;
                  }
                  else if (phi[k]==NODATA)
                      cnt_nodata++;
                  else if (phi[k]==INITIAL_VALUE)
                      cnt_initial++;
             }
             if (cnt>0) {
                 if ((cnt%2)==0)
                     kdp[i]=Calculate_lls_kdp(buf,cnt,BIN_SUPER/KM);
                 else 		/* cnt is an odd number */
                     kdp[i]=Calculate_kdp(buf,cnt,BIN_SUPER/KM);
             }
             else {
                 if (cnt_nodata>cnt_initial)
                     kdp[i]=NODATA;
                 else
                     kdp[i]=INITIAL_VALUE;
             }
	     }
	     else {				/* near the boundaries of array */
             for (j=-hw; j<=hw; j++) {
                  k=i+j;
                  if (k>=0 && k<MAXIMUM_RANGE_SUPER && phi[k]!=NODATA && phi[k]!=INITIAL_VALUE) {
                      buf[cnt]=phi[k];
		              cnt++;
		          }
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && phi[k]==NODATA)
                      cnt_nodata++;
                  else if (k>=0 && k<MAXIMUM_RANGE_SUPER && phi[k]==INITIAL_VALUE)
                      cnt_initial++;
		     }
		     if (cnt>0) {
                 if ((cnt%2)==0)
                     kdp[i]=Calculate_lls_kdp(buf,cnt,BIN_SUPER/KM);
                 else 		/* cnt is an odd number */
                     kdp[i]=Calculate_kdp(buf,cnt,BIN_SUPER/KM);
             }
             else {
                 if (cnt_nodata>cnt_initial)
                     kdp[i]=NODATA;
                 else
                     kdp[i]=INITIAL_VALUE;
             }
         }
    }
}

float Calculate_kdp(float *phi,
                    int m,
                    float g_size)
{
    static int prev_m=0;
    static float factor=0, st_r=0;
    float j, sum, kdp;
    int i;

    /* The following is correct only if m is odd and >= 3 */
    if (m!=prev_m) {
        factor=(float)6/(g_size*(float)(m*(m-1)*(m+1)));
        st_r=-((m-1)/2);
        prev_m=m;
    }

    sum=0;
    j=st_r;
    for (i=0; i<m; i++) {
         sum+=j*phi[i];
         j+=1.;
    }
    kdp=factor*sum;

    return(kdp);
}

float Calculate_lls_kdp(float *phi,
                        int m,
                        float g_size)
{
    int j;
    float sx, sy, sxx, sxy, x, y, hw, nd;

    if (m<=1)
	    return INITIAL_VALUE;

    sx=sy=sxx=sxy=0;
    hw=(m-1)*(float).5;

    for (j=0; j<m; j++) {
         x=((float)j-hw)*g_size;
         y=phi[j];
         sx+=x;
         sy+=y;
         sxy+=y*x;
         sxx+=x*x;
    }
    nd=(float)m;
    return (float).5*(nd*sxy-sx*sy)/(nd*sxx-sx*sx);
}

float computeAlpha(struct volume_info header,
                   short ***z,
                   int ***zdr,
                   int ***rho,
                   short ***melting,
                   int scale)
{
    FILE *ftemp;
    struct tm t;
    int azi, rng, i, j, k, npair, cpair, ni, nis, nsub,
        sample_rng_min=20*4, sample_rng_max=120*4, minpairs=3000, nsub_min=200,
        nstratiform=0, nvol, cepoch, pepoch, timevol[20], np[20], ivol, n, *nofpair, nzi;
    float rho_thresh=0.98, zfloat, zmin=20, zmax=50, zmaxs=38, val_zmax,
          minlocal, maxlocal, slope, alpha;
    float *zgroup, *zdrgroup, *zdrtemp, *zmean, *zdrmedian,
          ftemp1, ftemp2;

    /* arrays */
    zgroup=get1dFloatArray(MAXIMUM_AZIMUTH_SUPER*MAXIMUM_RANGE_SUPER*5);
    zdrgroup=get1dFloatArray(MAXIMUM_AZIMUTH_SUPER*MAXIMUM_RANGE_SUPER*5);
    zdrtemp=get1dFloatArray(MAXIMUM_AZIMUTH_SUPER*MAXIMUM_RANGE_SUPER*5);
    zdrmedian=get1dFloatArray((zmax-zmin)/2+1);
    zmean=get1dFloatArray((zmax-zmin)/2+1);
    nofpair=get1dIntArray((zmax-zmin)/2+1);

    /* open and read the temporary file containg z-zdr pairs ***************/
    cepoch=getEpochTime(header);
    npair=0;
    ivol=0;
    if ((ftemp=fopen("alpha.tmp","r"))!=NULL) {
        fscanf(ftemp,"%d",&nvol);  // number of timestamp saved
        for (i=0; i<nvol; i++) {
             fscanf(ftemp,"%d %d",&pepoch,&n);
             if (cepoch-pepoch>0 && cepoch-pepoch<=HALFHOURINSEC) { // within 30-min
                 timevol[ivol]=pepoch;
                 np[ivol]=n;
                 ivol++;
                 for (j=0; j<n; j++) {
                      fscanf(ftemp,"%f %f",&ftemp1,&ftemp2);
                      zgroup[npair]=ftemp1;
                      zdrgroup[npair]=ftemp2;
                      npair++;
                 }
             }
             else { // not within 30-min
                 for (j=0; j<n; j++) {
                      fscanf(ftemp,"%f %f",&ftemp1,&ftemp2);
                 }
             }
        }
        fclose(ftemp);
    }
    /************************************************************************/
    /* get Z-Zdr pairs for 20<=Z<=50 */
    cpair=0;
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=sample_rng_min; rng<sample_rng_max; rng++) {
              if (melting[azi][rng][0]!=MELTING_FLAG && melting[azi][rng][0]!=FREEZING_FLAG) {
                  if ((float)rho[azi][rng][0]/scale>rho_thresh) {
                      if (z[azi][rng][0]!=NODATA && z[azi][rng][0]!=INITIAL_VALUE)
                          zfloat=int_dBZ2float_dBZ(z[azi][rng][0]);
                      else
                          continue;

                      if (zfloat>=zmin-1 && zfloat<=zmax+1) {
                          zgroup[npair]=zfloat;
                          zdrgroup[npair]=(float)zdr[azi][rng][0]/scale;
                          npair++;
                          cpair++;
                      }
                   }
              }
         }
    }
    printf("current pair: %d\n",cpair);

    /* write z-zdr pairs **************************************/
    timevol[ivol]=cepoch;
    np[ivol]=cpair;
    ivol++;

    ftemp=fopen("alpha.tmp","w");
    if (cpair>=minpairs) { // store current pairs only
        fprintf(ftemp,"%d\n",1);
        fprintf(ftemp,"%d %d\n",cepoch,cpair);
        for (j=npair-cpair; j<npair; j++) {
             fprintf(ftemp,"%f %f\n",zgroup[j],zdrgroup[j]);
        }
    }
    else {
        fprintf(ftemp,"%d\n",ivol);
        npair=0;
        for (i=0; i<ivol; i++) {
             fprintf(ftemp,"%d %d\n",timevol[i],np[i]);
             for (j=0; j<np[i]; j++) {
                  fprintf(ftemp,"%f %f\n",zgroup[npair],zdrgroup[npair]);
                  npair++;
             }
        }
    }
    fclose(ftemp);
/***************************************************************/

    if (npair<minpairs) {
        alpha=0.015;
        return alpha;
    }

    /* calulate the number of pairs used for current time */
/*
    int dn;
    dn=0;
    if (cpair>=minpairs)
        dn=cpair;
    else {
        for (i=ivol-1; i<ivol; i--) {
             dn=dn+np[i];
             if (dn>=minpairs)
                 break;
        }
    }
*/
    /* discretize Z and compute median Zdr */
    ni=0;
    nis=0; // for stratiform rain (20-40 dBZ)
    nstratiform=0;
    nzi=0;
    for (i=zmin; i<=zmax; i+=2) {
         nsub=0;
//         for (j=npair-dn; j<npair; j++) {  // for current scan
         for (j=0; j<npair; j++) {
              minlocal=i-1;
              maxlocal=i+1;
              if (zgroup[j]>=minlocal && zgroup[j]<maxlocal) {
                  zdrtemp[nsub]=zdrgroup[j];
                  nsub++;
              }
         }

         nofpair[nzi]=nsub;
         if (i<=zmaxs)
             nstratiform=nstratiform+nsub;

         if (nsub>nsub_min) {
             zmean[nzi]=(float)i;
             sortAscending(zdrtemp,nsub);
             zdrmedian[nzi]=zdrtemp[nsub/2];
             printf("Z: %f, Zdr: %f, n: %d\n",zmean[nzi],zdrmedian[nzi],nsub);
             ni++;
             if (i>=zmin && i<=zmaxs)
                 nis++;
         }
         nzi++;
    }

    /* find maximum Z */
    val_zmax=(float)0.;
    for (i=0; i<npair; i++) {
         if (zgroup[i]>val_zmax && zgroup[i]<=zmax)
             val_zmax=zgroup[i];
    }

    free(zgroup);
    free(zdrgroup);
    free(zdrtemp);

    /* compute slope and alpha
    method 1: no stratiform classification with bi-linear
           2:                              with nonlinear (quadratic)
           3:                              with nonlinear (fully)
           4: with stratiform classification with bi-linear
           5:                                with nonlinear (quadratic)
           6:                                with nonlinear (fully)
    */
    int method;
    double a1, a2, a3;

    method=1;  // need to define method

    if (method==1 || method==2 || method==3) {
    /* no stratiform calssification */
        if (ni<8)
            alpha=0.015;
        else {
            slope=getLeastSquareSlope(zmean,zdrmedian,(zmax-zmin)/2+1);
            if (method==1) {
                if (slope<=0.045)
                    alpha=0.04875-0.75*slope;
                else
                    alpha=0.015;
            }
            else if (method==2)
                alpha=0.054-1.31*slope+10.9*slope*slope;
            else if (method==3) {
                a1=0.530-0.0188*val_zmax+0.000194*val_zmax*val_zmax;
                a2=25.0-0.983*val_zmax+0.0103*val_zmax*val_zmax;
                a3=303-12.5*val_zmax+0.133*val_zmax*val_zmax;
                alpha=a1-a2*slope+a3*slope*slope;
            }
        }
    }
    else if (method==4 || method==5 || method==6) {
    /* with stratiform calssification */
        if (ni>=13) {
            slope=getLeastSquareSlope(zmean,zdrmedian,(zmax-zmin)/2+1);
            /* slope=getWeightedLeastSquareSlope(zmean,zdrmedian,nofpair,(zmax-zmin)/2+1); -> was not effective */
            if (method==4) {
                if (slope<=0.045)
                    alpha=0.04875-0.75*slope;
                else
                    alpha=0.015;
            }
            else if (method==5)
                alpha=0.054-1.31*slope+10.9*slope*slope;
            else if (method==6) {
                a1=0.530-0.0188*val_zmax+0.000194*val_zmax*val_zmax;
                a2=25.0-0.983*val_zmax+0.0103*val_zmax*val_zmax;
                a3=303-12.5*val_zmax+0.133*val_zmax*val_zmax;
                alpha=a1-a2*slope+a3*slope*slope;
            }
        }
        else if (nis>=8) {
            slope=getLeastSquareSlope(zmean,zdrmedian,(zmaxs-zmin)/2+1);
            if (method==4) {
                if (slope<=0.045)
                    alpha=0.04875-0.75*slope;
                else
                    alpha=0.035;
            }
            else if (method==5)
                alpha=0.054-1.31*slope+10.9*slope*slope;
            else if (method==6) {
                a1=0.530-0.0188*val_zmax+0.000194*val_zmax*val_zmax;
                a2=25.0-0.983*val_zmax+0.0103*val_zmax*val_zmax;
                a3=303-12.5*val_zmax+0.133*val_zmax*val_zmax;
                alpha=a1-a2*slope+a3*slope*slope;
            }
        }
        else {
            if (nstratiform>minpairs*0.9)
                alpha=0.035;
            else
                alpha=0.015;
        }
    }

    if (alpha<0.010 || alpha>0.040)
        alpha=0.015;

    //printf("val_zmax: %f, a1: %f, a2: %f, a3: %f\n",val_zmax, a1, a2, a3);
    printf("slope: %f, alpha: %f\n",slope, alpha);
    free(zmean);
    free(zdrmedian);
    free(nofpair);

    return alpha;
}

void sortAscending(float *array,
                   int n)
{
    int i, j;
    float a;

    for (i=0; i<n; i++) {
         for (j=i+1; j<n; j++) {
              if (array[i]>array[j]) {
                  a=array[i];
                  array[i]=array[j];
                  array[j]=a;
               }
         }
    }
}

float getLeastSquareSlope(float *x,
                          float *y,
                          int n)
{
    int i, count;
    float sumx, sumy, xbar, ybar, sxy, sxx, slope;

    sumx=0;
    sumy=0;
    count=0;
    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE) {
             sumx=x[i]+sumx;
             sumy=y[i]+sumy;
             count=count+1;
         }
    }
    if (count>1) {
        xbar=sumx/count;
        ybar=sumy/count;
    }
    else
        return NODATA;

    sxy=0;
    sxx=0;
    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE) {
             sxy=(x[i]-xbar)*(y[i]-ybar)+sxy;
             sxx=(x[i]-xbar)*(x[i]-xbar)+sxx;
         }
    }
    slope=sxy/sxx;

    return slope;
}

float getWeightedLeastSquareSlope(float *x,
                                  float *y,
                                  int *nofpair,
                                  int n)

{
    int i, sum, count;
    float sumx, sumy, xbar, ybar, sxy, sxx, slope, w[40];

    /* weight computation */
    sum=0;
    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE)
             sum=nofpair[i]+sum;
    }

    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE)
             w[i]=(float)nofpair[i]/sum;
    }

    sumx=0;
    sumy=0;
    count=0;
    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE) {
             sumx=x[i]+sumx;
             sumy=y[i]+sumy;
             count=count+1;
         }
    }
    if (count>1) {
        xbar=sumx/count;
        ybar=sumy/count;
    }
    else
        return NODATA;

    sxy=0;
    sxx=0;
    for (i=0; i<n; i++) {
         if (x[i]!=NODATA && x[i]!=INITIAL_VALUE) {
             sxy=w[i]*(x[i]-xbar)*(y[i]-ybar)+sxy;
             sxx=w[i]*(x[i]-xbar)*(x[i]-xbar)+sxx;
         }
    }
    slope=sxy/sxx;

    return slope;
}

void copy3Darray(int ***in,
                 int ***out)
{
    int ang, azi, rng;

    for (ang=0; ang<MAXIMUM_ANGLE_SUPER; ang++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                   out[azi][rng][ang]=in[azi][rng][ang];
              }
         }
    }
}

float int_dBZ2float_dBZ(short int_dBZ)
{
    float float_dBZ;

    if (int_dBZ != NODATA)
        float_dBZ=(float)(int_dBZ-66)/2;
    else
        float_dBZ=NODATA;

    return float_dBZ;
}

short float_dBZ2int_dBZ(float float_dBZ)
{
    short int_dBZ;

    if (float_dBZ != NODATA)
        int_dBZ=(short)(2*float_dBZ+66+0.5);
    else
        int_dBZ=NODATA;

    return int_dBZ;
}

float getSystemPhidp(int ***phi,
                     int scale,
                     int ***flag)
{
    int azi, rng;
    float phi_min, phi_val, *minphi;

    /* array */
    minphi=get1dFloatArray(MAXIMUM_AZIMUTH_SUPER);

    /* get mininum phidp for each ray */
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         phi_min=10000.;
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (flag[azi][rng][0]==RAIN) {// only for rain
                  if (phi[azi][rng][0]!=NODATA && phi[azi][rng][0]!=INITIAL_VALUE) {
                      phi_val=(float)phi[azi][rng][0]/scale;
                      if (phi_val<phi_min)
                          phi_min=phi_val;
                  }
              }
         }
         minphi[azi]=phi_min;
    }

    /* sort minmun values with descending order*/
    sortAscending(minphi,MAXIMUM_AZIMUTH_SUPER);

    phi_min=minphi[14]; /* median of the lowest 30 values */
    if (phi_min>1000)
        phi_min=0;

    free(minphi);

    return phi_min;
}

void averageProduct(float **rr,
                    int dw)
{
    int i, j;
    float **temp;

    /* define temporary array */
    temp=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    for (i=0; i<MAXIMUM_AZIMUTH_SUPER; i++) {
         for (j=0; j<MAXIMUM_RANGE_SUPER; j++) {
              if (rr[i][j]!=NODATA)
                  temp[i][j]=getAverage(rr,i,j,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,dw);
         }
    }

    /* copy and return averaged product */
    for (i=0; i<MAXIMUM_AZIMUTH_SUPER; i++) {
         for (j=0; j<MAXIMUM_RANGE_SUPER; j++) {
              rr[i][j]=temp[i][j];
         }
    }
}

float getAverage(float **rr,
                 int azi,
                 int rng,
                 int max_azi,
                 int max_rng,
                 int dw)
{
    float sum, avg;
    int count, i, j, iazi, jrng;

    sum=0;
    count=0;
    for (i=azi-dw; i<=azi+dw; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+max_azi;
         if (iazi>=max_azi)
             iazi=iazi-max_azi;

         for (j=rng-dw; j<=rng+dw; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=max_rng)
                  continue;

              if (rr[iazi][jrng]!=NODATA && rr[iazi][jrng]!=INITIAL_VALUE) {
                  sum=sum+rr[iazi][jrng];
                  count=count+1;
              }
         }
    }

    //if (count>(dw+dw+1)*(dw+dw+1)*0.3) {
    if (count==(dw+dw+1)*(dw+dw+1)) {
        avg=sum/count;
    }
    else
        avg=rr[iazi][jrng];

    return avg;
}

int getEpochTime(struct volume_info header)
{
    struct tm t;

    /* compute epochtime */
    t.tm_year=header.year-1900;
    t.tm_mon=header.mon-1;
    t.tm_mday=header.day;
    t.tm_hour=header.hour;
    t.tm_min=header.min;
    t.tm_sec=header.sec;
    t.tm_isdst=-1;
    time_t epochtime=timegm(&t);

    return epochtime;
}

double computeBeamHeight(double angle, double slant_range)
{
    /* Beam height calculation (May 2006)

       h=(R*sin(PHI))+(R^2/2*IR*Re),

       h: radar beam height for a given elevation angle (km)
       R: slant range observed on radar (km)
       PHI: radar elevation angle
       IR: refractive index (1.21)
       Re: radius of the earth (6371 km)

  	# Source: (http://www.wdtb.noaa.gov/tools/misc/beamwidth/index.html) */

    double rad_angle, beam_height;

    rad_angle=angle*RAD;
	beam_height=slant_range*sin(rad_angle)+
                slant_range*slant_range/(2*REFRACTIVE_INDEX*EARTH_RADIUS);
    return beam_height;
}




