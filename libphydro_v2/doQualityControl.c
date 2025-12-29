/*************************************************************************
--> DoQualityControl.c

    : Manages radar data QC (quality control) and Phidp interpolation (for specific attenuation)

    # QC
    1. basic QC by rho (>0.85)
    2. QC by conditional thresholds (a high rho threshold (0.98) has been disabled)

    # Phidp interpolation
    1. hydrometeor classification (defineEcho)
    2. filling gaps and interpolation of Phidp

*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

void startQualityControl(struct volume_info header,
                         short ***z,
                         int ***rho,
                         int ***zdr,
                         int ***phi,
                         int ***echo,
                         int angle)
{
    int ang, nazi, azi, rng, avg_rho, std_phidp,
        rho_thresh_high=0.98*DP_SCALEING, rho_thresh_low=0.90*DP_SCALEING, rho_bad=0.85*DP_SCALEING,
        rho_perfect=1.00*DP_SCALEING, dbZ_higher=35, phidp_thresh=10*DP_SCALEING;
    float coverage;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    /* basic QC by rho
       remove when rho > 1 or rho <= 0.85 */
    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]==NODATA) {
                  echo[azi][rng][ang]=DONTKNOW;
                  continue;
              }
              if (rho[azi][rng][ang]==INITIAL_VALUE)
                  echo[azi][rng][ang]=NORAIN;
              else if (rho[azi][rng][ang]<=rho_bad || rho[azi][rng][ang]>rho_perfect) {
                  z[azi][rng][ang]=INITIAL_VALUE;
                  zdr[azi][rng][ang]=INITIAL_VALUE;
                  rho[azi][rng][ang]=INITIAL_VALUE;
                  phi[azi][rng][ang]=INITIAL_VALUE;
                  echo[azi][rng][ang]=CLUTTER;
              }
         }
    }

    /* conditinally mask Zh using RHO values for AP removal
       remove when ZH >= 35 dBZ & RHO < 0.98
       remove when ZH < 35 dBZ,  remove RHO < 0.90 and Std(PHIDP) > 10  */
    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]==NODATA)
                  continue;
              avg_rho=getMean(rho,azi,rng,ang,nazi);
              //avg_rho=getMedian(rho,azi,rng,ang,nazi);

              if (avg_rho==NODATA) {
                  echo[azi][rng][ang]=NORAIN;
              }
              else {
                  if (z[azi][rng][ang]>=float_dBZ2int_dBZ((float)dbZ_higher)) {
                      if (avg_rho<rho_thresh_high) {
                          coverage=getCoverage(z,rho,azi,rng,ang,nazi);
                          if (coverage<0.5)
                              echo[azi][rng][ang]=CLUTTER;
                          /*
                          std_phidp=getStd(phi,azi,rng,ang,nazi);
                          if (std_phidp>phidp_thresh)
                              echo[azi][rng][ang]=CLUTTER;
                          */
                      }
                  }
                  else if (z[azi][rng][ang]<float_dBZ2int_dBZ((float)dbZ_higher)) {
                      if (avg_rho<rho_thresh_low) {
                          std_phidp=getStd(phi,azi,rng,ang,nazi);
                          if (std_phidp>phidp_thresh)
                              echo[azi][rng][ang]=CLUTTER;
                      }
                  }
              }
         }
    }
}

float getCoverage(short ***z,
                  int ***rho,
                  int azi,
                  int rng,
                  int ang,
                  int nazi)
{
    int wbin, wazi, w_km; // window size
    int count, i, j, iazi, jrng, avg, rho_bad=0.70*DP_SCALEING;
    float R, theta, z_thresh=10., cover;

    w_km=5;
    wbin=((float)w_km/2)/0.25; // 10 km
    R=0.25*(rng+1); // current distance from the radar
    theta=M_PI*0.5/180; // half beam width (radian)
//    wazi=(int)((int)(R*theta/10)*0.5+0.5); // adaptive window size in azimuth direction
    wazi=(int)(0.5*w_km/(R*theta)+0.5); // adaptive window size in azimuth direction

    count=0;
    for (i=azi-wazi; i<=azi+wazi; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+nazi;
         if (iazi>=nazi)
             iazi=iazi-nazi;

         for (j=rng-wbin; j<=rng+wbin; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (z[iazi][jrng][ang]!=NODATA && z[iazi][jrng][ang]!=INITIAL_VALUE && rho[iazi][jrng][ang]>rho_bad) {
                  if (int_dBZ2float_dBZ(z[iazi][jrng][ang])>z_thresh)
                      count++;
              }
         }
    }

    if (count>0)
        cover=(float)count/((wbin+wbin+1)*(wazi+wazi+1));
    else
        cover=0.;

    return cover;
}

int getMedian(int ***polargrid,
              int azi,
              int rng,
              int ang,
              int nazi)
{
    int swin; // window size
    int count, i, j, iazi, jrng, avg;
    float *temp;

    swin=2;

    temp=get1dFloatArray((swin*2+1)*(swin*2+1));
    count=0;
    for (i=azi-swin; i<=azi+swin; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+nazi;
         if (iazi>=nazi)
             iazi=iazi-nazi;

         for (j=rng-swin; j<=rng+swin; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (polargrid[iazi][jrng][ang]!=NODATA) {
                  temp[count]=polargrid[iazi][jrng][ang]/DP_SCALEING;
                  count=count+1;
              }
         }
    }
    if (count>(swin+swin+1)*(swin+swin+1)*0.3) {
        sortAscending(temp,count);
        avg=temp[(int)(count*0.5)]*DP_SCALEING;
    }
    else
        avg=NODATA;

    free(temp);

    return avg;
}

int getMean(int ***polargrid,
            int azi,
            int rng,
            int ang,
            int nazi)
{
    int swin; // window size
    int sum, count, i, j, iazi, jrng, avg;

    swin=2;
    sum=0;
    count=0;
    for (i=azi-swin; i<=azi+swin; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+nazi;
         if (iazi>=nazi)
             iazi=iazi-nazi;

         for (j=rng-swin; j<=rng+swin; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (polargrid[iazi][jrng][ang]!=NODATA) {
                  sum=sum+polargrid[iazi][jrng][ang];
                  count=count+1;
              }
         }
    }

    if (count>(swin+swin+1)*(swin+swin+1)*0.3)
        avg=(int)((float)sum/count+0.5);
    else
        avg=NODATA;

    return avg;
}

int getStd(int ***polargrid,
           int azi,
           int rng,
           int ang,
           int nazi)
{
    int swin; // window size
    int sum, count, i, j, iazi, jrng, std;
    int gvalue[100];
    float avg, fsum, distance;

    swin=2;
    sum=0;
    count=0;
    for (i=azi-swin; i<=azi+swin; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+nazi;
         if (iazi>=nazi)
             iazi=iazi-nazi;

         for (j=rng-swin; j<=rng+swin; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (polargrid[iazi][jrng][ang]!=NODATA) {
                  gvalue[count]=polargrid[iazi][jrng][ang];
                  sum=sum+polargrid[iazi][jrng][ang];
                  count=count+1;
              }
         }
    }

    if (count>(swin+swin+1)*(swin+swin+1)*0.3) {
        avg=(float)sum/count;

        fsum=0;
        for (i=0; i<count; i++) {
             distance=(gvalue[i]-avg)*(gvalue[i]-avg);
             fsum=fsum+distance;
        }
        std=(int)(sqrt((double)fsum/count)+0.5);
    }
    else
        std=NODATA;

    return std;
}

void smoothClutter(struct volume_info header,
                   int ***echo,
                   short ***z,
                   int angle)

{
    int i, j, ang, azi, rng, nazi, **temp_echo, rain_w, count;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    /* define new array */
    temp_echo=get2dIntArray(nazi,MAXIMUM_RANGE_SUPER);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              //temp_echo[azi][rng]=defineClutter(echo,z[azi][rng][ang],azi,rng,ang); // do not smooth clutter
              temp_echo[azi][rng]=echo[azi][rng][ang];
         }
    }

    /* remove isolated rain echo along rays */
    rain_w=4;  // 1 km
    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (temp_echo[azi][rng]==RAIN) {
                  count=1;
                  for (i=rng+1; i<MAXIMUM_RANGE_SUPER; i++) { // count consequtive rain cells
                       if (temp_echo[azi][i]==RAIN)
                           count+=1;
                       else
                           break;
                  }
                  if (count<=rain_w) {
                      for (j=i-1; j>=i-count; j--) {
                           temp_echo[azi][j]=CLUTTER;
                      }
                  }
                  else
                      rng=i;
              }
         }
    }

    /* remove isolated clutter echo along rays */
    for (azi=0; azi<nazi; azi++) {
         for (rng=1; rng<MAXIMUM_RANGE_SUPER-3; rng++) {
              if (temp_echo[azi][rng]==CLUTTER) {
                  if (temp_echo[azi][rng-1]==RAIN) {
                      if (temp_echo[azi][rng+1]==RAIN || temp_echo[azi][rng+2]==RAIN || temp_echo[azi][rng+3]==RAIN)
                          temp_echo[azi][rng]=RAIN;
                  }
              }
         }
    }

    /* copy array */
    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              echo[azi][rng][ang]=temp_echo[azi][rng];
         }
    }
}

int defineClutter(int ***echo,
                  short zint,
                  int azi,
                  int rng,
                  int ang)
{
    int swin=2, count, i, j, iazi, jrng, avg;

    if (echo[azi][rng][ang]==CLUTTER) {
        count=0;
        for (i=azi-swin; i<=azi+swin; i++) {
             iazi=i;
             if (iazi<0)
                 iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
             if (iazi>=MAXIMUM_AZIMUTH_SUPER)
                 iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

             for (j=rng-swin; j<=rng+swin; j++) {
                  jrng=j;
                  if (jrng<0)
                      continue;
                  if (jrng>=MAXIMUM_RANGE_SUPER)
                      continue;

                  if (echo[iazi][jrng][ang]==CLUTTER)
                      count+=1;
             }
        }

        if (count>(swin*2+1)*(swin*2+1)*0.3)
            return CLUTTER;
        else {
            if (zint!=NODATA && zint!=INITIAL_VALUE) {
                if (int_dBZ2float_dBZ(zint)>0)
                    return RAIN;
                else
                    return NORAIN;
            }
            else
                return DONTKNOW;
        }
    }
    else
        return echo[azi][rng][ang];
}

void maskObservables(struct volume_info header,
                     int ***echo,
                     short ***z,
                     int ***rho,
                     int ***zdr,
                     int ***phi,
                     int angle)
{
    int ang, azi, rng, nazi;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (echo[azi][rng][ang]==CLUTTER) {
                  z[azi][rng][ang]=INITIAL_VALUE;
                  zdr[azi][rng][ang]=INITIAL_VALUE;
                  rho[azi][rng][ang]=INITIAL_VALUE;
                  phi[azi][rng][ang]=INITIAL_VALUE;
              }
         }
    }
}

void defineEcho(struct volume_info header,
                short ***z,
                int ***rho,
                int ***zdr,
                float ***T_top,
                float ***T_bottom,
                int ***flag,
                short ***melting,
                short ***rtype,
                int angle,
                int scale)
{
    /* unknown: -1, rain: 0, clutter: 1, norain: 2, hail: 3,
       dry snow: 4, wet snow: 5, snowice: 6 (dry snow above melting layer), ice: 7 */
    int ang, azi, rng, nazi, rng_nodata=8;
    short flag_type;
    float zfloat, zdrfloat, rhofloat, znorain=0., avgT;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (rng<rng_nodata) {
                  flag[azi][rng][ang]=DONTKNOW;
                  continue;
              }
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE) {
                  zfloat=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  if (zfloat<=znorain) {
                      flag[azi][rng][ang]=NORAIN;
                      continue;
                  }
                  zdrfloat=(float)zdr[azi][rng][ang]/scale;
                  rhofloat=(float)rho[azi][rng][ang]/scale;
                  flag_type=rtype[azi][rng][ang];
                  avgT=0.5*(T_top[azi][rng][ang]+T_bottom[azi][rng][ang]);

                  if (zfloat>=ZHAIL) {// hail
                  //if (zfloat>=zhail && zdrfloat<zdrhail) /* zdr value is not reliable (zdrhail was 0.5 dB) */
                      flag[azi][rng][ang]=HAIL;
                      continue;
                  }

                  if (melting[azi][rng][ang]==MELTING_FLAG || melting[azi][rng][ang]==FREEZING_FLAG) {
                      if (zfloat>znorain && zfloat<=20. && flag_type!=CONVECTIVE) {
                          if (zdrfloat>=1.) { // ice
                              if (avgT<MELTING_THRESHOLD_BELOW) // ice (t<5)
                                  flag[azi][rng][ang]=ICE;
                          }
                          else {
                              if (avgT<=MELTING_THRESHOLD_ABOVE) // ice (t<=-5)
                                  flag[azi][rng][ang]=ICE;
                              else if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW) // dry snow (-5<t<5)
                                  flag[azi][rng][ang]=DRYSNOW;
                          }
                      }
                      else if (zfloat>20. && zfloat<=25. && flag_type!=CONVECTIVE) {
                          if (zdrfloat>=1.) {
                              if (avgT<=MELTING_THRESHOLD_ABOVE) // ice (t<=-5)
                                  flag[azi][rng][ang]=ICE;
                              else if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW) // dry snow (-5<t<5)
                                  flag[azi][rng][ang]=DRYSNOW;
                          }
                          else {
                              if (avgT<=MELTING_THRESHOLD_ABOVE) // dry snow above melting layer
                                  flag[azi][rng][ang]=SNOWICE;
                              else if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW) // dry snow (-5<t<5)
                                  flag[azi][rng][ang]=DRYSNOW;
                          }
                      }
                      else if (zfloat>25. && zdrfloat<2. && flag_type!=CONVECTIVE) {
                          if (rhofloat>=0.95) {
                              if (avgT<=MELTING_THRESHOLD_ABOVE) {// dry snow above melting layer
                                  flag[azi][rng][ang]=SNOWICE;
                              }
                              else if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW) // dry snow (-5<t<5)
                                  flag[azi][rng][ang]=DRYSNOW;
                          }
                          else if (rhofloat>0.9 && rhofloat<0.95) {
                              if (avgT<=MELTING_THRESHOLD_ABOVE) // dry snow above melting layer
                                  flag[azi][rng][ang]=SNOWICE;
                              else if (avgT>MELTING_THRESHOLD_ABOVE && avgT<0) // dry snow (-5<t<0)
                                  flag[azi][rng][ang]=DRYSNOW;
                              else if (avgT>=0 && avgT<MELTING_THRESHOLD_BELOW) // wet snow (0=<t<5)
                                  flag[azi][rng][ang]=WETSNOW;
                          }
                          else if (rhofloat<=0.9)
                              flag[azi][rng][ang]=WETSNOW;
                      }
                  }
              }
              else if (z[azi][rng][ang]==INITIAL_VALUE && flag[azi][rng][ang]!=CLUTTER)
                  flag[azi][rng][ang]=NORAIN;
              else if (z[azi][rng][ang]==INITIAL_VALUE && flag[azi][rng][ang]==CLUTTER) // clutter has been replaced with 0
                  flag[azi][rng][ang]=CLUTTER;
              else if (z[azi][rng][ang]==NODATA)
                  flag[azi][rng][ang]=DONTKNOW;
         }
    }
}

void fillGaps(struct volume_info header,
              int ***phi,
              int ***flag,
              int w,
              int angle,
              int scale)
{
    int ang, azi, rng, nazi, ng;
    int *class_ray, *flag_ray, *ibegin, *iend;
    float *phi_ray, *temp_phi;

    /* arrays */
    phi_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    temp_phi=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    class_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    flag_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    ibegin=get1dIntArray(MAXIMUM_RANGE_SUPER);
    iend=get1dIntArray(MAXIMUM_RANGE_SUPER);

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              class_ray[rng]=flag[azi][rng][ang];
              flag_ray[rng]=flag[azi][rng][ang];
              if (phi[azi][rng][ang]!=NODATA)
                  phi_ray[rng]=(float)phi[azi][rng][ang]/scale;
              else if (phi[azi][rng][ang]==NODATA)
                  phi_ray[rng]=NODATA;
         }

         /* count the number of the meteorological block */
         ng=countMetGroups(flag_ray,ibegin,iend);

         /* interpolate */
         if (ng>0)
             interpolatePhi(phi_ray,temp_phi,class_ray,flag_ray,w,ng,ibegin,iend);
/*
if (azi==404) {
    for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
         printf("rng: %d, flag: %d, class: %d, phi: %f, phi_new: %f\n",rng,flag_ray[rng],class_ray[rng],phi_ray[rng],temp_phi[rng]);
    }
}
*/
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              phi[azi][rng][ang]=(int)(temp_phi[rng]*scale+0.5);
         }
    }
    free(phi_ray);
    free(temp_phi);
    free(flag_ray);
    free(ibegin);
    free(iend);
}

int countMetGroups(int *flag,
                   int *ibegin,
                   int *iend)
{
    int i, n, finish;

    /* unknown: -1, rain: 0, clutter: 1, norain: 2, hail: 3,
       dry snow: 4, wet snow: 5, snowice: 6 (dry snow above melting layer), ice: 7 */
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (i<8)
             flag[i]=NO_MET_ECHO;
         else {
             if (flag[i]==RAIN || flag[i]==HAIL || flag[i]==DRYSNOW || flag[i]==WETSNOW || flag[i]==SNOWICE || flag[i]==ICE)
                 flag[i]=MET_ECHO; // met echo
             else
                 flag[i]=NO_MET_ECHO;
         }
    }

    n=0;
    finish=0;
    if (flag[0]==MET_ECHO) {
        ibegin[n]=0;
        finish=1;
    }

    for (i=1; i<MAXIMUM_RANGE_SUPER; i++) {
         if (flag[i-1]==NO_MET_ECHO && flag[i]==MET_ECHO) {
             ibegin[n]=i;
             finish=1;
         }
         else if (flag[i-1]==MET_ECHO && flag[i]==NO_MET_ECHO) {
             iend[n]=i-1;
             finish=0;
             n++;
         }
    }

    if (finish) {
        iend[n]=MAXIMUM_RANGE_SUPER-1;
        n++;
    }

    return n;
}

void interpolatePhi(float *phi_ray,
                    float *new_ray,
                    int *class_ray,
                    int *flag_ray,
                    int w,
                    int ng,
                    int *ibegin,
                    int *iend)
{
    int i, j, n, hw, cnt, pre, mlength, count_flag;
    int beginx, endx;
    float phibegin, phiend, slope, phi_init, phi_min;

    /* get initial Phidp from rain areas */
    phi_min=1000.;
    count_flag=0;
    for (j=0; j<ng; j++) {
         mlength=iend[j]-ibegin[j]+1;
         if (mlength<w/2)
             continue;
         for (i=ibegin[j]; i<=iend[j]; i++) {
              count_flag=1;
              if (flag_ray[i]==MET_ECHO) {
                  if (phi_ray[i]!=NODATA && phi_ray[i]!=INITIAL_VALUE) {
                      if (phi_ray[i]<phi_min) {
                          phi_min=phi_ray[i];
                      }
                  }
              }
         }
         break; // stop at first Met segment
    }

    if (count_flag==0) { // if there is no sufficiently long segment
        for (j=0; j<ng; j++) {
             for (i=ibegin[j]; i<=iend[j]; i++) {
                  if (flag_ray[i]==MET_ECHO) {
                      if (phi_ray[i]!=NODATA && phi_ray[i]!=INITIAL_VALUE) {
                          if (phi_ray[i]<phi_min) {
                              phi_min=phi_ray[i];
                          }
                      }
                  }
             }
        }
    }

    if (phi_min<1000.)
        phi_init=phi_min;
    else
        phi_init=0.;

    hw=w/2;
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         new_ray[i]=phi_ray[i];
    }

    cnt=0;
    pre=0;
    for (i=0; i<=ng; i++) {
         if (i<ng) {
             mlength=iend[i]-ibegin[i]+1;
             if (mlength<w)
                 continue;

             if (cnt==0) { // first met group
                 for (j=0; j<=ibegin[i]; j++) {
                      new_ray[j]=phi_init;
                 }
                 pre=i;
                 cnt++;
                 continue;
             }
             else {
                 beginx=iend[pre]-hw;
                 endx=ibegin[i]+hw;
                 phibegin=phi_ray[beginx];
                 phiend=phi_ray[endx];
             }
         }
         else { // i==ng (will fill the last no echo block with constant phidp values)
             if (cnt==0)
                 break;
             beginx=iend[pre];
             endx=MAXIMUM_RANGE_SUPER-1;
             phibegin=new_ray[beginx];
             phiend=new_ray[beginx];
         }

         /* interpolat no echo block */
         if (endx-beginx>0) {
             slope=(phiend-phibegin)/(float)(endx-beginx);
             for (j=beginx; j<=endx; j++) {
                  new_ray[j]=slope*(j-beginx)+phibegin;
             }
         }
         pre=i;
         cnt++;
    }

    if (cnt==0) {
        for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
             new_ray[i]=phi_init;
        }
    }
}

