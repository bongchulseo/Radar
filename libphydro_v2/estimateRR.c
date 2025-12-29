/********************************************************************************************************
--> estimateRR.c

    : Functions for rain rate estimation

    1. R(Z): Z-R (conventional: 300 and 1.4)
    2. R(Z): Z-R with rain type classification
    3. R(Z, Zdr)
    4. R(Kdp)
    5. R(A)

Added:    getRayRA2 --> May 2020  five-point running median filtering (for different azimuths) in R(A)
                                  to mitigate radial patterns

Modified: corrected a way to combine rain and no rain segments in R(A) --> May 2020
          applies R(Z) if the farest pixel of the first rain segment is in 50 km --> June 2020
             to mitigate significatn radial features in R(A) near the radar site

***********************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

void estimateRA(struct volume_info header,
                short ***z,
                int ***phi,
                int ***kdp,
                int ***echo,
                short ***rtype,
                short ***melting,
                float **azimuth,
                float ***rr,
                float a,
                int angle,
                int scale)
{
    /* estimate rain rate using specific attenuation */
    short *melting_ray;
    int ang, azi, rng, nazi, *class_ray, *type_ray;
    float *z_ray, *phi_ray, *kdp_ray, *rate_ray;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    /* arrays */
    z_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    phi_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    kdp_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    rate_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    class_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    type_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    melting_ray=get1dShortArray(MAXIMUM_RANGE_SUPER);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              class_ray[rng]=echo[azi][rng][ang];
              type_ray[rng]=rtype[azi][rng][ang];
              melting_ray[rng]=melting[azi][rng][ang];

              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE)
                  z_ray[rng]=int_dBZ2float_dBZ(z[azi][rng][ang]);
              else
                  z_ray[rng]=z[azi][rng][ang];

              if (phi[azi][rng][ang]!=NODATA && phi[azi][rng][ang]!=INITIAL_VALUE)
                  phi_ray[rng]=(float)phi[azi][rng][ang]/scale;
              else
                  phi_ray[rng]=phi[azi][rng][ang];

              if (kdp[azi][rng][ang]!=NODATA && kdp[azi][rng][ang]!=INITIAL_VALUE)
                  kdp_ray[rng]=(float)kdp[azi][rng][ang]/scale;
              else
                  kdp_ray[rng]=kdp[azi][rng][ang];
         }
         initialize1dFloatArray(rate_ray,MAXIMUM_RANGE_SUPER,INITIAL_VALUE);

         //getRayRA(z_ray,phi_ray,kdp_ray,class_ray,type_ray,melting_ray,a,rate_ray);
         getRayRA2(z_ray,phi_ray,kdp_ray,class_ray,type_ray,melting_ray,a,rate_ray,phi,echo,ang,azi,nazi);  // five-point running median filter

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              rr[azi][rng][ang]=rate_ray[rng];
         }
    }
    free(z_ray);
    free(phi_ray);
    free(kdp_ray);
    free(class_ray);
    free(type_ray);
    free(rate_ray);
}

void estimateRAwVPR(struct volume_info header,
                    struct range_parameter evpr,
                    short ***z,
                    int ***phi,
                    int ***kdp,
                    int ***echo,
                    short ***rtype,
                    short ***melting,
                    float **vpr,
                    float **azimuth,
                    float a,
                    float *eangle,
                    int angle,
                    int scale,
                    float ***rr)
{
    /* estimate rain rate using specific attenuation */
    short *melting_ray;
    int ang, azi, rng, nazi, *class_ray, *type_ray;
    float *z_ray, *phi_ray, *kdp_ray, *rate_ray;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    /* arrays */
    z_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    phi_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    kdp_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    rate_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);
    class_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    type_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    melting_ray=get1dShortArray(MAXIMUM_RANGE_SUPER);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              class_ray[rng]=echo[azi][rng][ang];
              type_ray[rng]=rtype[azi][rng][ang];
              melting_ray[rng]=melting[azi][rng][ang];

              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]!=INITIAL_VALUE)
                  z_ray[rng]=int_dBZ2float_dBZ(z[azi][rng][ang]);
              else
                  z_ray[rng]=z[azi][rng][ang];

              if (phi[azi][rng][ang]!=NODATA && phi[azi][rng][ang]!=INITIAL_VALUE)
                  phi_ray[rng]=(float)phi[azi][rng][ang]/scale;
              else
                  phi_ray[rng]=phi[azi][rng][ang];

              if (kdp[azi][rng][ang]!=NODATA && kdp[azi][rng][ang]!=INITIAL_VALUE)
                  kdp_ray[rng]=(float)kdp[azi][rng][ang]/scale;
              else
                  kdp_ray[rng]=kdp[azi][rng][ang];
         }
         initialize1dFloatArray(rate_ray,MAXIMUM_RANGE_SUPER,INITIAL_VALUE);
         getRayRAwVPR(evpr,z_ray,phi_ray,kdp_ray,class_ray,type_ray,melting_ray,vpr,eangle[ang],a,rate_ray);

         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              rr[azi][rng][ang]=rate_ray[rng];
         }
    }
    free(z_ray);
    free(phi_ray);
    free(kdp_ray);
    free(class_ray);
    free(type_ray);
    free(melting_ray);
    free(rate_ray);
}
void getRayRA2(float *z,
               float *phi,
               float *kdp,
               int *eclass,
               int *rtype,
               short *melt,
               float alpha,
               float *rate,
               int ***phi3d,
               int ***echo3d,
               int iang,
               int iazi,
               int nazi)
{
    int i, j, k, last, *seg_start, *seg_end, *seg_flag, *nseg_flag, nseg, iseg, istart, iend, count, raincount, binmeltstart;
    float tI, subI, b=0.62, Fphi, up, down, phi_thresh=3., A, phi_delta, C1;
    int azi, phival, echoval, dwphi=5, rng_segend=50;
	float *phitemp, phi_median_start, phi_median_end;

    /* arrays */
    seg_start=get1dIntArray(MAXIMUM_RANGE_SUPER);
    seg_end=get1dIntArray(MAXIMUM_RANGE_SUPER);
    phitemp=get1dFloatArray(dwphi); // for five-point median filtering

    /* find a bin where the melting layer starts */
    binmeltstart=0;
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (melt[i]==MELTING_FLAG || melt[i]==FREEZING_FLAG) {
             binmeltstart=i;
             break;
         }
    }

    /****** find rain segment that does not contain hail and clutter ******/
    /* find the start and end of rain segment */
    count=0; // need to be removed
    nseg=0;
    istart=NODATA;
    iend=NODATA;
    for (i=0; i<binmeltstart; i++) {
         if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
             count=1;
             istart=i;
             break;
         }
    }

    for (i=binmeltstart-1; i>=0; i--) {
         if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
             iend=i;
             break;
         }
    }

    /* start finding rain segements */
    if (istart!=NODATA && iend!=NODATA && istart==iend) {
        seg_start[nseg]=istart;
        seg_end[nseg]=iend;
        nseg++;
    }
    else if (istart!=NODATA && iend!=NODATA && istart!=iend) {
        seg_start[nseg]=istart;
        for (i=istart; i<iend; i++) {
             if (eclass[i]==RAIN && eclass[i+1]!=RAIN && phi[i]!=INITIAL_VALUE) {
                 seg_end[nseg]=i;
                 nseg++;
             }
             else if (eclass[i]!=RAIN && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE) {
                 seg_start[nseg]=i+1;
             }
        }
        seg_end[nseg]=iend;
        nseg++;
    }

    /* to combine rain segments with no echo regions */
    if (nseg>1) {
        seg_flag=get1dIntArray(nseg-1);
        for (i=0; i<nseg-1; i++) {
             for (j=seg_end[i]+1; j<seg_start[i+1]; j++) {
                  if (eclass[j]==CLUTTER || eclass[j]==HAIL) {
                  //if (eclass[j]==HAIL) {
                      seg_flag[i]=NODATA;
                      break;
                  }
             }
        }

        /* assign new rain segment number */
        nseg_flag=get1dIntArray(nseg);
        count=0;
        for (i=0; i<nseg-1; i++) {
             nseg_flag[i]=count;
             if (seg_flag[i]==NODATA)
                 count++;
        }
        nseg_flag[nseg-1]=count;

        /* adjuat segment start and end */
        for (i=0; i<nseg-1; i++) {
             if (nseg_flag[i]!=nseg_flag[i+1]) {
                 seg_end[nseg_flag[i]]=seg_end[i];
                 seg_start[nseg_flag[i+1]]=seg_start[i+1];

             }
        }
        seg_end[nseg_flag[nseg-1]]=seg_end[nseg-1];
        nseg=count+1;
    }

    /*******  previous approach for splitting or combining rain segments (removed in May 2020)
    if (count!=0) {
        for (i=binmeltstart-1; i>=0; i--) {
             if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
                 iend=i;
                 break;
             }
        }

        for (i=istart; i<iend; i++) {
             if (i==istart) {
                 if (phi[i]!=INITIAL_VALUE)
                     seg_start[nseg]=i;
                 continue;
             }

             if (eclass[i]==RAIN && eclass[i+1]!=RAIN && phi[i]!=INITIAL_VALUE) {
                 k=0;
                 raincount=0;
                 for (j=i+1; j<iend; j++) {
                      if (eclass[j]==NORAIN || eclass[j]==DONTKNOW || eclass[j]==RAIN) {
                          if (eclass[j]==RAIN) {
                              raincount++;
                              last=j;
                          }
                          k++;
                          if (j==iend-1)
                              i=j;
                      }
                      else if (eclass[j]==CLUTTER || eclass[j]==HAIL) {
                          if (k==0 || raincount==0)
                              seg_end[nseg]=i;
                          else
                              seg_end[nseg]=last;
                          nseg++;
                          i=j;
                          break;
                      }
                 }
             }

             if (i!=iend-1) {
                 if ((eclass[i]==HAIL || eclass[i]==CLUTTER || eclass[i]==NORAIN || eclass[i]==DONTKNOW) && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE)
                     seg_start[nseg]=i+1;
             }
             else {
                 if ((eclass[i]==HAIL || eclass[i]==CLUTTER) && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE)
                     seg_start[nseg]=i+1;
             }
        }
        */

        /* check the last bin */
        /*
        if (eclass[iend]!=HAIL) {
            seg_end[nseg]=iend;
            nseg++;
        }
    }
    *******************************************/

    /* start computation for path-integrated attenuation (PIA) */
    for (i=0; i<nseg; i++) {
         phi_delta=phi[seg_end[i]]-phi[seg_start[i]];
         if (phi_delta>=phi_thresh) { // then, apply 5-point median filtering
             /* get five-running points median */
             /* seg_start */
             initialize1dFloatArray(phitemp,dwphi,NODATA);
             count=0;
             for (j=-2; j<=2; j++) {
                  azi=iazi+j;
                  if (azi<0)
                      azi=azi+nazi;
                  else if (azi>=nazi)
                      azi=azi-nazi;

                  phival=phi3d[azi][seg_start[i]][iang];
                  echoval=echo3d[azi][seg_start[i]][iang];
                  if (phival!=NODATA && phival!=INITIAL_VALUE && echoval==RAIN) {
                      phitemp[count]=(float)phival/DP_SCALEING;
                      count++;
                  }
             }
             sortAscending(phitemp,count);
             phi_median_start=phitemp[(count-1)/2];

             /* seg_end */
             initialize1dFloatArray(phitemp,dwphi,NODATA);
             count=0;
             for (j=-2; j<=2; j++) {
                  azi=iazi+j;
                  if (azi<0)
                      azi=azi+nazi;
                  else if (azi>=nazi)
                      azi=azi-nazi;

                  phival=phi3d[azi][seg_end[i]][iang];
                  echoval=echo3d[azi][seg_end[i]][iang];
                  if (phival!=NODATA && phival!=INITIAL_VALUE && echoval==RAIN) {
                      phitemp[count]=(float)phival/DP_SCALEING;
                      count++;
                  }
             }
             sortAscending(phitemp,count);
             phi_median_end=phitemp[(count-1)/2];

             phi_delta=phi_median_end-phi_median_start;
         }

         if (phi_delta>=phi_thresh && seg_end[i]>rng_segend*KM/BIN_SUPER) { // then, use R(A)
         //if (phi_delta>=phi_thresh) { // then, use R(A)
             /*
             added "seg_end[i]>rng_segend*KM/BIN_SUPER"
             to mitigate unstable rainfall feature at close range (e.g. due to ground clutter)
             in June 2020
             */

             /* calculate total integral for the segment */
             tI=0.;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE)
                      tI=0.46*b*pow(pow(10.,z[j]/10.),b)*((float)BIN_SUPER/1000)+tI;
             }

             Fphi=exp(0.23*alpha*b*(phi_delta))-1;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  subI=0.;
                  for (k=j; k<=seg_end[i]; k++) {
                       if (z[k]!=NODATA && z[k]>INITIAL_VALUE)
                           subI=0.46*b*pow(pow(10.,z[k]/10.),b)*((float)BIN_SUPER/1000)+subI;
                  }

                  /* specific attenuation A(r) */
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE) {
                      up=pow(pow(10.,z[j]/10.),b)*Fphi;
                      down=tI+Fphi*subI;
                      A=up/down;
                      if (eclass[j]==RAIN)
                          rate[j]=4120*pow(A,1.03);
                  }
             }
         }
         else { // use R(Z)
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE) {
                      if (eclass[j]==CLUTTER)
                          rate[j]=0.;
                      else { // there is no hail in rain segments
                          if (rtype[j]==CONVECTIVE)
                              rate[j]=1.7007*pow(10,-2)*pow(10.,0.0714*z[j]); // NEXRAD
                          else
                              rate[j]=3.6463*pow(10,-2)*pow(10.,0.0625*z[j]); // M-P
                      }
                  }
             }
         }
    }
    free(seg_start);
    free(seg_end);
    /*--------------------------------------------------------------------------------------------*/

    /* within and above melting layer */
    for (i=binmeltstart; i<MAXIMUM_RANGE_SUPER; i++) {
         if (eclass[i]!=DONTKNOW && eclass[i]!=CLUTTER && eclass[i]!=NORAIN && eclass[i]!=HAIL) {
             if (eclass[i]==RAIN) {
                 if (rtype[i]==CONVECTIVE)
                     rate[i]=1.7007*pow(10.,-2.)*pow(10.,0.0714*z[i]); // NEXRAD
                 else
                     rate[i]=3.6463*pow(10.,-2.)*pow(10.,0.0625*z[i]); // M-P
             }
             else {
                 if (rtype[i]==CONVECTIVE)
                     C1=1.;
                 else {
                     if (eclass[i]==DRYSNOW)
                         C1=1.;
                     else if (eclass[i]==WETSNOW)
                         C1=0.6;
                     else if (eclass[i]==SNOWICE || eclass[i]==ICE)
                         C1=2.8;
                 }
                 rate[i]=C1*1.7007*pow(10,-2)*pow(10.,0.0714*z[i]);
             }
         }
         else if (eclass[i]==HAIL) {
             if (kdp[i]>0)
                 rate[i]=29*pow(kdp[i],0.77);
             else {
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*ZHAIL);
             }
         }
    }

    /* fill hail cells below melting layer */
    for (i=0; i<binmeltstart; i++) {
         if (eclass[i]==HAIL) {
             if (kdp[i]>0)
                 rate[i]=29*pow(kdp[i],0.77);
             else {
                 if (z[i]>ZHAIL)
                     z[i]=ZHAIL;
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*z[i]); // hail is an indication of convection
             }
         }
    }
}

void getRayRA(float *z,
              float *phi,
              float *kdp,
              int *eclass,
              int *rtype,
              short *melt,
              float alpha,
              float *rate)
{
    int i, j, k, last, *seg_start, *seg_end, nseg, istart, iend, count, binmeltstart;
    float tI, subI, b=0.62, Fphi, up, down, phi_thresh=3., A, phi_delta, C1;

    /* arrays */
    seg_start=get1dIntArray(MAXIMUM_RANGE_SUPER);
    seg_end=get1dIntArray(MAXIMUM_RANGE_SUPER);

    /* find a bin where the melting layer starts */
    binmeltstart=0;
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (melt[i]==MELTING_FLAG || melt[i]==FREEZING_FLAG) {
             binmeltstart=i;
             break;
         }
    }

    /****** find rain segment that does not contain hail ******/
    /* find the start and end of rain segment by phidp */
    count=0;
    for (i=0; i<binmeltstart; i++) {
         if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
             count=1;
             istart=i;
             break;
         }
    }

    nseg=0;
    if (count!=0) {
        for (i=binmeltstart-1; i>=0; i--) {
             if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
                 iend=i;
                 break;
             }
        }

        for (i=istart; i<iend; i++) {
             if (i==istart) {
                 if (phi[i]!=INITIAL_VALUE)
                     seg_start[nseg]=i;
                 continue;
             }

             if (eclass[i]==RAIN && eclass[i+1]!=RAIN && phi[i]!=INITIAL_VALUE) {
                 k=0;
                 for (j=i+1; j<iend; j++) {
                      if (eclass[j]==NORAIN || eclass[j]==DONTKNOW || eclass[j]==RAIN) {
                          k++;
                          last=j;
                      }
                      if (eclass[j]==CLUTTER || eclass[j]==HAIL) {
                          if (k==0)
                              seg_end[nseg]=i;
                          else
                              seg_end[nseg]=last;
                          nseg++;
                          i=j;
                          break;
                      }
                 }
                 /*
                 if (j==iend) {
                     seg_end[nseg]=iend;
                     nseg++;
                     break;
                 }
                 */
             }

             if ((eclass[i]==HAIL || eclass[i]==CLUTTER || eclass[i]==NORAIN || eclass[i]==DONTKNOW) && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE)
                 seg_start[nseg]=i+1;

/*
             if ((eclass[i]==HAIL || eclass[i]==CLUTTER) && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE)
                 seg_start[nseg]=i+1;

             if ((eclass[i]!=CLUTTER && eclass[i]!=HAIL) && (eclass[i+1]==CLUTTER || eclass[i+1]==HAIL) && phi[i]!=INITIAL_VALUE) {
                 seg_end[nseg]=i;
                 nseg++;
             }
*/
        }

        /* check the last bin */
        if (eclass[iend]!=HAIL) {
            seg_end[nseg]=iend;
            nseg++;
        }

    }

    /* start computation for path-integrated attenuation (PIA) */
    for (i=0; i<nseg; i++) {
         phi_delta=phi[seg_end[i]]-phi[seg_start[i]];
         if (phi_delta>=phi_thresh) { // then, use R(A)
             /* calculate total integral for the segment */
             tI=0.;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE)
                      tI=0.46*b*pow(pow(10.,z[j]/10.),b)*((float)BIN_SUPER/1000)+tI;
             }

             Fphi=exp(0.23*alpha*b*(phi_delta))-1;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  subI=0.;
                  for (k=j; k<=seg_end[i]; k++) {
                       if (z[k]!=NODATA && z[k]>INITIAL_VALUE)
                           subI=0.46*b*pow(pow(10.,z[k]/10.),b)*((float)BIN_SUPER/1000)+subI;
                  }

                  /* specific attenuation A(r) */
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE) {
                      up=pow(pow(10.,z[j]/10.),b)*Fphi;
                      down=tI+Fphi*subI;
                      A=up/down;
                      if (eclass[j]==RAIN)
                          rate[j]=4120*pow(A,1.03);
                  }
             }
         }
         else { // use R(Z)
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]>INITIAL_VALUE) {
                      if (eclass[j]==CLUTTER)
                          rate[j]=0.;
                      else { // there is no hail in rain segments
                          if (rtype[j]==CONVECTIVE)
                              rate[j]=1.7007*pow(10,-2)*pow(10.,0.0714*z[j]); // NEXRAD
                          else
                              rate[j]=3.6463*pow(10,-2)*pow(10.,0.0625*z[j]); // M-P
                      }
                  }
             }
         }
    }
    free(seg_start);
    free(seg_end);
    /*--------------------------------------------------------------------------------------------*/

    /* within and above melting layer */
    for (i=binmeltstart; i<MAXIMUM_RANGE_SUPER; i++) {
         if (eclass[i]!=DONTKNOW && eclass[i]!=CLUTTER && eclass[i]!=NORAIN && eclass[i]!=HAIL) {
             if (eclass[i]==RAIN) {
                 if (rtype[i]==CONVECTIVE)
                     rate[i]=1.7007*pow(10.,-2.)*pow(10.,0.0714*z[i]); // NEXRAD
                 else
                     rate[i]=3.6463*pow(10.,-2.)*pow(10.,0.0625*z[i]); // M-P
             }
             else {
                 if (rtype[i]==CONVECTIVE)
                     C1=1.;
                 else {
                     if (eclass[i]==DRYSNOW)
                         C1=1.;
                     else if (eclass[i]==WETSNOW)
                         C1=0.6;
                     else if (eclass[i]==SNOWICE || eclass[i]==ICE)
                         C1=2.8;
                 }
                 rate[i]=C1*1.7007*pow(10,-2)*pow(10.,0.0714*z[i]);
             }
         }
         else if (eclass[i]==HAIL) {
             if (kdp[i]>0)
                 rate[i]=29*pow(kdp[i],0.77);
             else {
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*ZHAIL);
             }
         }
    }

    /* fill hail cells below melting layer */
    for (i=0; i<binmeltstart; i++) {
         if (eclass[i]==HAIL) {
             if (kdp[i]>0)
                 rate[i]=29*pow(kdp[i],0.77);
             else {
                 if (z[i]>ZHAIL)
                     z[i]=ZHAIL;
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*z[i]); // hail is an indication of convection
             }
         }
    }
}

float getIntercept(float t)
{
    /* R(A) parameter */
    float coeff[4]={2230., 3100., 4120., 5330.}, sT[4]={0.,10.,20.,30.}, m;
    int i;

    if (t<=0)
        return coeff[0];
    else if (t>=30)
        return coeff[3];
    else {
        i=(int)(t/10);
        m=(coeff[i+1]-coeff[i])/(sT[i+1]-sT[i])*(t-sT[i])+coeff[i];
        return m;
    }
}

void getRayRAwVPR(struct range_parameter evpr,
                  float *z,
                  float *phi,
                  float *kdp,
                  int *eclass,
                  int *rtype,
                  short *melt,
                  float **vpr,
                  float fangle,
                  float alpha,
                  float *rate)
{
    int i, j, k, last, *seg_start, *seg_end, nseg, istart, iend, count, binmeltstart;
    float tI, subI, b=0.62, Fphi, up, down, A, phi_thresh=3., phi_delta;


    /* arrays */
    seg_start=get1dIntArray(MAXIMUM_RANGE_SUPER);
    seg_end=get1dIntArray(MAXIMUM_RANGE_SUPER);

    /* find a bin where the melting layer starts */
    binmeltstart=0;
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (melt[i]==MELTING_FLAG || melt[i]==FREEZING_FLAG) {
             binmeltstart=i;
             break;
         }
    }

    /********** find rain segment that does not contain hail ***************************/
    /* find the start and end of rain segment by phidp */
    count=0;
    for (i=0; i<binmeltstart; i++) {
         if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
             count=1;
             istart=i;
             break;
         }
    }

    nseg=0;
    if (count!=0) {
        for (i=binmeltstart-1; i>=0; i--) {
             if (eclass[i]==RAIN && phi[i]!=INITIAL_VALUE) {
                 iend=i;
                 break;
             }
        }

        for (i=istart; i<=iend; i++) {
             if (i==istart) {
                 if (phi[i]!=INITIAL_VALUE)
                     seg_start[nseg]=i;
                 continue;
             }

             if (eclass[i]==RAIN && eclass[i+1]!=RAIN && phi[i]!=INITIAL_VALUE) {
                 k=0;
                 for (j=i+1; j<iend; j++) {
                      if (eclass[j]==NORAIN || eclass[j]==DONTKNOW || eclass[j]==RAIN) {
                          k++;
                          last=j;
                      }
                      if (eclass[j]==CLUTTER || eclass[j]==HAIL) {
                          if (k==0)
                              seg_end[nseg]=i;
                          else
                              seg_end[nseg]=last;
                          nseg++;
                          i=j;
                          break;
                      }
                 }
                 /*
                 if (j==iend) {
                     seg_end[nseg]=iend;
                     nseg++;
                     break;
                 } */
             }

             if ((eclass[i]==HAIL || eclass[i]==CLUTTER) && eclass[i+1]==RAIN && phi[i+1]!=INITIAL_VALUE)
                 seg_start[nseg]=i+1;
        }

        /* check the last bin */
        if (eclass[iend]!=HAIL) {
            seg_end[nseg]=iend;
            nseg++;
        }
    }
    /*---------------------------------------------------------------------------*/

    /* vpr correction */
    float hcen, htop, hbot, hgrid, gstd, gpdf, cdf, sumvp, deltaz;
    int itop, ibot, imax, q;

    imax=(int)(evpr.hmax/evpr.dh+0.5);
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (z[i]!=NODATA && z[i]>INITIAL_VALUE && rtype[i]==STRATIFORM) {
             htop=computeBeamHeight(fangle+0.25,i*0.25+0.125);
             hcen=computeBeamHeight(fangle,i*0.25+0.125);
             hbot=computeBeamHeight(fangle-0.25,i*0.25+0.125);
             itop=(int)(htop/evpr.dh-0.5);
             ibot=(int)(hbot/evpr.dh+0.5);
             if (ibot>=imax || itop>=imax)
                 continue;

             if (ibot<0)
                 ibot=0;
             gstd=0.85*(hcen-hbot);
             cdf=0.;
             sumvp=0.;
             for (q=ibot; q<=itop; q++) {
                  hgrid=q*evpr.dh+evpr.dh*0.5;
                  gpdf=1/(sqrt(2*M_PI)*gstd)*exp(-(hcen-hgrid)*(hcen-hgrid)/(2*gstd*gstd)); // pdf value
                  cdf+=gpdf;
                  if (vpr[q][0]!=0)
                      sumvp+=gpdf*vpr[q][0];
             }

             if (cdf==0 || sumvp==0)
                 continue;
             else
                 deltaz=10*log10(sumvp/cdf);

             z[i]-=deltaz;
         }
    }

    /********** start computation for path-integrated attenuation (PIA) ************************/
    for (i=0; i<nseg; i++) {
         phi_delta=phi[seg_end[i]]-phi[seg_start[i]];
         if (phi_delta>=phi_thresh) { // then, use R(A)
             /* calculate total integral for the segment */
             tI=0.;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]!=INITIAL_VALUE)
                      tI=0.46*b*pow(pow(10.,z[j]/10.),b)*((float)BIN_SUPER/1000)+tI;
             }

             Fphi=exp(0.23*alpha*b*(phi_delta))-1;
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  subI=0.;
                  for (k=j; k<=seg_end[i]; k++) {
                       if (z[k]!=NODATA && z[k]!=INITIAL_VALUE)
                           subI=0.46*b*pow(pow(10.,z[k]/10.),b)*((float)BIN_SUPER/1000)+subI;
                  }

                  /* specific attenuation A(r) */
                  if (z[j]!=NODATA && z[j]!=INITIAL_VALUE) {
                      up=pow(pow(10.,z[j]/10.),b)*Fphi;
                      down=tI+Fphi*subI;
                      A=up/down;
                      if (eclass[j]==RAIN)
                          rate[j]=4120*pow(A,1.03);
                  }
             }
         }
         else { // use R(Z)
             for (j=seg_start[i]; j<=seg_end[i]; j++) {
                  if (z[j]!=NODATA && z[j]!=INITIAL_VALUE) {
                      if (eclass[j]==CLUTTER)
                          rate[j]=0.;
                      else { // there is no hail in rain segments
                          if (rtype[j]==CONVECTIVE)
                              rate[j]=1.7007*pow(10,-2)*pow(10.,0.0714*z[j]); // NEXRAD
                          else
                              rate[j]=3.6463*pow(10,-2)*pow(10.,0.0625*z[j]); // M-P

                      }
                  }
             }
         }
    }
    free(seg_start);
    free(seg_end);
    /*--------------------------------------------------------------------------------------------*/

    /* above melting layer */
    for (i=binmeltstart; i<MAXIMUM_RANGE_SUPER; i++) {
         if (z[i]!=NODATA && z[i]>INITIAL_VALUE) {
             if (rtype[i]==CONVECTIVE)
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*z[i]); // NEXRAD
             else
                 rate[i]=3.6463*pow(10,-2)*pow(10.,0.0625*z[i]); // M-P
         }
    }

    /* fill hail cells */
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (eclass[i]==HAIL) {
             if (kdp[i]>0)
                 rate[i]=29*pow(kdp[i],0.77);
             else {
                 if (z[i]>ZHAIL)
                     z[i]=ZHAIL;
                 rate[i]=1.7007*pow(10,-2)*pow(10.,0.0714*z[i]); // NEXRAD (hail: convection)
             }
         }
    }
}

void estimateRZdefault(struct volume_info header,
                       short ***z,
                       int ***echo,
                       float ***rr,
                       int angle)
{
    /* estimate rain rate based on NEXRAD Z-R */
    int echo_id, ang, azi, rng, nazi;
    float zval;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  echo_id=echo[azi][rng][ang];
                  if (echo_id!=CLUTTER) {
                      if (zval>=ZHAIL) // hail
                          zval=ZHAIL;
                      rr[azi][rng][ang]=1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                  }
                  else
                      rr[azi][rng][ang]=0.;
              }
         }
    }
}

void estimateRZC(struct volume_info header,
                 short ***z,
                 int ***kdp,
                 int ***echo,
                 short ***rtype,
                 float ***rr,
                 int angle)
{
    /* estimate rain rate based on NEXRAD Z-R */
    int echo_id, ang, azi, rng, nazi;
    float zval, kdpval, C1;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  kdpval=kdp[azi][rng][ang]/DP_SCALEING;
                  echo_id=echo[azi][rng][ang];
                  if (echo_id!=DONTKNOW && echo_id!=CLUTTER && echo_id!=NORAIN && echo_id!=HAIL) {
                      if (echo_id==RAIN || rtype[azi][rng][ang]==CONVECTIVE)
                          C1=1.;
                      else {
                          if (echo_id==DRYSNOW)
                              C1=1.;
                          else if (echo_id==WETSNOW)
                              C1=0.6;
                          else if (echo_id==SNOWICE || echo_id==ICE)
                              C1=2.8;
                      }
                      rr[azi][rng][ang]=C1*1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                  }
                  else if (echo_id==HAIL) {
                      zval=ZHAIL;
                      rr[azi][rng][ang]=1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                  }
              }
         }
    }
}

void estimateRZCS(struct volume_info header,
                  short ***z,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle)
{
    /* estimate rain rate based on Z-Rs using rain type classification */
    int echo_id, ang, azi, rng, nazi;
    float zval, kdpval, C1;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  kdpval=kdp[azi][rng][ang]/DP_SCALEING;
                  echo_id=echo[azi][rng][ang];
                  if (echo_id!=DONTKNOW && echo_id!=CLUTTER && echo_id!=NORAIN && echo_id!=HAIL) {
                      if (echo_id==RAIN) {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                          else
                              rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P
                      }
                      else {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              C1=1.;
                          else {
                              if (echo_id==DRYSNOW)
                                  C1=1.;
                              else if (echo_id==WETSNOW)
                                  C1=0.6;
                              else if (echo_id==SNOWICE || echo_id==ICE)
                                  C1=2.8;
                          }
                          rr[azi][rng][ang]=C1*1.7007*pow(10.,-2.)*pow(10.,0.0714*zval);
                      }
                  }
                  else if (echo_id==HAIL) {
                      zval=ZHAIL;
                      rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval);
                  }
              }
         }
    }
}

void estimateRZDR(struct volume_info header,
                  short ***z,
                  int ***zdr,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle,
                  int scale)
{
    /* estimate rain rate based on R(Z,Zdr) */
    int echo_id, ang, azi, rng, nazi;
    float zval, zdrval, kdpval, C1;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  zdrval=(float)zdr[azi][rng][ang]/scale;
                  kdpval=kdp[azi][rng][ang]/DP_SCALEING;
                  echo_id=echo[azi][rng][ang];
                  if (echo_id!=DONTKNOW && echo_id!=CLUTTER && echo_id!=NORAIN && echo_id!=HAIL) {
                      if (echo_id==RAIN)
                          rr[azi][rng][ang]=6.7*pow(10.,-3.)*pow(10.,0.0927*zval)*pow(10.,-0.343*zdrval);
                      else {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              C1=1.;
                          else {
                              if (echo_id==DRYSNOW)
                                  C1=1.;
                              else if (echo_id==WETSNOW)
                                  C1=0.6;
                              else if (echo_id==SNOWICE || echo_id==ICE)
                                  C1=2.8;
                          }
                          rr[azi][rng][ang]=C1*1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                      }
                  }
                  else if (echo_id==HAIL) {
                      if (kdpval>0)
                          rr[azi][rng][ang]=29*pow(kdpval,0.77);
                      else {
                          zval=ZHAIL;
                          rr[azi][rng][ang]=1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                      }
                  }
              }
         }
    }
}

void estimateRKDP(struct volume_info header,
                  short ***z,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle,
                  int scale)
{
    /* estimate rain rate based on R(Z,Zdr) */
    int echo_id, ang, azi, rng, nazi;
    float zval, kdpval, C1;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  kdpval=(float)kdp[azi][rng][ang]/scale;
                  echo_id=echo[azi][rng][ang];
                  if (echo_id!=DONTKNOW && echo_id!=CLUTTER && echo_id!=NORAIN && echo_id!=HAIL) {
                      if (echo_id==RAIN) {
                          if (kdpval>0)
                              rr[azi][rng][ang]=40.5*pow(kdpval,0.85);
                          else
                              rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P (light rain)
                      }
                      else {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              C1=1.;
                          else {
                              if (echo_id==DRYSNOW)
                                  C1=1.;
                              else if (echo_id==WETSNOW)
                                  C1=0.6;
                              else if (echo_id==SNOWICE || echo_id==ICE)
                                  C1=2.8;
                          }
                          rr[azi][rng][ang]=C1*1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                      }
                  }
                  else if (echo_id==HAIL) {
                      if (kdpval>0)
                          rr[azi][rng][ang]=29*pow(kdpval,0.77);
                      else {
                          zval=ZHAIL;
                          rr[azi][rng][ang]=1.7007*pow(10,-2)*pow(10.,0.0714*zval);
                      }
                  }
              }
         }
    }
}

void estimateRZCwVPR(struct volume_info header,
                     struct range_parameter evpr,
                     short ***z,
                     int ***echo,
                     short ***rtype,
                     float **vpr,
                     float *eangle,
                     int angle,
                     float ***rr)
{
    /* estimate rain rate based on NEXRAD Z-R with VPR correction */
    int echo_id, ang, azi, rng, nazi;
    float zval;

    /* vpr parameters */
    float hcen, htop, hbot, hgrid, gstd, gpdf, cdf, sumvp, deltaz;
    int itop, ibot, imax, q;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    imax=(int)(evpr.hmax/evpr.dh+0.5);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  echo_id=echo[azi][rng][ang];
                  if (rtype[azi][rng][ang]==STRATIFORM) { /* vpr correction */
                      htop=computeBeamHeight(eangle[ang]+0.25,rng*0.25+0.125);
                      hcen=computeBeamHeight(eangle[ang],rng*0.25+0.125);
                      hbot=computeBeamHeight(eangle[ang]-0.25,rng*0.25+0.125);
                      itop=(int)(htop/evpr.dh-0.5);
                      ibot=(int)(hbot/evpr.dh+0.5);
                      if (ibot>=imax || itop>=imax)
                          continue;
                      if (ibot<0)
                          ibot=0;
                      gstd=0.85*(hcen-hbot);
                      cdf=0.;
                      sumvp=0.;
                      for (q=ibot; q<=itop; q++) {
                           hgrid=q*evpr.dh+evpr.dh*0.5;
                           gpdf=1/(sqrt(2*M_PI)*gstd)*exp(-(hcen-hgrid)*(hcen-hgrid)/(2*gstd*gstd)); // pdf value
                           cdf+=gpdf;
                           if (vpr[q][0]!=0)
                               sumvp+=gpdf*vpr[q][0];
                      }

                      if (cdf==0 || sumvp==0)
                          continue;
                      else
                          deltaz=10*log10(sumvp/cdf);

                      zval-=deltaz;
                  }
                  if (echo_id==RAIN || echo_id==DRYSNOW || echo_id==WETSNOW || echo_id==SNOWICE || echo_id==ICE)
                      rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                  else if (echo_id==HAIL) {
                      if (zval>=ZHAIL)
                          zval=ZHAIL;
                      rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                  }
              }
         }
    }
}

void estimateRZCSwVPR(struct volume_info header,
                      struct range_parameter evpr,
                      short ***z,
                      int ***echo,
                      short ***rtype,
                      float **vpr,
                      float *eangle,
                      int angle,
                      float ***rr)
{
    /* estimate rain rate based on Z-Rs using rain type classification with VPR correction*/
    int echo_id, ang, azi, rng, nazi;
    float zval;

    /* vpr parameters */
    float hcen, htop, hbot, hgrid, gstd, gpdf, cdf, sumvp, deltaz;
    int itop, ibot, imax, q;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    imax=(int)(evpr.hmax/evpr.dh+0.5);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  echo_id=echo[azi][rng][ang];
                  if (rtype[azi][rng][ang]==STRATIFORM) { /* vpr correction */
                      htop=computeBeamHeight(eangle[ang]+0.25,rng*0.25+0.125);
                      hcen=computeBeamHeight(eangle[ang],rng*0.25+0.125);
                      hbot=computeBeamHeight(eangle[ang]-0.25,rng*0.25+0.125);
                      itop=(int)(htop/evpr.dh-0.5);
                      ibot=(int)(hbot/evpr.dh+0.5);
                      if (ibot>=imax || itop>=imax)
                          continue;
                      if (ibot<0)
                          ibot=0;
                      gstd=0.85*(hcen-hbot);
                      cdf=0.;
                      sumvp=0.;
                      for (q=ibot; q<=itop; q++) {
                           hgrid=q*evpr.dh+evpr.dh*0.5;
                           gpdf=1/(sqrt(2*M_PI)*gstd)*exp(-(hcen-hgrid)*(hcen-hgrid)/(2*gstd*gstd)); // pdf value
                           cdf+=gpdf;
                           if (vpr[q][0]!=0)
                               sumvp+=gpdf*vpr[q][0];
                      }

                      if (cdf==0 || sumvp==0)
                          continue;
                      else
                          deltaz=10*log10(sumvp/cdf);

                      zval-=deltaz;
                  }
                  if (echo_id==RAIN || echo_id==DRYSNOW || echo_id==WETSNOW || echo_id==SNOWICE || echo_id==ICE) {
                      if (rtype[azi][rng][ang]==CONVECTIVE)
                          rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                      else
                          rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P
                   }
                   else if (echo_id==HAIL) {
                      if (zval>=ZHAIL)
                          zval=ZHAIL;
                      rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                  }
              }
         }
    }
}

void estimateRZDRwVPR(struct volume_info header,
                      struct range_parameter evpr,
                      short ***z,
                      int ***zdr,
                      int ***echo,
                      short ***rtype,
                      short ***melting,
                      float **vpr,
                      float *eangle,
                      int angle,
                      int scale,
                      float ***rr)
{
    /* estimate rain rate based on R(Z,Zdr) with VPR correction */
    int echo_id, ang, azi, rng, nazi;
    float zval, zdrval;

    /* vpr parameters */
    float hcen, htop, hbot, hgrid, gstd, gpdf, cdf, sumvp, deltaz;
    int itop, ibot, imax, q;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    imax=(int)(evpr.hmax/evpr.dh+0.5);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  zdrval=(float)zdr[azi][rng][ang]/scale;
                  echo_id=echo[azi][rng][ang];
                  if (rtype[azi][rng][ang]==STRATIFORM) { /* vpr correction */
                      htop=computeBeamHeight(eangle[ang]+0.25,rng*0.25+0.125);
                      hcen=computeBeamHeight(eangle[ang],rng*0.25+0.125);
                      hbot=computeBeamHeight(eangle[ang]-0.25,rng*0.25+0.125);
                      itop=(int)(htop/evpr.dh-0.5);
                      ibot=(int)(hbot/evpr.dh+0.5);
                      if (ibot>=imax || itop>=imax)
                          continue;
                      if (ibot<0)
                          ibot=0;
                      gstd=0.85*(hcen-hbot);
                      cdf=0.;
                      sumvp=0.;
                      for (q=ibot; q<=itop; q++) {
                           hgrid=q*evpr.dh+evpr.dh*0.5;
                           gpdf=1/(sqrt(2*M_PI)*gstd)*exp(-(hcen-hgrid)*(hcen-hgrid)/(2*gstd*gstd)); // pdf value
                           cdf+=gpdf;
                           if (vpr[q][0]!=0)
                               sumvp+=gpdf*vpr[q][0];
                      }

                      if (cdf==0 || sumvp==0)
                          continue;
                      else
                          deltaz=10*log10(sumvp/cdf);

                      zval-=deltaz;
                  }

                  if (melting[azi][rng][ang]!=MELTING_FLAG && melting[azi][rng][ang]!=FREEZING_FLAG) { /* below melting layer */
                      if (echo_id==RAIN)
                          rr[azi][rng][ang]=6.7*pow(10.,-3.)*pow(10.,0.0927*zval)*pow(10.,-0.343*zdrval);
                      else if (echo_id==HAIL) {
                          if (zval>=ZHAIL)
                              zval=ZHAIL;
                          rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval);
                      }
                  }
                  else { /* above melting layer -> use Z(R) */
                      if (echo_id==RAIN || echo_id==DRYSNOW || echo_id==WETSNOW || echo_id==SNOWICE || echo_id==ICE) {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                          else
                              rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P
                      }
                      else if (echo_id==HAIL) {
                          if (zval>=ZHAIL)
                              zval=ZHAIL;
                          rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval);
                      }
                  }
              }
         }
    }
}

void estimateRKDPwVPR(struct volume_info header,
                      struct range_parameter evpr,
                      short ***z,
                      int ***kdp,
                      int ***echo,
                      short ***rtype,
                      short ***melting,
                      float **vpr,
                      float *eangle,
                      int angle,
                      int scale,
                      float ***rr)
{
    /* estimate rain rate based on R(Kdp) with VPR correction */
    int echo_id, ang, azi, rng, nazi;
    float zval, kdpval;

    /* vpr parameters */
    float hcen, htop, hbot, hgrid, gstd, gpdf, cdf, sumvp, deltaz;
    int itop, ibot, imax, q;

    ang=angle-1;
    //if (ang<nsuper)
    if (header.typeangle[ang]==SUPER)
        nazi=MAXIMUM_AZIMUTH_SUPER;
    else
        nazi=MAXIMUM_AZIMUTH;

    imax=(int)(evpr.hmax/evpr.dh+0.5);

    for (azi=0; azi<nazi; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][ang]!=NODATA && z[azi][rng][ang]>INITIAL_VALUE) {
                  zval=int_dBZ2float_dBZ(z[azi][rng][ang]);
                  kdpval=(float)kdp[azi][rng][ang]/scale;
                  echo_id=echo[azi][rng][ang];
                  if (melting[azi][rng][ang]!=MELTING_FLAG  && melting[azi][rng][ang]!=FREEZING_FLAG) { /* below melting layer */
                      if (echo_id==RAIN) {
                          if (kdpval>0)
                              rr[azi][rng][ang]=40.5*pow(kdpval,0.85);
                          else
                              rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P (light rain)
                      }
                      else if (echo_id==HAIL) {
                          if (kdpval>0)
                              rr[azi][rng][ang]=29*pow(kdpval,0.77);
                          else {
                              if (zval>=ZHAIL)
                                  zval=ZHAIL;
                              rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval);
                          }
                      }
                  }
                  else { /* above melting layer -> use Z(R) */
                      if (rtype[azi][rng][ang]==STRATIFORM) { /* vpr correction */
                          htop=computeBeamHeight(eangle[ang]+0.25,rng*0.25+0.125);
                          hcen=computeBeamHeight(eangle[ang],rng*0.25+0.125);
                          hbot=computeBeamHeight(eangle[ang]-0.25,rng*0.25+0.125);
                          itop=(int)(htop/evpr.dh-0.5);
                          ibot=(int)(hbot/evpr.dh+0.5);
                          if (ibot>=imax || itop>=imax)
                              continue;
                          if (ibot<0)
                              ibot=0;
                          gstd=0.85*(hcen-hbot);
                          cdf=0.;
                          sumvp=0.;
                          for (q=ibot; q<=itop; q++) {
                               hgrid=q*evpr.dh+evpr.dh*0.5;
                               gpdf=1/(sqrt(2*M_PI)*gstd)*exp(-(hcen-hgrid)*(hcen-hgrid)/(2*gstd*gstd)); // pdf value
                               cdf+=gpdf;
                               if (vpr[q][0]!=0)
                                   sumvp+=gpdf*vpr[q][0];
                          }

                          if (cdf==0 || sumvp==0)
                              continue;
                          else
                              deltaz=10*log10(sumvp/cdf);

                          zval-=deltaz;
                      }

                      if (echo_id==RAIN || echo_id==DRYSNOW || echo_id==WETSNOW || echo_id==SNOWICE || echo_id==ICE) {
                          if (rtype[azi][rng][ang]==CONVECTIVE)
                              rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                          else
                              rr[azi][rng][ang]=3.6463*pow(10.,-2.)*pow(10.,0.0625*zval); // M-P
                      }
                      else if (echo_id==HAIL) {
                          if (zval>=ZHAIL)
                              zval=ZHAIL;
                          rr[azi][rng][ang]=1.7007*pow(10.,-2.)*pow(10.,0.0714*zval); // NEXRAD
                      }
                  }
              }
         }
    }
}

void get2DFixedGrid(float *razimuth,
                    float BWidth,
                    float **data)
{
    int azm, rng, radial, intAzC, intAzB, intAzA, status;
    float azAngle, weightA, weightB, weightC, fract_weighthresh, vdata,
          **sumWts, **sumWtdZ;

    /* get arrays */
    sumWts=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    sumWtdZ=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    // initialize the azimuthal weighting and fixed gird
    for (azm=0; azm<MAXIMUM_AZIMUTH_SUPER; azm++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              sumWtdZ[azm][rng] = INITIAL_VALUE;
              sumWts[azm][rng] = INITIAL_VALUE;
         }
	}

    /* Process each incoming radial of base data.  This loop defines the
    information needed to remap the reflectivity data to the whole degree
    Hybrid Scan, tests to see if the incoming radial is blocked or
    contaminated by AP/clutter, and implements the user selectable exclusion
    zones (if any have been defined). */
    for (radial=0; radial<MAXIMUM_AZIMUTH_SUPER; radial++){
         azAngle=razimuth[radial]*2;
         if (azAngle<0)
             continue;

    /* Determine the whole degree azimuth angles that overlap with this base data
    azimuth.  The whole degree angles will be used to store weighted
    reflectivity information into the temporary tilt scan.  Note, this
    technique is only valid for BEAM_WDTH <= 1.0 degrees. */

        intAzC=(int)(azAngle+BWidth);

        if (intAzC>MAXIMUM_AZIMUTH_SUPER-1) intAzC=0;
        if ((intAzB=intAzC-1)<0) intAzB=MAXIMUM_AZIMUTH_SUPER+intAzB;
        if ((intAzA=intAzB-1)<0) intAzA=MAXIMUM_AZIMUTH_SUPER+intAzA;

    /* Calculate the relative contribution (weight) of the beam for each of the
    overlapping whole degree azimuths. */

        if((weightC=azAngle+BWidth)>=MAXIMUM_AZIMUTH_SUPER)
            weightC=weightC-MAXIMUM_AZIMUTH_SUPER;

        weightC=weightC-intAzC;

        if (weightC>=2*BWidth-1) {
            weightB=2*BWidth-weightC;
            weightA=0;
        }
        else {
            weightB=1;
            weightA=2*BWidth-weightB;
        }

        weightA=weightA/2;
        weightB=weightB/2;
        weightC=weightC/2;

    /* For each 1 km range bin, ensure that the Hybrid Scan uses the lowest
    elevation reflectivity data that is uncontaminated and unblocked. */
        for (rng=0;rng<MAXIMUM_RANGE_SUPER;rng++){
             vdata=data[radial][rng];
             if (vdata != NODATA) { // && v_DBZ != INITIAL_VALUE){
                sumWts[intAzA][rng]+=weightA;
                sumWts[intAzB][rng]+=weightB;
                sumWts[intAzC][rng]+=weightC;
                sumWtdZ[intAzA][rng]+=(weightA*vdata);
                sumWtdZ[intAzB][rng]+=(weightB*vdata);
                sumWtdZ[intAzC][rng]+=(weightC*vdata);
             }
        }
    }
    fract_weighthresh=(float)0.50;
    for (azm=0;azm<MAXIMUM_AZIMUTH_SUPER;azm++){
         for (rng=0;rng<MAXIMUM_RANGE_SUPER;rng++){
              /* Compute the average power */
              if (sumWts[azm][rng] > fract_weighthresh) {
                  data[azm][rng]=sumWtdZ[azm][rng]/sumWts[azm][rng];
              }
         }
    }

//    status=free2dArray(sumWts,MAXIMUM_AZIMUTH_SUPER);
//    status=free2dArray(sumWtdZ,MAXIMUM_AZIMUTH_SUPER);
}
