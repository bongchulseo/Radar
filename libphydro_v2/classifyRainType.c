/*************************************************************************
--> ClassifyRainType.c

    Use VIL (vertically integrated liquid) and identify rain types (convective/stratiform)

    # Functions
      - computeVIL
        1. reference Z construction for the mitigation of bright band effect
        2.




    Developed by Bongchul Seo and Witold F. Krajewski
    The University of Iowa

    Created : Mar       2018

*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

void computeVIL(struct volume_info header,
                short ***z,
                int ***rho,
                int ***phi,
                float ***T_top,
                float ***T_bottom,
                float **azimuth,
                short **lookup,
                float **VIL,
                short **wf)
{
    FILE *fvil;
    int ang, azi, cazi, c2azi, rng, maxangle, dw=1, vcount,
        rho_bad=0.85*DP_SCALEING, rho_perfect=1.00*DP_SCALEING;
    float dh, ch, htotal, uVIL, h_high=20., bb_low_degree=5., bb_high_degree=-5., bb_over_degree=-30.,
          **z_rf, *z_ver, *t_ver, zdiff, bb_enhance=3., zval, t_cent;
    int i, imax, iclose, wf_count, rainmode;
    float tmax, tmin, tclose, z0, z10, z5, slope, zgrowing, raincover;
    short wf_flag, ok_flag, ok_flag2;

    /* parameters (thresholds) */
    float VIL_convective=6.5, VIL_low=2.;
    float Z_thresh_dBZ=10., Z_high_dBZ=50.;
    float raincover_thresh=0.5;

    maxangle=header.nangle+header.nangle_super+header.nangle_legacy2;

    /* limit the highest elevation for VIL calculation */
    if (maxangle>=9)
        maxangle=9;

    rainmode=NO;
    if (header.nangle_super>=3)
        rainmode=YES;

    /* Z reference array for bright band mitigation */
    z_rf=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    initialize2dFloatArray(z_rf,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,NODATA);
    z_ver=get1dFloatArray(maxangle);
    t_ver=get1dFloatArray(maxangle);

    /* reference value construction */
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              vcount=0;
              for (ang=0; ang<maxangle; ang++) {
                   /* corresponding azimuth */
                   if (ang==0)
                       cazi=azi;
                   else
                       cazi=lookup[azi][ang];

                   if (z[cazi][rng][ang]!=NODATA && z[cazi][rng][ang]!=INITIAL_VALUE) {
                       if (int_dBZ2float_dBZ(z[cazi][rng][ang])>Z_thresh_dBZ) {
                           //&& rho[cazi][rng][ang]>rho_bad && rho[cazi][rng][ang]<=rho_perfect) {
                           /* we do not include rho values when considering reference Z values
                              because rho has a range limit (bright band at higer elevations are also affected)
                              at certain elevation angles, which results in range rings in the map of VIL. */
                           if (ang==0 && T_bottom[cazi][rng][ang]<=-5.) {
                               z_rf[azi][rng]=int_dBZ2float_dBZ(z[cazi][rng][ang]);
                               break;
                           }

                           if (T_top[cazi][rng][ang]>bb_low_degree) { // bb_low_degree: 5 degrees
                               z_rf[azi][rng]=int_dBZ2float_dBZ(z[cazi][rng][ang]);
                               break;
                           }
                           else if (T_bottom[cazi][rng][ang]>=bb_over_degree && T_top[cazi][rng][ang]<=bb_low_degree) { // bb_over_degree: -30 degrees
                               z_ver[vcount]=int_dBZ2float_dBZ(z[cazi][rng][ang]);
                               t_ver[vcount]=0.5*(T_top[cazi][rng][ang]+T_bottom[cazi][rng][ang]);
                               vcount++;
                           }
                       }
                   }
              }

              /* look at vertical structure */
              if (z_rf[azi][rng]==NODATA && vcount>1) {
                  tmax=-100.;
                  tmin=100.;
                  for (i=0; i<vcount; i++) {
                       if (t_ver[i]>=tmax) {
                           tmax=t_ver[i];
                           imax=i;
                       }
                       if (t_ver[i]<=tmin)
                           tmin=t_ver[i];
                  }
                  if (tmax>0 && tmin<-10) {
                      z0=getinterpZ(t_ver,z_ver,vcount,0.);
                      z10=getinterpZ(t_ver,z_ver,vcount,-10.);
                      slope=(z0-z10)/10;
//if (azi==320 && rng==135)
//    printf("tmax: %f, tmin: %f, z0: %f, z10: %f, slope: %f\n",tmax,tmin,z0,z10,slope);

                      if (slope>=1.5) { /* stratiform (dB/degree) */
                          z5=getinterpZ(t_ver,z_ver,vcount,-5.); // interpolated reflectivity at -5 degrees
                          z_rf[azi][rng]=z5;
                      }
                      else
                          z_rf[azi][rng]=z_ver[0];
                  }
                  else {
                      z5=getinterpZ(t_ver,z_ver,vcount,-5.);
//if (azi==320 && rng==135)
//    printf("tmax: %f, tmin: %f, z5: %f\n",tmax,tmin,z5);

                      if (z5!=0)
                          z_rf[azi][rng]=z5;
                      else {  // -5 degrees are not in the temperature range
                          if (t_ver[0]<-5)
                              z_rf[azi][rng]=z_ver[0];
                          else {
                              tclose=100.;
                              for (i=0; i<vcount; i++) {
                                   if (fabs(t_ver[i]-(-5))<=tclose) {
                                       tclose=fabs(t_ver[i]-(-5));
                                       iclose=i;
                                   }
                              }
                              //if (tclose<2.) // if the difference is smaller than 2 degrees
                                  z_rf[azi][rng]=z_ver[iclose];
                              //else
                              //    z_rf[azi][rng]=z_ver[0]-20.;

                              //if (z_rf[azi][rng]<0)
                              //    z_rf[azi][rng]=Z_thresh_dBZ;
                          }
                      }
                  }
              }
         }
    }

    /* VIL calculation */
    wf_flag=FALSE;
    wf_count=0;
    initialize2dShortArray(wf,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,NOTHING);
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              htotal=0.;
              for (ang=0; ang<maxangle; ang++) {
                   ch=computeBeamHeight(header.angle[ang],rng*0.25+0.125);
                   if (ch>h_high)
                       break;

                   /* corresponding azimuth */
                   if (ang==0)
                       cazi=azi;
                   else
                       cazi=lookup[azi][ang];

                   /* vertical width */
                   if (ang==0)
                       dh=computeBeamHeight(0.5*(header.angle[ang]+header.angle[ang+1]),rng*0.25+0.125);
                   else if (ang>0 && ang<maxangle-1) {
                       dh=computeBeamHeight(0.5*(header.angle[ang+1]+header.angle[ang]),rng*0.25+0.125)
                          -computeBeamHeight(0.5*(header.angle[ang]+header.angle[ang-1]),rng*0.25+0.125);
                   }
                   else if(ang==maxangle-1) {
                       dh=computeBeamHeight(header.angle[ang]+header.beam_width/2,rng*0.25+0.125)
                          -computeBeamHeight(0.5*(header.angle[ang]+header.angle[ang-1]),rng*0.25+0.125);
                       /* integration up to 20 km */
                       dh=dh+h_high-computeBeamHeight(header.angle[ang]+header.beam_width/2,rng*0.25+0.125);
                   }
                   dh=dh*1000;
                   if (dh<=0)
                       continue;
                   htotal+=dh;

                   if (z[cazi][rng][ang]!=NODATA && z[cazi][rng][ang]!=INITIAL_VALUE) {
                       if (rho[cazi][rng][ang]>rho_bad && rho[cazi][rng][ang]<=rho_perfect) {
                           zval=int_dBZ2float_dBZ(z[cazi][rng][ang]);
                           if (ang==0 && zval<=Z_thresh_dBZ)
                               break;
                           /* bright band effect adjustment */
                           if (z_rf[azi][rng]!=NODATA) {
                               if (T_top[cazi][rng][ang]<bb_low_degree && T_bottom[cazi][rng][ang]>bb_high_degree) {
                                   zdiff=zval-z_rf[azi][rng];
                                   if (zdiff>bb_enhance) {
                                       //if (ang==maxangle-1)
                                       //    zval=10.; // 35 dBZ at 10 km yields VIL > 6.5
                                       //else
                                           zval=z_rf[azi][rng];
                                   }
                               }
                           }
                           if (ang==maxangle-1 && T_bottom[cazi][rng][ang]<20.) {
                               if (zval>30.) {
                                   for (i=0; i<header.nangle_super; i++) {
                                        /* corresponding azimuth */
                                        if (i==0)
                                            c2azi=azi;
                                        else
                                            c2azi=lookup[azi][i];

                                        if (T_bottom[c2azi][rng][i]>5. && int_dBZ2float_dBZ(z[c2azi][rng][i])>=zval)
                                            break;
                                        else {
                                            if (i==header.nangle_super-1) {
                                                zval=0.;
                                            }
                                            continue;
                                        }

                                   }
                               }
                           }
                           uVIL=(float)(3.44*0.001*pow(pow(10.,(double)zval/10),(double)4/7));
                           VIL[azi][rng]+=uVIL*dh/1000; // unit: kg/m^-2
//if (azi==575 && rng==49)
//    printf("ang: %d, azi: %f, z_ref: %f, z: %f, rho: %f, zval: %f, h: %f, VIL: %f, t_top: %f, t_bot: %f\n",ang,azimuth[cazi][ang],z_rf[azi][rng],int_dBZ2float_dBZ(z[cazi][rng][ang]),
//                                            (float)rho[cazi][rng][ang]/DP_SCALEING,zval,ch,VIL[azi][rng],T_top[cazi][rng][ang],T_bottom[cazi][rng][ang]);

                       }
                   }
              }

              /* quality check for convective cells */

              if (VIL[azi][rng]>VIL_convective) {
                  ok_flag=checkVeticalContinuity(z,rho,lookup,maxangle,header.nangle_super,azi,rng,
                                                 z_rf[azi][rng],T_top,T_bottom,Z_high_dBZ-15.,Z_thresh_dBZ);
                  if (ok_flag==FALSE) {
                      wf[azi][rng]=WINDTURBINE;
                      wf_count++;
                  }
              }

              else if (z[azi][rng][0]>=float_dBZ2int_dBZ(Z_high_dBZ) && VIL[azi][rng]<VIL_low) {
                  ok_flag=checkVeticalContinuity(z,rho,lookup,maxangle,header.nangle_super,azi,rng,
                                                 z_rf[azi][rng],T_top,T_bottom,Z_high_dBZ-15.,Z_thresh_dBZ);
                  if (ok_flag==FALSE) {
                      wf[azi][rng]=WINDTURBINE;
                      wf_count++;
                  }
              }

              if (z[azi][rng][0]>=float_dBZ2int_dBZ((float)25.) && VIL[azi][rng]<1.0) {
                  if (rainmode==YES) { // rain mode
                      raincover=getCoverage(z,rho,azi,rng,0,MAXIMUM_AZIMUTH_SUPER);
                      if (raincover>=raincover_thresh) {
                          ok_flag=checkVeticalContinuity(z,rho,lookup,maxangle,header.nangle_super,azi,rng,
                                                         z_rf[azi][rng],T_top,T_bottom,35.,Z_thresh_dBZ);
                      }
                      else {
                          ok_flag=checkVeticalContinuity(z,rho,lookup,maxangle,header.nangle_super,azi,rng,
                                                         z_rf[azi][rng],T_top,T_bottom,25.,Z_thresh_dBZ);
                      }

                      if (ok_flag==FALSE) {
                          wf[azi][rng]=WINDTURBINE;
                          wf_count++;

                      }
                  }
                  else { // no rain mode
                      ok_flag=checkVeticalContinuity(z,rho,lookup,maxangle,header.nangle_super,azi,rng,
                                                     z_rf[azi][rng],T_top,T_bottom,25.,Z_thresh_dBZ);
                      if (ok_flag==FALSE) {
                          wf[azi][rng]=WINDTURBINE;
                          wf_count++;

                      }
                  }
              }

         }
    }
    averageProduct(VIL,dw);


    /* seeded region growing for wind turbine cluster */
//    identifyWF(z,rho,z_rf,wf,40.,15.);

    free2dFloatArray(z_rf,MAXIMUM_AZIMUTH_SUPER);
    free(z_ver);
    free(t_ver);
//printf("wf_flag: %d\n",wf_flag);
//printf("wf_count: %d\n",wf_count);

/*
    fvil=fopen("VIL.out","w");
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         fprintf(fvil,"%.2f ",azimuth[azi][0]);
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              fprintf(fvil,"%f ",VIL[azi][rng]);
         }
         fprintf(fvil,"\n");
    }
    fclose(fvil);
*/
/*
    fvil=fopen("Z.out","w");
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         fprintf(fvil,"%.2f ",azimuth[azi][0]);
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (z[azi][rng][0]>0)
              //if (rho[azi][rng][0]>rho_bad && rho[azi][rng][0]<=rho_perfect)
                  fprintf(fvil,"%6.2f ",int_dBZ2float_dBZ(z[azi][rng][0]));
              else
                  fprintf(fvil,"%d ",INITIAL_VALUE);
         }
         fprintf(fvil,"\n");
    }
    fclose(fvil);
    exit(0);
*/
//    fvil=fopen("RCoverage.out","w");
//    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
//         fprintf(fvil,"%.2f ",azimuth[azi][0]);
//         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
//              fprintf(fvil,"%.2f ",getCoverage(z,rho,azi,rng,0,MAXIMUM_AZIMUTH_SUPER));
//         }
//         fprintf(fvil,"\n");
//    }
//    fclose(fvil);

//    exit(0);

}

int checkHorizontalContiuity(short ***z,
                             int ***rho,
                             int ***phi,
                             int azi,
                             int rng)
{
    int flag, ang=0, avg_rho, std_phidp, rho_thresh_high=0.98*DP_SCALEING, rho_thresh_low=0.90*DP_SCALEING,
        rho_bad=0.85*DP_SCALEING, rho_perfect=1.00*DP_SCALEING, dbZ_higher=35, phidp_thresh=20*DP_SCALEING;;

    std_phidp=getStd(phi,azi,rng,ang,MAXIMUM_AZIMUTH_SUPER);
    avg_rho=getMean(rho,azi,rng,ang,MAXIMUM_AZIMUTH_SUPER);

    flag=TRUE;
    if (int_dBZ2float_dBZ(z[azi][rng][ang])>40.) {
        if (avg_rho<rho_thresh_low && std_phidp>phidp_thresh)
            flag=FALSE;
    }
    else {
        if (avg_rho<rho_thresh_low && std_phidp>phidp_thresh)
            flag=FALSE;
    }
    return flag;
}

int checkVeticalContinuity(short ***z,
                           int ***rho,
                           short **lookup,
                           int nangle,
                           int nsplitcut,
                           int azi,
                           int rng,
                           float z_rf,
                           float ***T_top,
                           float ***T_bottom,
                           float z_value,
                           float z_thresh)
{
    int ang, cazi, flag;
    float bb_low_degree=5., bb_high_degree=-5., bb_enhance=5., zval, zdiff, *z_ver,
          rho_bad=0.85*DP_SCALEING, rho_perfect=1.00*DP_SCALEING;

    z_ver=get1dFloatArray(nangle);

    flag=TRUE;
    for (ang=0; ang<nangle; ang++) {
         if (ang==0)
             cazi=azi;
         else
             cazi=lookup[azi][ang];

         if (z[cazi][rng][ang]!=NODATA && z[cazi][rng][ang]!=INITIAL_VALUE)
             zval=int_dBZ2float_dBZ(z[cazi][rng][ang]);
         else if (z[cazi][rng][ang]==NODATA)
             zval=NODATA;
         else if (z[cazi][rng][ang]==INITIAL_VALUE)
             zval=INITIAL_VALUE;

//         if (zval!=NODATA && zval!=INITIAL_VALUE) {
             /* bright band effect adjustment */
/*
             if (z_rf!=NODATA) {
                 if (T_top[cazi][rng][ang]<bb_low_degree && T_bottom[cazi][rng][ang]>bb_high_degree) {
                     zdiff=zval-z_rf;
                     if (zdiff>bb_enhance)
                         zval=z_rf;
                 }
             }
         }
*/
         z_ver[ang]=zval;
    }

    if (z_ver[0]>z_value) {
        for (ang=0; ang<nsplitcut-1; ang++) {
             if (z_ver[ang]>45. && z_ver[ang+1]>z_thresh) {
                 zdiff=z_ver[ang]-z_ver[ang+1];
                 if (zdiff>20.)
                     flag=FALSE;
                     break;
             }

             if (z_ver[ang]>z_value && z_ver[ang+1]<z_thresh) {
                 flag=FALSE;
                 break;
             }
        }
    }
    else {
        for (ang=0; ang<nsplitcut-1; ang++) {
             if (z_ver[ang]>z_value && z_ver[ang+1]<z_thresh) {
                 flag=FALSE;
                 break;
             }
        }
    }

    /* further check for all above cuts */
    if (flag==FALSE && z_value<40.) {
        for (ang=nsplitcut; ang<nangle; ang++) {
             if (z_ver[ang]>=z_thresh) {
                 flag=TRUE;
                 break;
             }
        }
    }
    free (z_ver);

    return flag;
}

void identifyWF(short ***z,
                int ***rho,
                float **z_rf,
                short **wf,
                float zrain,
                float znorain)
{
    short **new, **dummy;
    int azi, rng, new_count, dw, change, cover;

    /* new array */
    new=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    dummy=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    new_count=0;
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              dummy[azi][rng]=wf[azi][rng];
              if (wf[azi][rng]==WINDTURBINE)
                  new_count++;
         }
    }

    /* seeded region growing (Adams and Bischof, 1994) */
    dw=1;
    //updateIndicator(dummy,wf,wf,dw,WINDTURBINE,NOTHING);
    while (new_count>0) {
        initialize2dShortArray(new,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,INITIAL_VALUE);
        new_count=0;
        for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
             for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                  if (dummy[azi][rng]==NOTHING) {
                      change=defineWF(wf,z,rho,z_rf,azi,rng,dw,zrain,znorain);
                      if (change==TRUE) {
                          wf[azi][rng]=WINDTURBINE;
                          new[azi][rng]=WINDTURBINE;
                          new_count++;
                      }
                  }
             }
        }
//printf("new_count: %d\n",new_count);

        /* assing "NOTHING" into the surrounding pixels of new "WINDTURBINE" */
        updateIndicator(dummy,wf,new,dw,WINDTURBINE,NOTHING);
    }
    free2dShortArray(dummy,MAXIMUM_AZIMUTH_SUPER);
    free2dShortArray(new,MAXIMUM_AZIMUTH_SUPER);
}

int defineWF(short **type,
             short ***z,
             int ***rho,
             float **z_rf,
             int azi,
             int rng,
             int dw,
             float zrain,
             float znorain)
{
    short zsample[30];
    int count, i, j, iazi, jrng, flag_count, zcount;
    float sum, zmean, zgrowing, cover, raincover_thresh=0.5;

    count=0;
    flag_count=0;
    for (i=azi-dw; i<=azi+dw; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
         if (iazi>=MAXIMUM_AZIMUTH_SUPER)
             iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

         for (j=rng-dw; j<=rng+dw; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              zsample[count]=z[iazi][jrng][0];
              count++;

              if (type[iazi][jrng]==WINDTURBINE)
                  flag_count++;
         }
    }

    if (flag_count==count-1)
        return TRUE;
    else if (flag_count==0)
        return FALSE;
    else {/* calculate Z mean */
        sum=0.0;
        zcount=0;
        for (i=0; i<count; i++) {
             if (zsample[i]!=NODATA && zsample[i]!=INITIAL_VALUE) {
                 sum+=pow(10.,(double)int_dBZ2float_dBZ(zsample[i])/10);
                 //sum+=int_dBZ2float_dBZ(zsample[i]);
                 zcount++;
             }
        }

        if (zcount>0) {
            zmean=10*log10((double)sum/zcount);
            cover=getCoverage(z,rho,azi,rng,0,MAXIMUM_AZIMUTH_SUPER);
            if (cover>=raincover_thresh)
                zgrowing=zrain;
            else
                zgrowing=znorain;

            //zmean=sum/zcount;
            if (zmean>=zgrowing)
                return TRUE;
            else
                return FALSE;
        }
        else
            return FALSE;
    }
}

int checkSurroundingZvalues(short ***z,
                            int azi,
                            int rng,
                            int ang,
                            int dw,
                            float z_thresh)

{
    int i, j, iazi, jrng, count;

    count=0;
    for (i=azi-dw; i<=azi+dw; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
         if (iazi>=MAXIMUM_AZIMUTH_SUPER)
             iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

         for (j=rng-dw; j<=rng+dw; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (z[iazi][jrng][ang]!=NODATA && z[iazi][jrng][ang]!=INITIAL_VALUE) {
                  if (int_dBZ2float_dBZ(z[iazi][jrng][ang])<=z_thresh)
                      count++;
              }
              else
                  count++;
         }
    }
//printf("count: %d\n",count);
    if (count>=(dw+dw+1)*(dw+dw+1)*0.7)
        return TRUE;
    else
        return FALSE;
}

float getinterpZ(float *t,
                 float *z,
                 int count,
                 float tx)
{
    int i;
    float zx;

    zx=0.;
    for (i=0; i<count-1; i++) {
         if ((t[i]-tx)*(t[i+1]-tx)<0) {
             zx=(tx-t[i])*(z[i+1]-z[i])/(t[i+1]-t[i])+z[i];
             break;
         }
    }

    return zx;
}

void classifyRainType(short ***z,
                      int ***rho,
                      float **vil,
                      short **wf,
                      short **type)
{
    short **new, **dummy;
    int azi, rng, new_count, dw, change;
    float vil_thresh=6.5,cover_convection;

    /* new array */
    new=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    dummy=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    new_count=0;
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (vil[azi][rng]>=vil_thresh) {
                  type[azi][rng]=CONVECTIVE;
                  dummy[azi][rng]=CONVECTIVE;
                  new_count++;
              }
              else if (vil[azi][rng]>0 && vil[azi][rng]<vil_thresh) {
                  type[azi][rng]=STRATIFORM;
                  dummy[azi][rng]=STRATIFORM;
              }
         }
    }

    /* seeded region growing (Adams and Bischof, 1994) */
    dw=1;
    updateIndicator(dummy,type,type,dw,CONVECTIVE,STRATIFORM);
    while (new_count>0) {
        initialize2dShortArray(new,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,INITIAL_VALUE);
        new_count=0;
        for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
             for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                  if (dummy[azi][rng]==STRATIFORM) {
                      change=defineConvective(type,vil,azi,rng,dw);
                      if (change==TRUE) {
                          type[azi][rng]=CONVECTIVE;
                          new[azi][rng]=CONVECTIVE;
                          new_count++;
                      }
                  }
             }
        }

        /* assing "STRATIFORM" into the surrounding pixels of new "CONVECTIVE" */
        updateIndicator(dummy,type,new,dw,CONVECTIVE,STRATIFORM);
    }

    /* insert wind farm echoes */
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (wf[azi][rng]==WINDTURBINE){
                  cover_convection=getConvection(type,azi,rng);
                  if (cover_convection<0.5)
                      type[azi][rng]=WINDTURBINE;
              }
         }
    }

    /* check surrounding pixels - smoothing twice */
    checkSurroundingTurbine(type,new);
    checkSurroundingTurbine(type,new);

    free2dShortArray(dummy,MAXIMUM_AZIMUTH_SUPER);
    free2dShortArray(new,MAXIMUM_AZIMUTH_SUPER);
}

void checkSurroundingTurbine(short **type,
                             short **new)
{

    int azi, rng, count, i, j, iazi, jrng, dw;

    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              new[azi][rng]=type[azi][rng];
         }
    }

    dw=1;
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              count=0;
              for (i=azi-dw; i<=azi+dw; i++) {
                   iazi=i;
                   if (iazi<0)
                       iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
                   if (iazi>=MAXIMUM_AZIMUTH_SUPER)
                       iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

                   for (j=rng-dw; j<=rng+dw; j++) {
                        jrng=j;
                        if (jrng<0)
                            continue;
                        if (jrng>=MAXIMUM_RANGE_SUPER)
                            continue;

                        if (new[iazi][jrng]==WINDTURBINE)
                            count++;
                   }
              }

              if (count>4)
                  type[azi][rng]=WINDTURBINE;
         }
    }
}

float getConvection(short **type,
                    int azi,
                    int rng)
{
    int wbin, wazi, w_km; // window size
    int count, i, j, iazi, jrng;
    float R, theta, cover;

    w_km=5;
    wbin=((float)w_km/2)/0.25; // 10 km
    R=0.25*(rng+1); // current distance from the radar
    theta=M_PI*0.5/180; // half beam width (radian)
//    wazi=(int)((int)(R*theta/10)*0.5+0.5); // adaptive window size in azimuth direction
    wazi=(int)(0.5*w_km/(R*theta)+0.5); // adaptive window size in azimuth direction

//    printf("rng: %d, wazi: %d\n",rng,wazi);

    count=0;
    for (i=azi-wazi; i<=azi+wazi; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
         if (iazi>=MAXIMUM_AZIMUTH_SUPER)
             iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

         for (j=rng-wbin; j<=rng+wbin; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              if (type[iazi][jrng]==CONVECTIVE) {
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

int defineConvective(short **type,
                     float **vil,
                     int azi,
                     int rng,
                     int dw)
{
    int count, i, j, iazi, jrng, flag_count, vilcount;
    float vilsample[30], sum, vilmean, vil_thresh=4.0;

    count=0;
    flag_count=0;
    for (i=azi-dw; i<=azi+dw; i++) {
         iazi=i;
         if (iazi<0)
             iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
         if (iazi>=MAXIMUM_AZIMUTH_SUPER)
             iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

         for (j=rng-dw; j<=rng+dw; j++) {
              jrng=j;
              if (jrng<0)
                  continue;
              if (jrng>=MAXIMUM_RANGE_SUPER)
                  continue;

              vilsample[count]=vil[iazi][jrng];
              count++;

              /* if one of surrounding pixel is "convective cell" */
              if (type[iazi][jrng]==CONVECTIVE)
                  flag_count++;
         }
    }

    if (flag_count==count-1)
        return TRUE;
    else if (flag_count==0)
        return FALSE;
    else {/* calculate mean */
        sum=0.0;
        vilcount=0;
        for (i=0; i<count; i++) {
             if (vilsample[i]!=NODATA && vilsample[i]!=INITIAL_VALUE) {
                 sum+=vilsample[i];
                 vilcount++;
             }

        }

        if (vilcount>0) {
            vilmean=sum/vilcount;
            if (vilmean>=vil_thresh)
                return TRUE;
            else
                return FALSE;
        }
        else
            return FALSE;
    }
}

void updateIndicator(short **dummy,
                     short **type,
                     short **new,
                     int dw,
                     int const_convective,
                     int const_stratiform)
{
    int azi, rng, i, j, iazi, jrng;

    initialize2dShortArray(dummy,MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,INITIAL_VALUE);
    for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
         for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
              if (new[azi][rng]==const_convective) {
                  for (i=azi-dw; i<=azi+dw; i++) {
                       iazi=i;
                       if (iazi<0)
                           iazi=iazi+MAXIMUM_AZIMUTH_SUPER;
                       if (iazi>=MAXIMUM_AZIMUTH_SUPER)
                           iazi=iazi-MAXIMUM_AZIMUTH_SUPER;

                       for (j=rng-dw; j<=rng+dw; j++) {
                            jrng=j;
                            if (jrng<0)
                                continue;
                            if (jrng>=MAXIMUM_RANGE_SUPER)
                                continue;

                            if (iazi!=azi && jrng!=rng) {
                                if (type[iazi][jrng]==const_stratiform) /* if it is already convective, don't do it! */
                                    dummy[iazi][jrng]=const_stratiform;
                            }
                       }
                  }
              }
         }
    }
}
