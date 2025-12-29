/*************************************************************************
--> correctSuperRangeEffect.c

    : Estimates hourly VPRs (Vertical Profile of Reflectivity)
      and corrects volume reflectivity data

	1. Creates and updates a temporary file ('vpr.tmp') to store
       hourly reflectivity information
	2. Estimates reflectivity ratios and
       hourly average and azimuth-dependent VPRs
	3. Corrects reflectivity volume data based on the estimated VPRs

    # Function inputs
      - volume_info: volume header information
                     (date & time, beam width, elevation angle, and radar ID)
      - bounding: bounding box information for a specific basin of interest
      - range_parameter: parameters needed to estimate VPR
      - fixedgrid_Z_3d: 3-dimensional reflectivity data based on fixed grid

    # Function output
      - fixedgrid_Z_3d: corrected 3-dimensional reflectivity data

    # Dimension information
      - fixedgrid_Z_3d[MAXIMUM_AZIMUTH][MAXIMUM_RANGE][MAXIMUM_ANGLE]

    Developed by Bongchul Seo and Witold F. Krajewski
    The University of Iowa

    Created : August     2008
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

void buildVPR(struct volume_info header,
              struct bounding box,
              struct range_parameter evpr,
              float **azimuth,
              short **lookup,
              short **rtype,
              short ***z,
              float **vp,
              int *flag_vp)
{
	FILE *ftemp;
    int i, j, k, ns, azi, bin, str, status, nhmax, number_range, number_sector, number_angle, ival;
	int z_thresh=15, number_thresh=5;
    float ***sample;

    number_range=(evpr.end_range-evpr.start_range+1)*4;
    number_sector=360/evpr.sector_move; // number of sector
    number_angle=header.nangle+header.nangle_super+header.nangle_legacy2;
    nhmax=(int)(evpr.hmax/evpr.dh+0.5);

    /* Define array to store data for a specified range to build Vertical Structure */
    sample=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,number_range,number_angle);
    initialize3dFloatArray(sample,MAXIMUM_AZIMUTH_SUPER,number_range,number_angle,NODATA);

    /* Store the data for a specified range */
    for (ns=0; ns<number_angle; ns++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              for (bin=0; bin<number_range; bin++) {
                   if (ns==0)
                       ival=z[azi][bin+(evpr.start_range-1)*4][ns];
                   else
                       ival=z[lookup[azi][ns]][bin+(evpr.start_range-1)*4][ns];

                   if (ival!=NODATA && ival!=INITIAL_VALUE && rtype[azi][bin]==STRATIFORM)
                       sample[azi][bin][ns]=int_dBZ2float_dBZ(ival);
                   else if (ival==NODATA)
                       sample[azi][bin][ns]=NODATA;
                   else if (ival==INITIAL_VALUE)
                       sample[azi][bin][ns]=INITIAL_VALUE;
              }
         }
    }

    /* Set the data to make Vertical Profile ------------------------------------------- */
    int **rhicount;
    double **rhi;

    /* Define array */
    rhi=get2dDoubleArray(nhmax,number_sector+1);
    rhicount=get2dIntArray(nhmax,number_sector+1);

    /* retrieve variability */
/*
    float **rhivar;
    rhivar=get2dFloatArray(nhmax,50000);
*/

    /* assign values to RHI and count */
    for (bin=0; bin<number_range; bin++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              if (sample[azi][bin][0]>=z_thresh && sample[azi][bin][1]>=z_thresh) {
                  for (ns=0; ns<number_angle; ns++) {
                       if (sample[azi][bin][ns]!=NODATA && sample[azi][bin][ns]!=0) {
                           assignValue2RHI(header,evpr,azimuth,pow(10,sample[azi][bin][ns]/10),
                                           azi,bin,ns,rhi,rhicount);
/*                           assignValue2RHIvar(header,evpr,azimuth,pow(10,sample[azi][bin][ns]/10),
                                              azi,bin,ns,rhi,rhicount,rhivar);
*/
                       }
                  }
              }
         }
    }

/*
    FILE *ft;
    ft=fopen("Var_RHI_z.out","w");
    k=rhicount[0][0];
    if (k>50000)
        k=50000;

    for (i=0; i<=k; i++) {
         for (j=0; j<nhmax; j++) {
              if (rhivar[j][i]!=0)
                  fprintf(ft,"%f ",10*log10(rhivar[j][i]));
              else
                  fprintf(ft,"%d ",0);
         }
         fprintf(ft,"\n");
    }
    fclose(ft);
*/
    /* -----------------------------------------------------------------------------------*/

    /* handle temporary file to store data within 1-h ----------------------------- */
    int current_epoch, past_epoch, nscan, timevol[20], ivol, start, **count_archive;
    float fval;
    double **sum_archive;

    /* get current time of volume scan */
    current_epoch=getEpochTime(header);

    /* read temporary file */
    ivol=0;
    if ((ftemp=fopen("vpr.tmp", "r"))!=NULL) {
        fscanf(ftemp,"%d",&nscan);

        /* define arrays */
        count_archive=get2dIntArray(nhmax,nscan*(number_sector+1));
        sum_archive=get2dDoubleArray(nhmax,nscan*(number_sector+1));

        for (i=0; i<nscan; i++) {
             fscanf(ftemp,"%d",&past_epoch);
             if (current_epoch-past_epoch>0 && current_epoch-past_epoch<=HOURINSEC) { // within 1-h
                 timevol[ivol]=past_epoch;
                 ivol++;

                 /* read counts */
                 start=i*(number_sector+1);
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%d",&ival);
                           count_archive[j][start+str]=ival;
                      }
                 }

                 /* read data */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%f",&fval);
                           if (fval!=0.)
                               sum_archive[j][start+str]=pow(10.0,(double)fval/10.0);
                      }
                 }
             }
             else { // not within 1-h domain
                 /* read counts */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%d",&ival);
                      }
                 }

                 /* read data */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%f",&fval);
                      }
                 }
             }
        } // end of loop: number of scans
        fclose(ftemp);
    }

    /* write temporary file */
    ftemp=fopen("vpr.tmp", "w");
    fprintf(ftemp,"%d\n",ivol+1);

    /* write archive arrays */
    if (ivol!=0) {
        for (i=0; i<ivol; i++) {
             fprintf(ftemp,"%d\n",timevol[i]);

             /* write counts of pairs between 1st and 2nd elevation angles */
             start=i*(number_sector+1);
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       fprintf(ftemp,"%d ",count_archive[j][start+str]);
                  }
                  fprintf(ftemp,"\n");
             }

             /* write sum */
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       if (sum_archive[j][start+str]==0)
                           fprintf(ftemp,"%d ",INITIAL_VALUE);
                       else
                           fprintf(ftemp,"%g ",10*log10(sum_archive[j][start+str]));
                  }
                  fprintf(ftemp,"\n");
             }
        }
    }

    /* write current one */
    fprintf(ftemp,"%d\n",current_epoch);
    for (str=0; str<number_sector+1; str++) {
         for (j=0; j<nhmax; j++) {
              fprintf(ftemp,"%d ",rhicount[j][str]);
         }
         fprintf(ftemp,"\n");
    }

    for (str=0; str<number_sector+1; str++) {
         for (j=0; j<nhmax; j++) {
              if (rhi[j][str]==0)
                  fprintf(ftemp,"%d ",INITIAL_VALUE);
              else
                  fprintf(ftemp,"%g ",10*log10(rhi[j][str]));
         }
         fprintf(ftemp,"\n");
    }
    fclose(ftemp);
    /* -----------------------------------------------------------------------------------*/

    /* combine current and past (within 1-h) ones */
    if (ivol!=0) {
        for (i=0; i<ivol; i++) {
             /* aggreate counts */
             start=i*(number_sector+1);
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       rhicount[j][str]+=count_archive[j][start+str];
                  }
             }

             /* aggreate sum */
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       rhi[j][str]+=sum_archive[j][start+str];
                  }
             }
        } // end of loop: number of scans
        status=free2dIntArray(count_archive,nhmax);
        status=free2dDoubleArray(sum_archive,nhmax);
    }

    /* compute ratio --------------------------------------------------------------*/
    float rf;

    /* ratio */
    for (i=0; i<nhmax; i++) {
         for (str=0; str<number_sector+1; str++) {
              rf=rhi[0][str]/rhicount[0][str];
              if (rhicount[i][str]>number_thresh*(ivol+1))
                  vp[i][str]=(rhi[i][str]/rhicount[i][str])/rf;
              else
                  vp[i][str]=0.;
         }
    }
    status=free2dDoubleArray(rhi,nhmax);
    status=free2dIntArray(rhicount,nhmax);

    /* fill missing gap below 6 km --------------------------------------------*/
    int nh6, idown, iup;

    nh6=(int)(6./evpr.dh+0.5);

    /* check values if vpr is available */
    for (i=nh6-1; i<nhmax; i++) {
         if (vp[i][0]!=0) {
             *flag_vp=TRUE;
             break;
         }
    }

    if (*flag_vp==TRUE) {
        for (i=1; i<nh6; i++) {
             if (vp[i][0]==0) {
                 /* find the nearest lower bin that contains a vpr value */
                 idown=0;
                 for (j=i-1; j>=0; j--) {
                      if (vp[j][0]!=0) {
                          idown=j;
                          break;
                      }
                 }

                 /*find the nearest higher bin that contains a vpr value */
                 for (j=i+1; j<nhmax; j++) {
                      if (vp[j][0]!=0) {
                          iup=j;
                          break;
                      }
                 }
                 vp[i][0]=(i-idown)*(vp[iup][0]-vp[idown][0])/(iup-idown)+vp[idown][0];
             }
        }
    }

    /* file out */
    char fvpr[50]="", stemp[30]="";

    sprintf(stemp,"%d",current_epoch);
    strcpy(fvpr,"vpr_");
    strcat(fvpr,stemp);
    strcat(fvpr,".out");
    ftemp=fopen(fvpr,"w");
    for (i=0; i<nhmax; i++) {
         if (vp[i][0]!=0)
             fprintf(ftemp,"%.2f, %f\n",evpr.dh*i+evpr.dh*0.5,10*log10(vp[i][0]));
         else
             fprintf(ftemp,"%.2f, %d\n",evpr.dh*i+evpr.dh*0.5,0);
    }
    fclose(ftemp);
}

void buildVPDP(struct volume_info header,
               struct bounding box,
               struct range_parameter evpr,
               float **azimuth,
               short **lookup,
               short **rtype,
               short ***z,
               int ***dp,
               int dpname,
               float **vp,
               int *flag_vp)
{
	FILE *ftemp;
    int i, j, k, ns, azi, bin, str, status, nhmax, number_range, number_sector, number_angle, ival;
	int z_thresh=15, number_thresh=5;
    float ***sample, ***zsample;

    number_range=(evpr.end_range-evpr.start_range+1)*4;
    number_sector=360/evpr.sector_move; // number of sector
    number_angle=header.nangle+header.nangle_super+header.nangle_legacy2;
    nhmax=(int)(evpr.hmax/evpr.dh+0.5);

    /* Define array to store data for a specified range to build Vertical Structure */
    sample=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,number_range,number_angle);
    zsample=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,number_range,number_angle);
    initialize3dFloatArray(sample,MAXIMUM_AZIMUTH_SUPER,number_range,number_angle,NODATA);

    /* Store reflectivity data for a specified range */
    for (ns=0; ns<number_angle; ns++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              for (bin=0; bin<number_range; bin++) {
                   if (ns==0)
                       ival=z[azi][bin+(evpr.start_range-1)*4][ns];
                   else
                       ival=z[lookup[azi][ns]][bin+(evpr.start_range-1)*4][ns];

                   if (ival!=NODATA && ival!=INITIAL_VALUE && rtype[azi][bin]==STRATIFORM)
                       zsample[azi][bin][ns]=int_dBZ2float_dBZ(ival);
                   else if (ival==NODATA)
                       zsample[azi][bin][ns]=NODATA;
                   else if (ival==INITIAL_VALUE)
                       zsample[azi][bin][ns]=INITIAL_VALUE;
              }
         }
    }

    /* Store DP data for a specified range */
    for (ns=0; ns<number_angle; ns++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              for (bin=0; bin<number_range; bin++) {
                   if (ns==0)
                       ival=dp[azi][bin+(evpr.start_range-1)*4][ns];
                   else
                       ival=dp[lookup[azi][ns]][bin+(evpr.start_range-1)*4][ns];

                   if (ival!=NODATA && ival!=INITIAL_VALUE)
                       sample[azi][bin][ns]=(float)ival/DP_SCALEING;
                   else if (ival==NODATA)
                       sample[azi][bin][ns]=NODATA;
                   else if (ival==INITIAL_VALUE)
                       sample[azi][bin][ns]=INITIAL_VALUE;
              }
         }
    }

    /* Set the data to make Vertical Profile ------------------------------------------- */
    int **rhicount;
    double **rhi;

    /* Define array */
    rhi=get2dDoubleArray(nhmax,number_sector+1);
    rhicount=get2dIntArray(nhmax,number_sector+1);

    /* retrieve variability */
/*
    float **rhivar;
    rhivar=get2dFloatArray(nhmax,50000);
*/

    /* assign values to RHI and count */
    for (bin=0; bin<number_range; bin++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
              if (zsample[azi][bin][0]>=z_thresh && zsample[azi][bin][1]>=z_thresh) {
                  for (ns=0; ns<number_angle; ns++) {
                       if (sample[azi][bin][ns]!=NODATA && sample[azi][bin][ns]!=0) {
                           assignValue2RHI(header,evpr,azimuth,sample[azi][bin][ns],
                                           azi,bin,ns,rhi,rhicount);
/*                           assignValue2RHIvar(header,evpr,azimuth,sample[azi][bin][ns],
                                              azi,bin,ns,rhi,rhicount,rhivar);
*/
                       }
                  }
              }
         }
    }

/*
    FILE *ft;
    if (dpname==RH_INDEX)
        ft=fopen("Var_RHI_rho.out","w");
    else if (dpname==DR_INDEX)
        ft=fopen("Var_RHI_zdr.out","w");
    else if (dpname==PH_INDEX)
        ft=fopen("Var_RHI_kdp.out","w");

    k=rhicount[0][0];
    if (k>50000)
        k=50000;

    for (i=0; i<=k; i++) {
         for (j=0; j<nhmax; j++) {
              if (rhivar[j][i]!=0)
                  fprintf(ft,"%f ",rhivar[j][i]);
              else
                  fprintf(ft,"%d ",0);
         }
         fprintf(ft,"\n");
    }
    fclose(ft);
*/

    /* -----------------------------------------------------------------------------------*/

    /* handle temporary file to store data within 1-h ----------------------------- */
    int current_epoch, past_epoch, nscan, timevol[20], ivol, start, **count_archive;
    float fval;
    double **sum_archive;
    char nametemp[20]="";

    /* temporary file name */
    if (dpname==RH_INDEX)
        strcpy(nametemp,"vprho.tmp");
    else if (dpname==DR_INDEX)
        strcpy(nametemp,"vpzdr.tmp");
    else if (dpname==PH_INDEX)
        strcpy(nametemp,"vpkdp.tmp");

    /* get current time of volume scan */
    current_epoch=getEpochTime(header);

    /* read temporary file */
    ivol=0;
    if ((ftemp=fopen(nametemp, "r"))!=NULL) {
        fscanf(ftemp,"%d",&nscan);

        /* define arrays */
        count_archive=get2dIntArray(nhmax,nscan*(number_sector+1));
        sum_archive=get2dDoubleArray(nhmax,nscan*(number_sector+1));

        for (i=0; i<nscan; i++) {
             fscanf(ftemp,"%d",&past_epoch);
             if (current_epoch-past_epoch>0 && current_epoch-past_epoch<=HOURINSEC) { // within 1-h
                 timevol[ivol]=past_epoch;
                 ivol++;

                 /* read counts */
                 start=i*(number_sector+1);
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%d",&ival);
                           count_archive[j][start+str]=ival;
                      }
                 }

                 /* read data */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%f",&fval);
                           if (fval!=0.)
                               sum_archive[j][start+str]=fval;
                      }
                 }
             }
             else { // not within 1-h domain
                 /* read counts */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%d",&ival);
                      }
                 }

                 /* read data */
                 for (str=0; str<number_sector+1; str++) {
                      for (j=0; j<nhmax; j++) {
                           fscanf(ftemp,"%f",&fval);
                      }
                 }
             }
        } // end of loop: number of scans
        fclose(ftemp);
    }

    /* write temporary file */
    ftemp=fopen(nametemp, "w");
    fprintf(ftemp,"%d\n",ivol+1);

    /* write archive arrays */
    if (ivol!=0) {
        for (i=0; i<ivol; i++) {
             fprintf(ftemp,"%d\n",timevol[i]);

             /* write counts of pairs between 1st and 2nd elevation angles */
             start=i*(number_sector+1);
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       fprintf(ftemp,"%d ",count_archive[j][start+str]);
                  }
                  fprintf(ftemp,"\n");
             }

             /* write sum */
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       if (sum_archive[j][start+str]==0)
                           fprintf(ftemp,"%d ",INITIAL_VALUE);
                       else
                           fprintf(ftemp,"%g ",sum_archive[j][start+str]);
                  }
                  fprintf(ftemp,"\n");
             }
        }
    }

    /* write current one */
    fprintf(ftemp,"%d\n",current_epoch);
    for (str=0; str<number_sector+1; str++) {
         for (j=0; j<nhmax; j++) {
              fprintf(ftemp,"%d ",rhicount[j][str]);
         }
         fprintf(ftemp,"\n");
    }

    for (str=0; str<number_sector+1; str++) {
         for (j=0; j<nhmax; j++) {
              if (rhi[j][str]==0)
                  fprintf(ftemp,"%d ",INITIAL_VALUE);
              else
                  fprintf(ftemp,"%g ",rhi[j][str]);
         }
         fprintf(ftemp,"\n");
    }
    fclose(ftemp);
    /* -----------------------------------------------------------------------------------*/

    /* combine current and past (within 1-h) ones */
    if (ivol!=0) {
        for (i=0; i<ivol; i++) {
             /* aggreate counts */
             start=i*(number_sector+1);
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       rhicount[j][str]+=count_archive[j][start+str];
                  }
             }

             /* aggreate sum */
             for (str=0; str<number_sector+1; str++) {
                  for (j=0; j<nhmax; j++) {
                       rhi[j][str]+=sum_archive[j][start+str];
                  }
             }
        } // end of loop: number of scans
        status=free2dIntArray(count_archive,nhmax);
        status=free2dDoubleArray(sum_archive,nhmax);
    }

    /* averaging */
    for (i=0; i<nhmax; i++) {
         for (str=0; str<number_sector+1; str++) {
              if (rhicount[i][str]>number_thresh*(ivol+1))
                  vp[i][str]=rhi[i][str]/rhicount[i][str];
              else
                  vp[i][str]=0.;
         }
    }
    status=free2dDoubleArray(rhi,nhmax);
    status=free2dIntArray(rhicount,nhmax);

    /* fill missing gap below 6 km --------------------------------------------*/
    int nh6, idown, iup;

    nh6=(int)(6./evpr.dh+0.5);

    /* check values if vpr is available */
    for (i=nh6-1; i<nhmax; i++) {
         if (vp[i][0]!=0) {
             *flag_vp=TRUE;
             break;
         }
    }

    if (*flag_vp==TRUE) {
        for (i=1; i<nh6; i++) {
             if (vp[i][0]==0) {
                 /* find the nearest lower bin that contains a vpr value */
                 idown=0;
                 for (j=i-1; j>=0; j--) {
                      if (vp[j][0]!=0) {
                          idown=j;
                          break;
                      }
                 }

                 /*find the nearest higher bin that contains a vpr value */
                 for (j=i+1; j<nhmax; j++) {
                      if (vp[j][0]!=0) {
                          iup=j;
                          break;
                      }
                 }
                 vp[i][0]=(i-idown)*(vp[iup][0]-vp[idown][0])/(iup-idown)+vp[idown][0];
             }
        }
    }

    for (i=0; i<nhmax; i++) {
         if (vp[i][0]!=0)
             printf("i: %d,h: %.1f, %f\n",i,evpr.dh*i,vp[i][0]);
         else
             printf("i: %d, h: %.1f, %d\n",i,evpr.dh*i,0);
    }
}

void assignValue2RHI(struct volume_info header,
                     struct range_parameter evpr,
                     float **azimuth,
                     double dval,
                     int azi,
                     int bin,
                     int ns,
                     double **rhi,
                     int **rhicount)
{
    int i, j, str, str2, noverlap, number_sector, nhmax, ith, ibh;
    float upangle, lowangle, htop, hbot;

    number_sector=360/evpr.sector_move; // number of sector
    nhmax=(int)(evpr.hmax/evpr.dh+0.5);

    /* find sector numbers */
    str=(int)(azimuth[azi][0]/evpr.sector_move)+1; // start from 1
    noverlap=evpr.sector_width/evpr.sector_move;

    /* assign value to corresponding vertical grids */
    upangle=header.angle[ns]+header.beam_width/2; // top
    htop=computeBeamHeight(upangle,bin*0.25+evpr.start_range+0.125);
    lowangle=header.angle[ns]-header.beam_width/2; // bottom
    hbot=computeBeamHeight(lowangle,bin*0.25+evpr.start_range+0.125);
    ith=(int)(htop/evpr.dh+0.5);
    if (ith>nhmax)
	    ith=nhmax;
    ibh=(int)(hbot/evpr.dh+0.5)+1;
	if (ibh<=0)
	    ibh=1;

    if (ns==0) { // base scan
        for (i=0; i<ith; i++) {
             rhi[i][0]+=dval;
             rhicount[i][0]++;
             for (j=str; j>str-noverlap; j--) {
                  str2=j;
                  if (str2<=0) {
                      str2+=number_sector;
                      rhi[i][str2]+=dval;
                      rhicount[i][str2]++;
                  }
             }
        }
    }
    else {
	    if (ibh-1<nhmax) {
	        for (i=ibh-1; i<ith; i++) {
                 rhi[i][0]+=dval;
                 rhicount[i][0]++;
                 for (j=str; j>str-noverlap; j--) {
                      str2=j;
                      if (str2<=0) {
                          str2+=number_sector;
                          rhi[i][str2]+=dval;
                          rhicount[i][str2]++;
                      }
                 }
	        }
	    }
    }
}

void assignValue2RHIvar(struct volume_info header,
                        struct range_parameter evpr,
                        float **azimuth,
                        double dval,
                        int azi,
                        int bin,
                        int ns,
                        double **rhi,
                        int **rhicount,
                        float **varrhi)
{
    int i, j, str, str2, noverlap, number_sector, nhmax, ith, ibh;
    float upangle, lowangle, htop, hbot;

    number_sector=360/evpr.sector_move; // number of sector
    nhmax=(int)(evpr.hmax/evpr.dh+0.5);

    /* find sector numbers */
    str=(int)(azimuth[azi][0]/evpr.sector_move)+1; // start from 1
    noverlap=evpr.sector_width/evpr.sector_move;

    /* assign value to corresponding vertical grids */
    upangle=header.angle[ns]+header.beam_width/2; // top
    htop=computeBeamHeight(upangle,bin*0.25+evpr.start_range+0.125);
    lowangle=header.angle[ns]-header.beam_width/2; // bottom
    hbot=computeBeamHeight(lowangle,bin*0.25+evpr.start_range+0.125);
    ith=(int)(htop/evpr.dh+0.5);
    if (ith>nhmax)
	    ith=nhmax;
    ibh=(int)(hbot/evpr.dh+0.5)+1;
	if (ibh<=0)
	    ibh=1;

    if (ns==0) { // base scan
        for (i=0; i<ith; i++) {
             if (rhicount[i][0]<50000)
                 varrhi[i][rhicount[i][0]]=dval;

             rhi[i][0]+=dval;
             rhicount[i][0]++;
             for (j=str; j>str-noverlap; j--) {
                  str2=j;
                  if (str2<=0) {
                      str2+=number_sector;
                      rhi[i][str2]+=dval;
                      rhicount[i][str2]++;
                  }
             }
        }
    }
    else {
	    if (ibh-1<nhmax) {
	        for (i=ibh-1; i<ith; i++) {
                 if (rhicount[i][0]<50000)
                     varrhi[i][rhicount[i][0]]=dval;

                 varrhi[i][rhicount[i][0]]=dval;
                 rhi[i][0]+=dval;
                 rhicount[i][0]++;
                 for (j=str; j>str-noverlap; j--) {
                      str2=j;
                      if (str2<=0) {
                          str2+=number_sector;
                          rhi[i][str2]+=dval;
                          rhicount[i][str2]++;
                      }
                 }
	        }
	    }
    }
}
