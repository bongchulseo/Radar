/********************************************************************************************************
--> getRain.c

    : Manages subfuctions to generate a rain rate product for an individual radar domain

    # Procedures
    1. read Level II data: 'readLevel2'
    2. build melting layer information based on temperature soundings: 'getMeltinglayer'
    3. rain type classification (convective/stratiform)
       : 'computeVIL' and 'classifyRainType'
    4. data QC (noise filtering and averaging)
    5. Kdp and VPR (optional, but currently disabled) estiamtion
    6. rainfall estimation using a variety of estimators (default: specific attenuation)
    7. build a PPI/CAPPI product based on a fixed grid (0.5 degrees by 250 m --> 720*920)
       : 'getRainRateCAPPI' and 'get2DFixedGrid'

Added:    Rkdp--> October 2018
          HCA --> February 2019
          default R-Z (when no temperature sounding is available) -> March 2019

Modified: radial feature near the radar site was corrected --> May 2020
          segfault due to missing rays (caused error in 'getMeltinglayer') was corrected --> June 2020

***********************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "libphydro.h"

void getRain(char *fname,
             char *nwppath,
             int angle_number,
             float cappi_h,
             float tower_h,
             int qc_flag,
             int range_flag,
             int norain,
             struct volume_info* header,
             struct bounding box,
             struct range_parameter evpr,
             float **rainrate,
             int estimator,
             int *status)
{
    float **azimuth, ***polargrid_rain_3d, alpha;
    short ***polargrid_Z_3d, **azi_lookup;
    int ***polargrid_RHO_3d, ***polargrid_ZDR_3d, ***polargrid_PHI_3d, ***polargrid_PHI_short,
        ***polargrid_PHI_long, ***polargrid_KDP_short, ***polargrid_KDP_long, ***flag_echo_3d,
        max_nangle, cappi_nangle=9, zwin=3, dpwin=5, shortwin=9, longwin=25, azi, cazi, rng, ang;
    short ***melting, **Rtype, ***Rtype_3d, **wf;
    float **VIL, **vpr, **vprho, **vpzdr, **vpphi, ***air_top, ***air_bottom, *temp_azi, bw=0.95;
    int nhmax, number_sector, *flag_vpr, *flag_vprho, *flag_vpzdr, *flag_vpphi, *flag_nwp,
        ini_vpr=FALSE, ini_vpzdr=FALSE, ini_vprho=FALSE, ini_vpphi=FALSE, ini_nwp=TRUE;

    /* initialize vp flags */
    flag_vpr=&ini_vpr;
    flag_vprho=&ini_vprho;
    flag_vpzdr=&ini_vpzdr;
    flag_vpphi=&ini_vpphi;

    /* initialize nwp (temperature sounding) */
    flag_nwp=&ini_nwp;

    /* define arrays ******************************************************/
    azimuth=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_ANGLE_SUPER);
    VIL=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    Rtype=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    Rtype_3d=get3dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    wf=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);
    polargrid_Z_3d=get3dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_RHO_3d=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_ZDR_3d=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_PHI_3d=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_PHI_short=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_PHI_long=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_KDP_short=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_KDP_long=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    polargrid_rain_3d=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    flag_echo_3d=get3dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    melting=get3dShortArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    air_top=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    air_bottom=get3dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER);
    /*--------------------------------------------------------------------------------------------------*/

    /* read Level II data */
    readLevel2(fname,azimuth,polargrid_Z_3d,polargrid_RHO_3d,
               polargrid_ZDR_3d,polargrid_PHI_3d,header);

    /* find maximum nuber of elevation angles */
    max_nangle=header->nangle+header->nangle_super+header->nangle_legacy2;
    if (max_nangle>MAXIMUM_ANGLE_SUPER)
        max_nangle=MAXIMUM_ANGLE_SUPER;

    /* create lookup for a vertical alignment of azimuth */
    azi_lookup=get2dShortArray(MAXIMUM_AZIMUTH_SUPER,max_nangle);
    getvaziLookup(*header,max_nangle,azimuth,azi_lookup);

    /* get melting layer and air temperature information */
    getMeltinglayer(*header,nwppath,azimuth,melting,air_top,air_bottom,flag_nwp);

    /* print melting layer information */
    float sum;
    int count;

    if (*flag_nwp==TRUE) {
        sum=0.;
        count=0;
        for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
             for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                  if (melting[azi][rng][0]==MELTING_FLAG || melting[azi][rng][0]==FREEZING_FLAG) {
                      sum=sum+(float)(rng-1)/4+0.125;
                      count++;
                      break;
                  }
                  else if (rng==MAXIMUM_RANGE_SUPER-1) {
                      sum=sum+(float)rng/4+0.125;
                      count++;
                      break;
                  }
             }
        }
        printf("average melting range: %f km\n",sum/count);
    }

    if (*flag_nwp==TRUE) {
        computeVIL(*header,polargrid_Z_3d,polargrid_RHO_3d,polargrid_PHI_3d,air_top,air_bottom,azimuth,azi_lookup,VIL,wf);
        classifyRainType(polargrid_Z_3d,polargrid_RHO_3d,VIL,wf,Rtype);
        for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
             for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                  for (ang=0; ang<max_nangle; ang++) {
                       if (ang==0)
                           cazi=azi;
                       else
                           cazi=azi_lookup[azi][ang];
                       Rtype_3d[cazi][rng][ang]=Rtype[azi][rng];

                       if (wf[azi][rng]==WINDTURBINE)
                           flag_echo_3d[cazi][rng][ang]=CLUTTER;

                  }
             }
        }
    }

    /* generate classification file */
/*
    char temp_str[50]="",epoch_str[30]="";
    int current_epoch;

    current_epoch=getEpochTime(*header);
    sprintf(epoch_str,"%d",current_epoch);
    strcpy(temp_str,"RainType_");
    strcat(temp_str,header->radar_id);
    strcat(temp_str,"_");
    strcat(temp_str,epoch_str);
    strcat(temp_str,".out");
*/

    if (max_nangle!=0) {
        /* copy differential phase data */
        copy3Darray(polargrid_PHI_3d,polargrid_PHI_short);
	    copy3Darray(polargrid_PHI_3d,polargrid_PHI_long);

        if (range_flag==TRUE || cappi_h>0) { /* multiple elevation processing */
            for (ang=1; ang<=max_nangle; ang++) {
                 /* 1. unfold Phidp */
                 unfoldPhidp(*header,polargrid_PHI_3d,polargrid_RHO_3d,ang,DP_SCALEING);

                 /* 2. peform QC */
                 if (qc_flag==TRUE) {
                     startQualityControl(*header,polargrid_Z_3d,polargrid_RHO_3d,polargrid_ZDR_3d,
                                         polargrid_PHI_3d,flag_echo_3d,ang);
                     smoothClutter(*header,flag_echo_3d,polargrid_Z_3d,ang);
                     maskObservables(*header,flag_echo_3d,polargrid_Z_3d,polargrid_RHO_3d,
                                     polargrid_ZDR_3d,polargrid_PHI_3d,ang);
                 }

                 /* 3. filter DP observables */
                 averageZ(*header,polargrid_Z_3d,zwin,ang); // average Z
                 averageDPvariables(*header,polargrid_RHO_3d,dpwin,ang,DP_SCALEING); // average Rho
                 averageDPvariables(*header,polargrid_ZDR_3d,dpwin,ang,DP_SCALEING); // average Zdr
                 medianDPvariables(*header,polargrid_PHI_3d,dpwin,ang,DP_SCALEING); // median filter Phidp
                 averageDPvariables(*header,polargrid_PHI_3d,dpwin,ang,DP_SCALEING); // average Phidp

                 /* 4. split Phidp into short and long averaging windows */
                 medianDPvariables(*header,polargrid_PHI_short,shortwin,ang,DP_SCALEING); // median filter Phidp
                 medianDPvariables(*header,polargrid_PHI_long,longwin,ang,DP_SCALEING); // median filter Phidp
                 averageDPvariables(*header,polargrid_PHI_short,shortwin,ang,DP_SCALEING);
                 averageDPvariables(*header,polargrid_PHI_long,longwin,ang,DP_SCALEING);

                 /* 5. flag echo and fill gaps in Phidp that have no echo */
                 if (*flag_nwp==TRUE) {
                     defineEcho(*header,polargrid_Z_3d,polargrid_RHO_3d,polargrid_ZDR_3d,air_top,air_bottom,
                                flag_echo_3d,melting,Rtype_3d,ang,DP_SCALEING);
                     fillGaps(*header,polargrid_PHI_short,flag_echo_3d,shortwin,ang,DP_SCALEING);
                     fillGaps(*header,polargrid_PHI_long,flag_echo_3d,longwin,ang,DP_SCALEING);
                 }

                 /* 6. Kdp estimation */
                 computeKdp(*header,polargrid_PHI_short,polargrid_RHO_3d,polargrid_Z_3d,polargrid_KDP_short,
                            shortwin,ang,DP_SCALEING);
                 computeKdp(*header,polargrid_PHI_long,polargrid_RHO_3d,polargrid_Z_3d,polargrid_KDP_long,
                            longwin,ang,DP_SCALEING);
                 adjustKdp(*header,ang,polargrid_Z_3d,polargrid_KDP_short,polargrid_KDP_long);
            }

            if (range_flag==TRUE && *flag_nwp==TRUE) { /* build vertical structure */

                /* vertical profile */
                nhmax=(int)(evpr.hmax/evpr.dh+0.5);
                number_sector=360/evpr.sector_move;
                vpr=get2dFloatArray(nhmax,number_sector+1); // first row: overall
                vprho=get2dFloatArray(nhmax,number_sector+1);
                vpzdr=get2dFloatArray(nhmax,number_sector+1);
                vpphi=get2dFloatArray(nhmax,number_sector+1);

                /* vertical profile */
                buildVPR(*header,box,evpr,azimuth,azi_lookup,Rtype,polargrid_Z_3d,vpr,flag_vpr);
//                buildVPDP(*header,box,evpr,azimuth,azi_lookup,Rtype,polargrid_Z_3d,polargrid_RHO_3d,RH_INDEX,vprho,flag_vprho);
//                buildVPDP(*header,box,evpr,azimuth,azi_lookup,Rtype,polargrid_Z_3d,polargrid_ZDR_3d,DR_INDEX,vpzdr,flag_vpzdr);
//                buildVPDP(*header,box,evpr,azimuth,azi_lookup,Rtype,polargrid_Z_3d,polargrid_KDP_long,PH_INDEX,vpphi,flag_vpphi);
            }
        }
        else { /* single elevation processing */
			
            /* 1. unfold Phidp */
            unfoldPhidp(*header,polargrid_PHI_3d,polargrid_RHO_3d,angle_number,DP_SCALEING);

            /* 2. peform QC */
            if (qc_flag==TRUE) {
                startQualityControl(*header,polargrid_Z_3d,polargrid_RHO_3d,polargrid_ZDR_3d,
                                    polargrid_PHI_3d,flag_echo_3d,angle_number);
                smoothClutter(*header,flag_echo_3d,polargrid_Z_3d,angle_number);
                maskObservables(*header,flag_echo_3d,polargrid_Z_3d,polargrid_RHO_3d,
                                polargrid_ZDR_3d,polargrid_PHI_3d,angle_number);
            }

            /* 3. filter DP observables */
            averageZ(*header,polargrid_Z_3d,zwin,angle_number); // average Z
            averageDPvariables(*header,polargrid_RHO_3d,dpwin,angle_number,DP_SCALEING); // average Rho
            averageDPvariables(*header,polargrid_ZDR_3d,dpwin,angle_number,DP_SCALEING); // average Zdr
            medianDPvariables(*header,polargrid_PHI_3d,dpwin,angle_number,DP_SCALEING); // median filter Phidp
            averageDPvariables(*header,polargrid_PHI_3d,dpwin,angle_number,DP_SCALEING); // average Phidp

            /* 4. split Phidp into short and long averaging windows */
            medianDPvariables(*header,polargrid_PHI_short,shortwin,angle_number,DP_SCALEING); // median filter Phidp
            medianDPvariables(*header,polargrid_PHI_long,longwin,angle_number,DP_SCALEING); // median filter Phidp
            averageDPvariables(*header,polargrid_PHI_short,shortwin,angle_number,DP_SCALEING);
            averageDPvariables(*header,polargrid_PHI_long,longwin,angle_number,DP_SCALEING);

            /* 5. flag echo and fill gaps in Phidp that have no echo */
            if (*flag_nwp==TRUE) {
                defineEcho(*header,polargrid_Z_3d,polargrid_RHO_3d,polargrid_ZDR_3d,air_top,air_bottom,
                           flag_echo_3d,melting,Rtype_3d,angle_number,DP_SCALEING);
                fillGaps(*header,polargrid_PHI_short,flag_echo_3d,shortwin,angle_number,DP_SCALEING);
                fillGaps(*header,polargrid_PHI_long,flag_echo_3d,longwin,angle_number,DP_SCALEING);
            }

            /* 6. Kdp estimation */
            computeKdp(*header,polargrid_PHI_short,polargrid_RHO_3d,polargrid_Z_3d,polargrid_KDP_short,
                       shortwin,angle_number,DP_SCALEING);
            computeKdp(*header,polargrid_PHI_long,polargrid_RHO_3d,polargrid_Z_3d,polargrid_KDP_long,
                       longwin,angle_number,DP_SCALEING);
            adjustKdp(*header,ang,polargrid_Z_3d,polargrid_KDP_short,polargrid_KDP_long);
        }

        /* 7. compute specific attenuation and estimate rain rate */
        if (estimator==RA && *flag_nwp==TRUE) {
            alpha=computeAlpha(*header,polargrid_Z_3d,polargrid_ZDR_3d,polargrid_RHO_3d,melting,DP_SCALEING);
            printf("alpha: %f\n",alpha);
        }

        if (*flag_nwp==FALSE) {// no NWP, then apply simple Z-R (no VPR correction)
            printf("Will apply a default estimator!!\n");
            if (angle_number>=1 && angle_number<=9) { /* for angle scan */
                estimateRZdefault(*header,polargrid_Z_3d,flag_echo_3d,polargrid_rain_3d,angle_number);
            }
            else if (angle_number==10) { /* for CAPPI */
                if (max_nangle<cappi_nangle)
                    cappi_nangle=max_nangle;
                for (ang=1; ang<=cappi_nangle; ang++) {
                     estimateRZdefault(*header,polargrid_Z_3d,flag_echo_3d,polargrid_rain_3d,ang);
                }
            }
        }

        if (*flag_nwp==TRUE && *flag_vpr==FALSE) {
            if (angle_number>=1 && angle_number<=9) { /* for angle scan */
                if (estimator==RA) {
                    estimateRA(*header,polargrid_Z_3d,polargrid_PHI_long,polargrid_KDP_long,flag_echo_3d,Rtype_3d,
                               melting,azimuth,polargrid_rain_3d,alpha,angle_number,DP_SCALEING);
                }
                else if (estimator==RZC) {
                    estimateRZC(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,angle_number);
                }
                else if (estimator==RZCS) {
                    estimateRZCS(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,angle_number);
                }
                else if (estimator==RZDR) {
                    estimateRZDR(*header,polargrid_Z_3d,polargrid_ZDR_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,angle_number,DP_SCALEING);
                }
                else if (estimator==RKDP) {
                    estimateRKDP(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,angle_number,DP_SCALEING);
                }
            }
            else if (angle_number==10) { /* for CAPPI */
                if (max_nangle<cappi_nangle)
                    cappi_nangle=max_nangle;
                for (ang=1; ang<=cappi_nangle; ang++) {
                     if (estimator==RA) {
                         estimateRA(*header,polargrid_Z_3d,polargrid_PHI_long,polargrid_KDP_long,flag_echo_3d,Rtype_3d,
                                    melting,azimuth,polargrid_rain_3d,alpha,ang,DP_SCALEING);
                     }
                     else if (estimator==RZC) {
                         estimateRZC(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,ang);
                     }
                     else if (estimator==RZCS) {
                         estimateRZCS(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,ang);
                     }
                     else if (estimator==RZDR) {
                         estimateRZDR(*header,polargrid_Z_3d,polargrid_ZDR_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,ang,DP_SCALEING);
                     }
                     else if (estimator==RKDP) {
                         estimateRKDP(*header,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,polargrid_rain_3d,ang,DP_SCALEING);
                     }
                }
            }
        }
        else if (*flag_nwp==TRUE && range_flag==TRUE && *flag_vpr==TRUE) {
            if (angle_number>=1 && angle_number<=9) { /* for angle scan */
                if (estimator==RA) {
                    estimateRAwVPR(*header,evpr,polargrid_Z_3d,polargrid_PHI_long,polargrid_KDP_long,flag_echo_3d,Rtype_3d,
                                   melting,vpr,azimuth,alpha,header->angle,angle_number,DP_SCALEING,polargrid_rain_3d);
                }
                else if (estimator==RZC) {
                    estimateRZCwVPR(*header,evpr,polargrid_Z_3d,flag_echo_3d,Rtype_3d,vpr,header->angle,
                                    angle_number,polargrid_rain_3d);
                }
                else if (estimator==RZCS) {
                    estimateRZCSwVPR(*header,evpr,polargrid_Z_3d,flag_echo_3d,Rtype_3d,vpr,header->angle,
                                     angle_number,polargrid_rain_3d);
                }
                else if (estimator==RZDR) {
                    estimateRZDRwVPR(*header,evpr,polargrid_Z_3d,polargrid_ZDR_3d,flag_echo_3d,Rtype_3d,melting,vpr,
                                     header->angle,angle_number,DP_SCALEING,polargrid_rain_3d);
                }
                else if (estimator==RKDP) {
                    estimateRKDPwVPR(*header,evpr,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,melting,vpr,
                                     header->angle,angle_number,DP_SCALEING,polargrid_rain_3d);
                }
            }
            else if (angle_number==10) { /* for CAPPI */
                if (max_nangle<cappi_nangle)
                    cappi_nangle=max_nangle;
                for (ang=1; ang<=cappi_nangle; ang++) {
                     if (estimator==RA) {
                         estimateRAwVPR(*header,evpr,polargrid_Z_3d,polargrid_PHI_long,polargrid_KDP_long,flag_echo_3d,Rtype_3d,
                                        melting,vpr,azimuth,alpha,header->angle,ang,DP_SCALEING,polargrid_rain_3d);
                     }
                     else if (estimator==RZC) {
                         estimateRZCwVPR(*header,evpr,polargrid_Z_3d,flag_echo_3d,Rtype_3d,vpr,header->angle,
                                         ang,polargrid_rain_3d);
                     }
                     else if (estimator==RZCS) {
                         estimateRZCSwVPR(*header,evpr,polargrid_Z_3d,flag_echo_3d,Rtype_3d,vpr,header->angle,
                                          ang,polargrid_rain_3d);
                     }
                     else if (estimator==RZDR) {
                         estimateRZDRwVPR(*header,evpr,polargrid_Z_3d,polargrid_ZDR_3d,flag_echo_3d,Rtype_3d,melting,vpr,
                                          header->angle,ang,DP_SCALEING,polargrid_rain_3d);
                     }
                     else if (estimator==RKDP) {
                         estimateRKDPwVPR(*header,evpr,polargrid_Z_3d,polargrid_KDP_long,flag_echo_3d,Rtype_3d,melting,vpr,
                                          header->angle,ang,DP_SCALEING,polargrid_rain_3d);
                     }
                }
            }
        }

        /* final product */
        temp_azi=get1dFloatArray(MAXIMUM_AZIMUTH_SUPER);
        if (angle_number==10)
            ang=1;
        else
            ang=angle_number;

        for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
             temp_azi[azi]=azimuth[azi][ang-1];
        }

        if (angle_number>=1 && angle_number<=9) { /* angle scan */
            for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {
                 for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                      rainrate[azi][rng]=polargrid_rain_3d[azi][rng][angle_number-1];
                 }
            }
            /* Fixed grid (only for super-resolution angles) */
            get2DFixedGrid(temp_azi,bw,rainrate);
            *status=TRUE;
        }
        else if (angle_number==10) { /* CAPPI */
            getRainRateCAPPI(polargrid_rain_3d,melting,azi_lookup,cappi_nangle,
                             header->angle,cappi_h,tower_h,rainrate);
            /* Fixed grid */
            get2DFixedGrid(temp_azi,bw,rainrate);
            *status=TRUE;
        }
    }
    else
        *status=FALSE;
}
