/**************************************************************************
--> readLevel2.c

    : Reads Polarimetric Level-II file using NASA RSL

    # Function input
      - fname: level2 super-resolution file

    # Function outputs
      - azimuth: azimuth angle data for super resolution
      - polargrid_Z_3d: 3-dimensional Z data for super resolution
                        (0-255 integer values by dBZ*2+66, missing: '-99')
      - polargrid_RHO_3d: 3-dimensional Rhohv data scaled by 10000 (DP_SCALEING)
      - polargrid_ZDR_3d: 3-dimensional Zdr data scaled by 10000 (DP_SCALEING)
      - polargrid_PHI_3d: 3-dimensional Phidp data scaled by 10000 (DP_SCALEING)
      - volume_info: volume header information for both super and legacy resolution
                     (date & time, beam width, elevation angle, and radar ID)

    # Dimension information
      - azimuth: [MAXIMUM_AZIMUTH_SUPER][MAXIMUM_ANGLE_SUPER]
      - polargrid_Z_3d, polargrid_RHO_3d, polargrid_ZDR_3d, polargrid_PHI_3d
        : [MAXIMUM_AZIMUTH_SUPER][MAXIMUM_RANGE_SUPER][MAXIMUM_ANGLE_SUPER]

    Developed by Bongchul Seo, Anton Kruger and Witold F. Krajewski
    The University of Iowa

    Created : August    2009
    Modified: August    2011  Adapted resolution change (1 degree by 0.25 km)
              July      2012  Modified for partially filled rays
              February  2013  Added AP algorithm using Dual Pol. RHO variable
              September 2013  Modified AP algorithm (RHO + PHIDP)
              January   2016  Refined DualPol. quality control
              April     2016  Modified QC with RHO
              July      2017  Modified DP QC
              August    2017  Added Phidp unfolding and Kdp estimation
              September 2017  Adapted Dual-pol preprocessing from "NWS CODE"
              October   2017  Added Specific Attenuation
***************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "libphydro.h"

void readLevel2(char *fname,
                float **azimuth,
                short ***polargrid_Z_3d,
                int ***polargrid_RHO_3d,
                int ***polargrid_ZDR_3d,
                int ***polargrid_PHI_3d,
                struct volume_info* vol_info)
{
    /* Data structure: File > Radar > Volume > Sweep > Ray */
    Radar *radar;
    Volume *volume_rho, *volume_zdr, *volume_phidp;

    char path[]="/", *iname, *pname, fin[100]="", callID[5]="", temp_str[2]="";
    int i;

    /* copy file name */
    for (i=0; i<100; i++) {
        fin[i]=fname[i];
    }

    /* Separate file name from path */
    pname=strtok(fin,path);
    iname=pname;
    while (pname != NULL) {
        pname=strtok(NULL,path);
        if (pname != NULL)
            iname=pname;
    }

    /* Get radar ID */
    temp_str[0]=iname[0];
    temp_str[1]=='\0';
    for (i=0; i<4; i++) {
         if (strcmp(temp_str,"K")==0)
             callID[i]=iname[i]; // for NCDC
         else if (strcmp(temp_str,"L")==0)
             callID[i]=iname[i+7]; // for IFC (from LDM)
    }
    callID[4]='\0';

    /* Initialization */

    /* Read super-resolution file */
    radar=RSL_wsr88d_to_radar(fname,callID);
    if (radar == NULL) {
        fprintf(stderr,"Error: could not read input file: %s\n",fname);
        exit(-1);
    }

/*
	printf("\nRadar header data :: \n");
	printf("----------------------\n");
	printf("-> Month %d, Day %d, Year %d, Hour %d, Minute %d, Sec %f \n",
			radar->h.month, radar->h.day, radar->h.year, radar->h.hour, radar->h.minute, radar->h.sec);
	printf("-> Radar Type: %s \n", radar->h.radar_type);
	printf("-> Nexrad Site Name: %s \n", radar->h.name);
	printf("-> Radar Name: %s \n", radar->h.radar_name);
	printf("-> Project: %s, City: %s, State: %c%c \n", radar->h.project, radar->h.city, radar->h.state[0], radar->h.state[1]);
	printf("-> No. of Volumes (array lenght only) %d \n", radar->h.nvolumes);
	printf("-> Latitude:  %d : %d : %d\n   Longitude: %d : %d : %d, \n",
			radar->h.latd, radar->h.latm, radar->h.lats, radar->h.lond, radar->h.lonm, radar->h.lons);
*/

    if (radar->h.nvolumes <1) {
        fprintf(stderr, "Error: empty volume file: %s\n",fname);
        exit(-1);
    }

    /* Read reflectivity and DP volumes */
    getZvolume(radar,DZ_INDEX,vol_info,azimuth,polargrid_Z_3d);
    getDPvolume(radar,RH_INDEX,vol_info,polargrid_RHO_3d,DP_SCALEING);
    getDPvolume(radar,DR_INDEX,vol_info,polargrid_ZDR_3d,DP_SCALEING);
    getDPvolume(radar,PH_INDEX,vol_info,polargrid_PHI_3d,DP_SCALEING);

    RSL_free_radar(radar);
}

void getZvolume(Radar *radar,
                int index,
                struct volume_info* vol_info,
                float **azimuth,
                short ***polargrid_Z_3d)
{
    Volume *volume_zh;
    Sweep *sweep;
    Ray *ray;
    int ang, azi, rng;
    int nsuper, nlegacy, nlegacy2, max_nangle, actual_ang, max_nbin, max_nbin_SUPER, firstgate;

    /* initialize arrays with nodata value for missing rays */
    initialize2dFloatArray(azimuth,MAXIMUM_AZIMUTH_SUPER,
                           MAXIMUM_ANGLE_SUPER,NODATA);
    initialize3dShortArray(polargrid_Z_3d,MAXIMUM_AZIMUTH_SUPER,
                           MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER,NODATA);

    volume_zh=radar->v[index];
    if (volume_zh == NULL) {
        fprintf(stderr, "Error: reflectivity volume is NULL... \n");
        exit(-1);
    }

    /* Read sweep and ray data */
    nsuper=0;  // the number of super-resolution sweeps (0.5degree by 0.25 km)
    nlegacy=0; // the number of legacy-resolution sweeps (1 degree by 1 km)
    nlegacy2=0; // the number of legacy2-resolution sweeps (1 degree by 0.25 km)

    if (MAXIMUM_ANGLE_SUPER>volume_zh->h.nsweeps)
        max_nangle=volume_zh->h.nsweeps;
    else
        max_nangle=MAXIMUM_ANGLE_SUPER;

    actual_ang=0;
    for (ang=0; ang<max_nangle; ang++) {
         sweep=volume_zh->sweep[ang];
         if (sweep != NULL) {
             if (ang<MAXIMUM_ANGLE+MAXIMUM_ANGLE_SUPER) {
                 if (sweep->h.elev<=vol_info->angle[actual_ang-1])
                     continue;
                 else
                     vol_info->angle[actual_ang]=sweep->h.elev;
             }
             if (ang==0) {
                 vol_info->beam_width=(sweep->h.beam_width);
             }
             for (azi=0; azi<sweep->h.nrays; azi++) {
                  ray=sweep->ray[azi];
// if (azi==0)
//    printf("num: %d, angle: %f, nray: %d, nbin: %d, size: %d, firstgate: %d\n",sweep->h.sweep_num,sweep->h.elev,sweep->h.nrays,ray->h.nbins,ray->h.gate_size, ray->h.range_bin1);

                  max_nbin_SUPER=MAXIMUM_RANGE_SUPER;
                  if (ray->h.nbins<MAXIMUM_RANGE_SUPER)
                      max_nbin_SUPER=ray->h.nbins;

                  max_nbin=MAXIMUM_RANGE;
                  if (ray->h.nbins<MAXIMUM_RANGE)
                      max_nbin=ray->h.nbins;

                  firstgate=(ray->h.range_bin1-ray->h.gate_size/2)/ray->h.gate_size;
                  if (ray != NULL) {
                      if (ray->h.gate_size == BIN_SUPER && sweep->h.nrays>MAXIMUM_AZIMUTH) {
                          if (azi<MAXIMUM_AZIMUTH_SUPER) {
                              azimuth[azi][nsuper+nlegacy2]=ray->h.azimuth;
                              for (rng=firstgate; rng<max_nbin_SUPER; rng++) {
                                   if (ray->h.f(ray->range[rng-firstgate])==BADVAL)
                                       polargrid_Z_3d[azi][rng][nsuper+nlegacy2]=INITIAL_VALUE;
                                   else
                                       polargrid_Z_3d[azi][rng][nsuper+nlegacy2]
                                              =float_dBZ2int_dBZ(ray->h.f(ray->range[rng-firstgate]));
                              }
                          }
                          if (azi==0)
                              actual_ang=actual_ang+1;
                      }
                      else if (ray->h.gate_size == BIN_SUPER && sweep->h.nrays==MAXIMUM_AZIMUTH) {
                          if (azi<MAXIMUM_AZIMUTH) {
                              azimuth[azi][nsuper+nlegacy2]=ray->h.azimuth;
                              for (rng=firstgate; rng<max_nbin_SUPER; rng++) {
                                   if (ray->h.f(ray->range[rng-firstgate])==BADVAL)
                                       polargrid_Z_3d[azi][rng][nsuper+nlegacy2]=INITIAL_VALUE;
                                   else
                                       polargrid_Z_3d[azi][rng][nsuper+nlegacy2]
                                              =float_dBZ2int_dBZ(ray->h.f(ray->range[rng-firstgate]));
                              }
                          }
                          if (azi==0)
                              actual_ang=actual_ang+1;
                      }
                  }
             }

             if (ray->h.gate_size == BIN_SUPER && sweep->h.nrays>600) {
                 vol_info->typeangle[actual_ang-1]=SUPER;
                 nsuper=nsuper+1;
             }
             else if (ray->h.gate_size == BIN_SUPER && sweep->h.nrays==MAXIMUM_AZIMUTH) {
                 vol_info->typeangle[actual_ang-1]=LEGACY2;
                 nlegacy2=nlegacy2+1;
             }
             else if (ray->h.gate_size == BIN_LEGACY && sweep->h.nrays==MAXIMUM_AZIMUTH) {
                 vol_info->typeangle[actual_ang-1]=LEGACY;
                 nlegacy=nlegacy+1;
             }

             if (nsuper+nlegacy2>=MAXIMUM_ANGLE_SUPER)
                 continue;
             if (nlegacy>=MAXIMUM_ANGLE)
                 break;
         }
    }
    vol_info->year=radar->h.year;
    vol_info->mon=radar->h.month;
    vol_info->day=radar->h.day;
    vol_info->hour=radar->h.hour;
    vol_info->min=radar->h.minute;
    vol_info->sec=(int)(radar->h.sec);
    vol_info->nangle=nlegacy;
    vol_info->nangle_super=nsuper;
    vol_info->nangle_legacy2=nlegacy2;
    strcpy(vol_info->radar_id,radar->h.radar_name);
//printf("nsuper: %d, nlegacy: %d, nlegacy2: %d\n", nsuper,nlegacy,nlegacy2);
}

void getDPvolume(Radar *radar,
                 int index,
                 struct volume_info* vol_info,
                 int ***polargrid_DP_3d,
                 int scale)
{
    Volume *volume_dp;
    Sweep *sweep;
    Ray *ray;
    int ang, azi, rng, max_nbin, firstgate, max_nangle;

    /* initialize array */
    initialize3dIntArray(polargrid_DP_3d,MAXIMUM_AZIMUTH_SUPER,
                         MAXIMUM_RANGE_SUPER,MAXIMUM_ANGLE_SUPER,NODATA);

    volume_dp=radar->v[index];
    if (volume_dp == NULL) {
        fprintf(stderr, "Error: DP volume %d is NULL. \n",index);
        exit(-1);
    }

    if (MAXIMUM_ANGLE_SUPER>volume_dp->h.nsweeps)
        max_nangle=volume_dp->h.nsweeps;
    else
        max_nangle=MAXIMUM_ANGLE_SUPER;

    for (ang=0; ang<max_nangle; ang++) {
         sweep=volume_dp->sweep[ang];
         if (sweep != NULL) {
//printf("ang: %d, angle: %f\n",ang,sweep->h.elev);
             for (azi=0; azi<sweep->h.nrays; azi++) {
                  ray=sweep->ray[azi];

                  max_nbin=MAXIMUM_RANGE_SUPER;
                  if (ray->h.nbins<MAXIMUM_RANGE_SUPER)
                      max_nbin=ray->h.nbins;

                  firstgate=(ray->h.range_bin1-ray->h.gate_size/2)/ray->h.gate_size;
                  if (ray != NULL) {
                      for (rng=firstgate; rng<max_nbin; rng++) {
                           if (ray->h.f(ray->range[rng-firstgate])==BADVAL || ray->h.f(ray->range[rng-firstgate])==RFVAL)
                               polargrid_DP_3d[azi][rng][ang]=INITIAL_VALUE;
                           else
                               polargrid_DP_3d[azi][rng][ang]=(int)(ray->h.f(ray->range[rng-firstgate])*DP_SCALEING+0.5);
                      }
                  }
             }
         }
    }
}

