/*************************************************************************
--> getMeltinglayer.c

    : Defines melting layer altitude (-5 to 5 degrees) and above

    (real-time)     -> '-d' argument gets an exact name of Tsounding data
    (retrospective) ->                    a path whtere the Tsounding data exist

*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "libphydro.h"

void getMeltinglayer(struct volume_info header,
                     char *nwppath,
                     float **azimuth,
                     short ***melting,
                     float ***T_top,
                     float ***T_bottom,
                     int *flag_nwp)
{
    FILE *ft;
    int i, j, ang, nazi, ngrid, azi, iazi, rng, max_nangle, ival, status, **lookup, *lookup_ray;
    float **T, **h, *T_ray, avgT, fval;
    char date_string[12], year_string[5], mon_string[3], day_string[3], hour_string[3],
         zero_string[2]="0", temp_string[3]="", id_string[5]="", name_Tprofile[50]="", name_lookup[50]="", temp_lookup[50]="";

    strcpy(temp_lookup,nwppath);
    if (strlen(nwppath)<20) { // in case where only the data path is provide
        sprintf(year_string, "%d",header.year);
        sprintf(mon_string, "%d",header.mon);
        if (header.mon<10){
            strcpy(temp_string,mon_string);
            strcpy(mon_string,zero_string);
            strcat(mon_string,temp_string);
        }
        sprintf(day_string, "%d", header.day);
        if (header.day<10){
            strcpy(temp_string,day_string);
            strcpy(day_string,zero_string);
            strcat(day_string,temp_string);
        }
        strcpy(date_string,year_string);
        strcat(date_string,mon_string);
        strcat(date_string,day_string);
        strcpy(name_Tprofile,nwppath);
        strcat(name_Tprofile,"/Tprofile_");
/*        strcpy(temp_string,header.radar_id); -> fixed: temp_strng>3 */
        strcpy(id_string,header.radar_id);
        strcat(name_Tprofile,id_string);
        strcat(name_Tprofile,"_");
        strcat(name_Tprofile,date_string);
        strcat(name_Tprofile,"_");
        
        sprintf(hour_string, "%d",header.hour);
        if (header.hour<10){
            strcpy(temp_string,hour_string);
            strcpy(hour_string,zero_string);
            strcat(hour_string,temp_string);
        }
        strcat(name_Tprofile,hour_string);
        strcat(name_Tprofile,"00.dat");
    }
    else { // in case where the full data name is provided
        strcpy(name_Tprofile,nwppath);
    }
    /* read NWP temperature profile and altitude ------------------------*/
    if ((ft=fopen(name_Tprofile,"r"))==NULL) {
         fprintf(stderr, "Error.....: no NWP Temperature sounding, %s\n",name_Tprofile);
         *flag_nwp=FALSE;
         return;
    }
    fscanf(ft,"%d",&ngrid);

    /* Arrays */
    T=get2dFloatArray(ngrid,NTEMP);
    h=get2dFloatArray(ngrid,NTEMP);

    for (i=0; i<ngrid; i++) {
         for (j=0; j<NTEMP; j++) {
              fscanf(ft,"%f",&fval);
              T[i][j]=fval;
         }
    }
    for (i=0; i<ngrid; i++) {
         for (j=0; j<NTEMP; j++) {
              fscanf(ft,"%f",&fval);
              h[i][j]=fval;
         }
    }
    fclose(ft);

    /* read lookup table (NWP grids with radar polar grids) ---------------------*/
    if (strlen(nwppath)<20) {
        strcpy(name_lookup,nwppath);
    }
    else {
        for (i=0; i<strlen(nwppath)-31-1; i++) {  // 31: the number of strings in the data file
             name_lookup[i]=temp_lookup[i];
        }
    }
    strcat(name_lookup,"/Lookup_Tprofile-");
    strcat(name_lookup,header.radar_id);
    strcat(name_lookup,".dat");

    /* Array */
    lookup=get2dIntArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    if ((ft=fopen(name_lookup,"r"))==NULL) {
         fprintf(stderr, "Error.....: no NWP lookup file, %s\n",name_lookup);
    }

    for (i=0; i<MAXIMUM_AZIMUTH_SUPER; i++) {
         for (j=0; j<MAXIMUM_RANGE_SUPER; j++) {
              fscanf(ft,"%d",&ival);
              lookup[i][j]=ival;
         }
    }
    fclose(ft);

    /* find melting layer height and corresponding range ------------------------*/
    /* Arrays */
    lookup_ray=get1dIntArray(MAXIMUM_RANGE_SUPER);
    T_ray=get1dFloatArray(MAXIMUM_RANGE_SUPER);

    max_nangle=header.nangle+header.nangle_super+header.nangle_legacy2;
    for (ang=0; ang<max_nangle; ang++) {
         //if (ang<header.nangle_super)
         if (header.typeangle[ang]==SUPER)
             nazi=MAXIMUM_AZIMUTH_SUPER;
         else
             nazi=MAXIMUM_AZIMUTH;

         for (azi=0; azi<nazi; azi++) {
              if (azimuth[azi][ang]!=NODATA) { // to avoid segfault from incomplete scans (added June 2020)
                  iazi=(int)(2*azimuth[azi][ang]);
                  for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                       lookup_ray[rng]=lookup[iazi][rng];
                  }
              }

              /* beam top temperature */
              initialize1dFloatArray(T_ray,MAXIMUM_RANGE_SUPER,NODATA);
              computeTemperature(header.angle[ang]+header.beam_width/2,lookup_ray,T,h,T_ray);
              for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                   T_top[azi][rng][ang]=T_ray[rng];
              }

              /* beam bottom temperature */
              initialize1dFloatArray(T_ray,MAXIMUM_RANGE_SUPER,NODATA);
              computeTemperature(header.angle[ang]-header.beam_width/2,lookup_ray,T,h,T_ray);
              for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
                   if (T_ray[rng]!=NODATA) { // due to missing rays
                       T_bottom[azi][rng][ang]=T_ray[rng];
                       if (T_bottom[azi][rng][ang]<=MELTING_THRESHOLD_ABOVE)
                           melting[azi][rng][ang]=FREEZING_FLAG;
                       else if (T_bottom[azi][rng][ang]>MELTING_THRESHOLD_ABOVE && T_bottom[azi][rng][ang]<MELTING_THRESHOLD_BELOW) {
                           if (T_top[azi][rng][ang]>MELTING_THRESHOLD_ABOVE)
                               melting[azi][rng][ang]=MELTING_FLAG;
                           else {
                               avgT=0.5*(T_top[azi][rng][ang]+T_bottom[azi][rng][ang]);
                               if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW)
                                   melting[azi][rng][ang]=MELTING_FLAG;
                               else if (avgT<=MELTING_THRESHOLD_ABOVE)
                                   melting[azi][rng][ang]=FREEZING_FLAG;
                           }
                       }
                   }
                   /*
                   else if (T_bottom[azi][rng][ang]>=MELTING_THRESHOLD_BELOW) {
                       if (T_top[azi][rng][ang]>=MELTING_THRESHOLD_BELOW)
                           continue;
                       else {
                           avgT=0.5*(T_top[azi][rng][ang]+T_bottom[azi][rng][ang]);
                           if (avgT>MELTING_THRESHOLD_ABOVE && avgT<MELTING_THRESHOLD_BELOW)
                               melting[azi][rng][ang]=MELTING_FLAG;
                           else if (avgT<=MELTING_THRESHOLD_ABOVE)
                               melting[azi][rng][ang]=FREEZING_FLAG;
                       }
                   }
                   */
              }
         }
    }
    free(lookup_ray);
    free(T_ray);
    status=free2dFloatArray(T,ngrid);
    status=free2dFloatArray(h,ngrid);
    status=free2dIntArray(lookup,MAXIMUM_AZIMUTH_SUPER);
}

void computeTemperature(float angle,
                        int *lookup_ray,
                        float **T,
                        float **h,
                        float *T_ray)
{
    int i, j, ilookup;
    float bheight, h1, h2, T1, T2;

    /* bin bottom height caluation and temperature (linear) interpolation */
    for (i=0; i<MAXIMUM_RANGE_SUPER; i++) {
         if (lookup_ray[i]!=NODATA) { // missing rays
             bheight=computeBeamHeight((double)angle,i*0.25+0.125);
             ilookup=lookup_ray[i]-1;

             /* find corresponding altitude */
             if (bheight>h[ilookup][NTEMP-1]) {
                 T_ray[i]=T[ilookup][NTEMP-1]-273.15;
                 continue;
             }
             else if (bheight<=h[ilookup][0]) {
                 T_ray[i]=T[ilookup][0]-273.15;
                 continue;
             }
             else {
                 for (j=0; j<NTEMP-1; j++) {
                      if (bheight>h[ilookup][j] && bheight<=h[ilookup][j+1]) {
                          h1=h[ilookup][j];
                          h2=h[ilookup][j+1];
                          T1=T[ilookup][j];
                          T2=T[ilookup][j+1];
                          T_ray[i]=(T2-T1)*(bheight-h1)/(h2-h1)+T1-273.15;
                          break;
                      }
                 }
             }
         }
    }
}
