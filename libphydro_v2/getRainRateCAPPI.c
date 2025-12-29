/*************************************************************************
--> getRainRateCAPPI.c

    : Builds CAPPI from polar-based rain rate maps for multiple elevation angles
      (using log-N kernel)


Added:    ML --> December 2019 (supress a weight for an altitude above the ML

*************************************************************************/
#include <stdio.h>
#include <math.h>
#include "libphydro.h"

static void computeSuperWeight(int num_angle,
                               float *elev_angle,
                               float cappi_height,
                               float tower_height,
                               float **wgt_smooth);

void getRainRateCAPPI(float ***rain_3d,
                      short ***melting,
                      short **lookup,
                      int cappi_nangle,
                      float *elevation_angle,
                      float cappi_h,
                      float tower_h,
                      float **rr)
{
    int ang, rng, azi, cazi, count, countpair, wgt_count, melt, status;
    float **wgt_smooth, *wgt, *wgtR, wt, wgt_sum, temp_R, wgt_thresh=0.0001;

    /* define array */
    wgt_smooth=get2dFloatArray(cappi_nangle,MAXIMUM_RANGE_SUPER);
    wgt=get1dFloatArray(cappi_nangle);
    wgtR=get1dFloatArray(cappi_nangle);

    /* weight computation for smoothing */
    computeSuperWeight(cappi_nangle,elevation_angle,cappi_h,tower_h,wgt_smooth);

    /* build CAPPI */
    for (rng=0; rng<MAXIMUM_RANGE_SUPER; rng++) {
         for (azi=0; azi<MAXIMUM_AZIMUTH_SUPER; azi++) {

              /* check missing */
              count=0;
              countpair=0;
              for (ang=0; ang<cappi_nangle; ang++) {
                   if (ang==0)
                       cazi=azi;
                   else
                       cazi=lookup[azi][ang];
                   temp_R=rain_3d[cazi][rng][ang];
                   melt=melting[cazi][rng][ang];
                   wt=wgt_smooth[ang][rng];

                   //if (temp_R==NODATA || temp_R<0)
                   //    count++;
                   if (temp_R>=0)
                       count++;
                   if (temp_R>=0 && wt>wgt_thresh && melt!=MELTING_FLAG && melt!=FREEZING_FLAG)
                       countpair++;
              }

              if (count==cappi_nangle && count==countpair) { /* when there is no missing in the vertical column */
                  for (ang=0; ang<cappi_nangle; ang++) {
                       if (ang==0)
                           cazi=azi;
                       else
                           cazi=lookup[azi][ang];
                       temp_R=rain_3d[cazi][rng][ang];
                       rr[azi][rng]+=temp_R*wgt_smooth[ang][rng];
//if (rng==46 && azi==602)
//    printf("ang: %d, temp_R: %f, wgt: %f, R: %f\n",ang,temp_R,wgt_smooth[ang][rng],rr[azi][rng]);
                  }
              }
              else { /* adjust weights */
                  wgt_sum=0.;
                  wgt_count=0;
                  for (ang=0; ang<cappi_nangle; ang++) {
                       if (ang==0)
                           cazi=azi;
                       else
                           cazi=lookup[azi][ang];
                       temp_R=rain_3d[cazi][rng][ang];
                       melt=melting[cazi][rng][ang];

                       if (countpair>0) {
                           if (temp_R>=0 && melt!=MELTING_FLAG && melt!=FREEZING_FLAG) {
                               if (wgt_smooth[ang][rng]>wgt_thresh) {
                                   wgt[wgt_count]=wgt_smooth[ang][rng];
                                   wgt_sum+=wgt_smooth[ang][rng];
                                   wgtR[wgt_count]=temp_R;
                                   wgt_count++;
                               }
                           }
                       }
                       else {
                           if (temp_R>=0) {
                               if (wgt_smooth[ang][rng]>wgt_thresh) {
                                   wgt[wgt_count]=wgt_smooth[ang][rng];
                                   wgt_sum+=wgt_smooth[ang][rng];
                                   wgtR[wgt_count]=temp_R;
                                   wgt_count++;
                               }
                           }
                       }
                  }

                  if (wgt_sum>0) {
                      for (ang=0; ang<wgt_count; ang++) {
                           rr[azi][rng]+=wgtR[ang]*wgt[ang]/wgt_sum;
                      }
                  }
              }
         }
    }

    status=free2dFloatArray(wgt_smooth,cappi_nangle);
    //free(wgt);
    //free(wgtR);
}

static void computeSuperWeight(int num_angle,
                               float *elev_angle,
                               float cappi_height,
                               float tower_height,
                               float **wgt_smooth)
{
    int ang, bin, status;
    double beam_height, real_height, sum_wgt,
           wgt_thresh=0.0001, **pdf;
    float mu, sigma;

    /* Dynamic allocation */
    pdf=get2dDoubleArray(num_angle,MAXIMUM_RANGE_SUPER);
    tower_height=tower_height/1000;

    /* Gaussian (normal) distribution */
/*
    for (ang=0; ang<num_angle; ang++) {
         for (bin=0; bin<MAXIMUM_RANGE_SUPER; bin++) {
              beam_height=computeBeamHeight(elev_angle[ang],(float)bin/4+0.25/2);
              real_height=beam_height+tower_height;
              pdf[ang][bin]=1/(sqrt(2*M_PI)*gaussian_width)*
                    exp(-(real_height-cappi_height)*(real_height-cappi_height)/
                    (2*gaussian_width*gaussian_width));
         }
    }
*/

    /* Parameter estimation */
    /* Log-normal distribution
       - CAPPI height 1.5 km with mu (0.50) and sigma (0.31)
       - CAPPI height 2.0 km with mu (0.75) and sigma (0.24) */

    if (cappi_height<1.5)
        cappi_height=1.5;

    if (cappi_height==1.5) {
        mu=(float)0.5;
        sigma=(float)0.31;
    }
    else {
        sigma=(float)0.31*(1-cappi_height/11);
        mu=(float)(sigma*sigma+log(cappi_height));
    }

    for (ang=0; ang<num_angle; ang++) {
         for (bin=0; bin<MAXIMUM_RANGE_SUPER; bin++) {
              beam_height=computeBeamHeight(elev_angle[ang],(float)bin/4+0.25/2);
              real_height=beam_height+tower_height;
              pdf[ang][bin]=1/(real_height*sigma*sqrt(2*M_PI))*
                    exp(-(log(real_height)-mu)*(log(real_height)-mu)/
                    (2*sigma*sigma));
         }
    }

    /* Normalization of weight */
    for (bin=0; bin<MAXIMUM_RANGE_SUPER; bin++) {
         sum_wgt=INITIAL_VALUE;
         for (ang=0; ang<num_angle; ang++) {
             if (pdf[ang][bin]>wgt_thresh)
                 sum_wgt=sum_wgt+pdf[ang][bin];
         }

         // Calculate weight
         if (sum_wgt==0) {
             if (bin<20) // within 5 km
                 wgt_smooth[num_angle-1][bin]=1;
             else if (bin>400)
                 wgt_smooth[0][bin]=1;
         }
         else {
             for (ang=0; ang<num_angle; ang++) {
                  if (pdf[ang][bin]>wgt_thresh)
                      wgt_smooth[ang][bin]=(float)(pdf[ang][bin]/sum_wgt);
             }
         }
    }
    status=free2dDoubleArray(pdf,num_angle);
}
