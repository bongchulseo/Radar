#include "rsl.h"

/* For radar processing */
#define TRUE            1
#define FALSE           0
#define YES             1
#define NO              0
#define UNKNOWN        -1         /* unknown field in sweep header */
#ifndef BADVAL
#define BADVAL          320       /* data below SNR threshold */
#define RFVAL          (BADVAL-1) /* meaning unclear. A.K. 10/18/94 */
#define APFLAG         (BADVAL-2) /* AP mask: keep pixel */
#define APFLAG2        (BADVAL-3) /* AP mask: remove pixel */
#endif

#define BUFFERSIZE 8192
#define MAXIMUM_AZIMUTH 360
#define MAXIMUM_RANGE   230 // We will use the data within 230km from the radar
#define MAXIMUM_RAY     400
#define MAXIMUM_RANGE_HALF 115 // for 2km bin
#define MAXIMUM_ANGLE    6
#define NONE 0
#define NODATA -99
#define INITIAL_VALUE 0
#define MISSING_FLAG -1
#define REFRACTIVE_INDEX 1.21
#define EARTH_RADIUS 6371.0
#define lowZLimit_val 0.0
#define RATE_SCALING 10.0
#define RAIN_SCALING 100.0
#define DP_SCALEING 10000
#define ZHAIL   53.
#define RAD (double) (M_PI/180)
#define BIN_SUPER   250
#define BIN_LEGACY 1000
#define KM 1000.
#define HOURINSEC   3600
#define HALFHOURINSEC   1800
#define NTEMP 37

/* Melting layer */
#define MELTING_THRESHOLD_BELOW 5.
#define MELTING_THRESHOLD_ABOVE -5.
#define MELTING_FLAG 1
#define FREEZING_FLAG 2

/* Rain Type */
#define CONVECTIVE  3
#define STRATIFORM  2
#define WINDTURBINE 1
#define NOTHING 2

/* For super resolution data */
#define MAXIMUM_AZIMUTH_SUPER 720
#define MAXIMUM_RANGE_SUPER   920
#define MAXIMUM_RANGE_SUPER_HALF 460
#define MAXIMUM_ANGLE_SUPER   15

/* For multiple radar data merging */
#define MAXIMUM_LATITUDE 900
#define MAXIMUM_LONGITUDE 1500
#define MAXIMUM_RADAR 30

/* Radar data type */
#define REFLECTIVITY 0
#define RATE 1
#define ACCUMULATION 2
#define ELSE 3

/* Merging type */
#define DATA_BASED 0
#define PRODUCT_BASED 1

/* DP processing */
#define DONTKNOW -1
#define RAIN 0
#define CLUTTER 1
#define NORAIN 2
#define HAIL 3
#define DRYSNOW 4  // dry snow within melting layer
#define WETSNOW 5
#define SNOWICE 6 // dry snow above melting layer
#define ICE 7
#define MET_ECHO 1
#define NO_MET_ECHO 0

/* QPE method */
#define RZC      1
#define RZCS     2
#define RZDR     3
#define RKDP     4
#define RA       5

struct volume_info {
    int year;
    int mon;
    int day;
    int hour;
    int min;
    int sec;
    float beam_width;
    int nangle;
    int nangle_legacy2;
    int nangle_super; /* super resolution */
    float angle[MAXIMUM_ANGLE+MAXIMUM_ANGLE_SUPER];
    int typeangle[MAXIMUM_ANGLE+MAXIMUM_ANGLE_SUPER];
    char radar_id[5];
};
/* type angle */
#define LEGACY 0
#define SUPER 1
#define LEGACY2 2

struct bounding {
    int azi1;
    int azi2;
    int bin1;
    int bin2;
};
struct range_parameter {
    int start_range;      /* VPR analysis range */
    int end_range;        /* VPR analysis range */
    int sector_width;     /* sector width to estimate local VPR */
    int sector_move;      /* azimuthal movement to estimate local VPR */
    int hmax;             /* maximum height (km) to estimate VPR */
    double dh;            /* vertical resolution (km)*/
};
struct asciiheader { /* Ascii grid file header information */
    int dt;
    int number_column;
    int number_row;
    int corner_x;
    int corner_y;
    float xllcorner;
    float yllcorner;
    float cellsize;
    float nodata_value;
    char radar_id[5];
    char month[4];
    int year;
    int mon;
    int day;
    int hour;
    int min;
    int sec;
};
struct polarheader { /* Polar grid file header information */
    int year;
    int mon;
    int day;
    int hour;
    int min;
    int sec;
};
struct mergedradar { /* radar to be merged */
    char inf[1500];  /* input file name to be merged */
    char dist[100];   /* input distance file name */
    int  dt;         /* integration time for accumulation */
};

void getAngleScan(FILE *fp,
                  int angle_number,
                  int range_flag,
                  int norain,
                  struct volume_info* header,
                  struct bounding box,
                  struct range_parameter evpr,
                  float (*dBZ_product)[MAXIMUM_RANGE],
                  int *status);
void getCAPPI(FILE *fp,
              int range_flag,
              float cappi_height,
              float tower_height,
              float gaussian_width,
              int norain,
              struct volume_info* header,
              struct bounding box,
              struct range_parameter evpr,
              float (*dBZ_product)[MAXIMUM_RANGE],
              int *status);
void getSCAPPI(FILE *fp,
              int range_flag,
              float cappi_height,
              float tower_height,
              float site_elevation,
              float gaussian_width,
              int norain,
              struct volume_info* header,
              struct bounding box,
              struct range_parameter evpr,
              float (*dBZ_product)[MAXIMUM_RANGE],
              int *status);
void getHybridScan(FILE *fp,
                   int range_flag,
                   int norain,
                   struct volume_info* header,
                   struct bounding box,
                   struct range_parameter evpr,
                   float (*dBZ_product)[MAXIMUM_RANGE],
                   int *status);
void readRLE(FILE *fdata,
             float **azimuth,
             short ***polargrid_Z_3d,
             struct volume_info* volume);
void make3Dfixedgrid(float **azimuth_volume,
                     short ***polargrid_Z_3d,
                     float beam_width,
                     int number_angle,
                     float ***fixedgrid_Z_3d);
/*void make2DFixedGrid(float *razimuth,
                     float BWidth,
                     short **stored_Z,
                     short **fixedgrid_Z_2d);
*/
void correctRangeEffect(struct volume_info header,
                        struct bounding box,
                        struct range_parameter evpr,
                        float ***fixedgrid_Z_3d);
void buildCAPPI(float ***volume_data,
                float *elevation_angle,
                int total_num_angle,
                float cappi_height,
                float tower_height,
                float gaussian_width,
                int norain,
                float (*dBZ_product)[MAXIMUM_RANGE]);
void removeClutter(int norain,
                   float (*dBZ_product)[MAXIMUM_RANGE]);
void buildSCAPPI(float ***volume_data,
                float *elevation_angle,
                int total_num_angle,
                float cappi_height,
                float tower_height,
                float site_elevation,
                float gaussian_width,
                int norain,
                float (*dBZ_product)[MAXIMUM_RANGE]);
void buildAngleScan(float ***fixedgrid_Z_3d,
                    int elevation,
                    struct bounding box,
                    int norain,
                    float (*dBZ_product)[MAXIMUM_RANGE]);
void buildHybridScan(float ***fixedgrid_Z_3d,
                     struct bounding box,
                     float *elevation_angle,
                     float width_beam,
                     int norain,
                     float (*dBZ_product)[MAXIMUM_RANGE]);
double computeBeamHeight(double angle, double slant_range);

/* Dynamic allocation */
char* get1dCharArray(int dim_x);
short* get1dShortArray(int dim_x);
int* get1dIntArray(int dim_x);
float* get1dFloatArray(int dim_x);
double* get1dDoubleArray(int dim_x);
short** get2dShortArray(int dim_x,int dim_y);
int** get2dIntArray(int dim_x,int dim_y);
float** get2dFloatArray(int dim_x,int dim_y);
double** get2dDoubleArray(int dim_x,int dim_y);
short*** get3dShortArray(int dim_x,int dim_y,int dim_z);
int*** get3dIntArray(int dim_x,int dim_y,int dim_z);
float*** get3dFloatArray(int dim_x,int dim_y,int dim_z);
double*** get3dDoubleArray(int dim_x,int dim_y,int dim_z);
int free1dArray(void *array);
int free2dArray(void **array,int dim_x);
int free2dIntArray(int **array,int dim_x);
int free2dFloatArray(float **array,int dim_x);
int free2dDoubleArray(double **array,int dim_x);
int free3dArray(void ***array,int dim_x,int dim_y);
int free3dIntArray(int ***array,int dim_x,int dim_y);
int free3dFloatArray(float ***array,int dim_x,int dim_y);
int free3dDoubleArray(double ***array,int dim_x,int dim_y);

/* Initialize Arrays */
void initialize3dShortArray(short ***array,int dim_x,int dim_y,int dim_z,short data_value);
void initialize2dShortArray(short **array,int dim_x,int dim_y,short data_value);
void initialize1dShortArray(short *array,int dim_x,short data_value);
void initialize3dIntArray(int ***array,int dim_x,int dim_y,int dim_z,int data_value);
void initialize2dIntArray(int **array,int dim_x,int dim_y,int data_value);
void initialize1dIntArray(int *array,int dim_x,int data_value);
void initialize3dFloatArray(float ***array,int dim_x,int dim_y,int dim_z,float data_value);
void initialize2dFloatArray(float **array,int dim_x,int dim_y,float data_value);
void initialize1dFloatArray(float *array,int dim_x,float data_value);
void initialize3dDoubleArray(double ***array,int dim_x,int dim_y,int dim_z,double data_value);
void initialize2dDoubleArray(double **array,int dim_x,int dim_y,double data_value);
void initialize1dDoubleArray(double *array,int dim_x,double data_value);

/* Data-based merging */
void mergeData(short **Z_merged,
               short ***Z_grid,
               short ***dist_grid,
               struct asciiheader header,
               int number_radar,
               float para,
               float gamma);

/* Product-based merging */
void mergeProduct(short **R_merged,
                  short ***R_grid,
                  short ***dist_grid,
                  struct asciiheader header,
                  int number_radar);

void makeMergedName(char *out_name,
                    char *inf,
                    int number_radar,
                    int data_type);
void writeMergedOutput(FILE *fout,
                       short **merged,
                       struct asciiheader header,
                       struct mergedradar *nexrad,
                       int number_radar,
                       int data_type);

/* Read ascii grid file */
void readAsciiHeader(FILE *fin,
                     struct asciiheader* header);
int readDistanceHeader(FILE *fin,
                       struct asciiheader header);
void readAsciiData(FILE *fin,
                   short **data,
                   int number_rows,
                   int number_columns,
                   int void_value,
                   int data_type);
void readMergedHeader(FILE *fin,
                      char *fname,
                      int merging_type,
                      struct asciiheader* header);
void store3dAsciiData(FILE *fin,
                      short ***data,
                      int number_rows,
                      int number_columns,
                      int void_value,
                      int data_type,
                      int index);

/* Read polar file */
void readPolarHeader(FILE *fin,
                     char *fname,
                     int *dt,
                     struct polarheader* header);
void readPolarData(FILE *fin,
                   short **data,
                   int number_rows,
                   int number_columns,
                   int void_value,
                   int data_type);
void store3dPolarData(FILE *fin,
                      short ***data,
                      int number_rows,
                      int number_columns,
                      int void_value,
                      int data_type,
                      int index);

/* Super resolution */
void readLevel2(char *fname,
                float **azimuth_super,
                short ***polargrid_Z_3d_super,
                int ***polargrid_RHO_3d_super,
                int ***polargrid_ZDR_3d_super,
                int ***polargrid_PHI_3d_super,
                struct volume_info* vol_info);
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
             int *status);
void getRainRateCAPPI(float ***rain_3d,
                      short ***melting,
                      short **lookup,
                      int cappi_nangle,
                      float *elevation_angle,
                      float cappi_h,
                      float tower_h,
                      float **rr);








void makeSuper3DFixedGrid(float **azimuth,
                          float **azimuth_super,
                          short ***polargrid_Z_3d,
                          short ***polargrid_Z_3d_super,
                          float beam_width,
                          int number_angle,
                          int number_angle_super,
                          int number_angle_legacy2,
                          float ***fixedgrid_Z_3d);
void makeSuper2DFixedGrid(float *razimuth,
                     float BWidth,
                     short **stored_Z,
                     short **fixedgrid_Z_2d);
void getDualAngleScan(char *fname,
                      int angle_number,
                      int ap_flag,
                      int range_flag,
                      int norain,
                      struct volume_info* header,
                      struct bounding box,
                      struct range_parameter evpr,
                      float (*dBZ_product)[MAXIMUM_RANGE_SUPER],
                      int *status);
void getDualCAPPI(char *fname,
                  int ag_flag,
                  int range_flag,
                  float cappi_height,
                  float tower_height,
                  float gaussian_width,
                  int norain,
                  struct volume_info* header,
                  struct bounding box,
                  struct range_parameter evpr,
                  float (*dBZ_product)[MAXIMUM_RANGE_SUPER],
                  int *status);
void getSuperCAPPI2(char *fname,
                    int ap_flag,
                    int range_flag,
                    float cappi_height,
                    float tower_height,
                    float gaussian_width,
                    int norain,
                    struct volume_info* header,
                    struct bounding box,
                    struct range_parameter evpr,
                    float (*dBZ_product)[MAXIMUM_RANGE_SUPER],
                    int *status);
void buildSuperAngleScan(float ***fixedgrid_Z_3d,
                         int elevation,
                         struct bounding box,
                         int norain,
                         float (*dBZ_product)[MAXIMUM_RANGE_SUPER]);
void buildSuperCAPPI(float ***volume_data,
                     float *elevation_angle,
                     int total_num_angle,
                     int total_num_angle_super,
                     int total_num_angle_legacy2,
                     float cappi_height,
                     float tower_height,
                     float gaussian_width,
                     int norain,
                     float (*dBZ_product)[MAXIMUM_RANGE_SUPER]);
void removeSuperClutter(int norain,
                        float (*dBZ_product)[MAXIMUM_RANGE_SUPER]);
void correctSuperAP(struct volume_info header,
                    float ***fixedgrid_Z_3d);

void recombineSuperZ(float **azimuth,
                     float **azimuth_super,
                     short ***polargrid_Z_3d,
                     short ***polargrid_Z_3d_super,
                     float beam_width,
                     int number_angle,
                     int number_angle_super,
                     float ***fixedgrid_Z_3d);
void getRecombinedDualCAPPI(char *fname,
                            int ap_flag,
                            int range_flag,
                            float cappi_height,
                            float tower_height,
                            float gaussian_width,
                            int norain,
                            struct volume_info* header,
                            struct bounding box,
                            struct range_parameter evpr,
                            float (*dBZ_product)[MAXIMUM_RANGE],
                            int *status);
void getRecombinedDualAngleScan(char *fname,
                                int angle_number,
                                int ap_flag,
                                int range_flag,
                                int norain,
                                struct volume_info* header,
                                struct bounding box,
                                struct range_parameter evpr,
                                float (*dBZ_product)[MAXIMUM_RANGE],
                                int *status);
void correctRecombineAP(struct volume_info header,
                        float ***fixedgrid_Z_3d);
void getDualHybridScan(char *fname,
                       int ap_flag,
                       int range_flag,
                       int norain,
                       struct volume_info* header,
                       struct bounding box,
                       struct range_parameter evpr,
                       float (*dBZ_product)[MAXIMUM_RANGE],
                       int *status);
void buildHybridScan(float ***fixedgrid_Z_3d,
                     struct bounding box,
                     float *elevation_angle,
                     float width_beam,
                     int norain,
                     float (*RdBZ_product)[MAXIMUM_RANGE]);

/* Super resolution - merging */
void getDualSCAPPI(char *fname,
                   int ap_flag,
                   int range_flag,
                   float cappi_height,
                   float tower_height,
                   float site_elevation,
                   float gaussian_width,
                   int norain,
                   struct volume_info* header,
                   struct bounding box,
                   struct range_parameter evpr,
                   float (*dBZ_product)[MAXIMUM_RANGE_SUPER],
                   int *status);
void buildSuperSCAPPI(float ***volume_data,
                      float *elevation_angle,
                      int total_num_angle,
                      int total_num_angle_super,
                      int total_num_angle_legacy2,
                      float cappi_height,
                      float tower_height,
                      float site_elevation,
                      float gaussian_width,
                      int norain,
                      float (*dBZ_product)[MAXIMUM_RANGE_SUPER]);

/* Reflectivity comparison - super */
void extractSuperAPmask(struct volume_info header,
                        float ***fixedgrid_Z_3d,
                        int **APmask);

/* Dual Polarization */
void getDualPolScan(char *fname,
                    int product_type,
                    struct volume_info* header,
                    struct bounding box,
                    float (*dBZ_product)[MAXIMUM_RANGE_SUPER],
                    int *status);
void readDualPol(char *fname,
                 int product_type,
                 float **azimuth,
                 float **azimuth_super,
                 short ***polargrid_Z_3d,
                 short ***polargrid_Z_3d_super,
                 struct volume_info* vol_info);
void makeDualPol3DFixedGrid(int product_type,
                            float **azimuth,
                            float **azimuth_super,
                            short ***polargrid_Z_3d,
                            short ***polargrid_Z_3d_super,
                            float beam_width,
                            int number_angle,
                            int number_angle_super,
                            int number_angle_legacy2,
                            float ***fixedgrid_Z_3d);
void buildDualPolScan(float ***fixedgrid_Z_3d,
                      int elevation,
                      struct bounding box,
                      float (*dBZ_product)[MAXIMUM_RANGE_SUPER]);






/* Rainfall type classification */
void computeVIL(struct volume_info header,
                short ***z,
                int ***rho,
                int ***phi,
                float ***T_top,
                float ***T_bottom,
                float **azimuth,
                short **lookup,
                float **VIL,
                short **wf);
void classifyRainType(short ***z,
                 int ***rho,
                 float **vil,
                 short **wf,
                 short **type);
void identifyWF(short ***z,
                int ***rho,
                float **z_rf,
                short **wf,
                float zrain,
                float znorain);
int checkHorizontalContiuity(short ***z,
                             int ***rho,
                             int ***phi,
                             int azi,
                             int rng);
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
                           float z_thresh);
int checkSurroundingZvalues(short ***z,
                            int azi,
                            int rng,
                            int ang,
                            int dw,
                            float z_thresh);
float getinterpZ(float *t,
                  float *z,
                  int count,
                  float tx);
int defineConvective(short **type,
                     float **vil,
                     int azi,
                     int rng,
                     int dw);
int defineWF(short **type,
             short ***z,
             int ***rho,
             float **z_rf,
             int azi,
             int rng,
             int dw,
             float zrain,
             float znorain);
void updateIndicator(short **dummy,
                     short **type,
                     short **new,
                     int dw,
                     int const_convective,
                     int const_stratiform);
float getConvection(short **type,
                    int azi,
                    int rng);
void checkSurroundingTurbine(short **type,
                             short **new);

/* QC*/
void startQualityControl(struct volume_info header,
                         short ***z,
                         int ***rho,
                         int ***zdr,
                         int ***phi,
                         int ***echo,
                         int angle);
float getCoverage(short ***z,
                  int ***rho,
                  int azi,
                  int rng,
                  int ang,
                  int nazi);
int getMean(int ***polargrid,
            int azi,
            int rng,
            int ang,
            int nazi);
int getStd(int ***polargrid,
           int azi,
           int rng,
           int ang,
           int nazi);
void smoothClutter(struct volume_info header,
                   int ***echo,
                   short ***z,
                   int angle);
int defineClutter(int ***echo,
                  short zint,
                  int azi,
                  int rng,
                  int ang);
void maskObservables(struct volume_info header,
                     int ***echo,
                     short ***z,
                     int ***rho,
                     int ***zdr,
                     int ***phi,
                     int angle);
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
                int scale);
void fillGaps(struct volume_info header,
              int ***phi,
              int ***flag,
              int w,
              int angle,
              int scale);
int countMetGroups(int *flag,
                   int *ibegin,
                   int *iend);
void interpolatePhi(float *phi_ray,
                    float *new_ray,
                    int *class_ray,
                    int *flag_ray,
                    int w,
                    int ng,
                    int *ibegin,
                    int *iend);

/* vertical profile */
void getvaziLookup(struct volume_info header,
                   int nangle,
                   float **azimuth,
                   short **lookup);
void buildVPR(struct volume_info header,
              struct bounding box,
              struct range_parameter evpr,
              float **azimuth,
              short **lookup,
              short **rtype,
              short ***z,
              float **vp,
              int *flag_vp);
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
               int *flag_vp);
void assignValue2RHI(struct volume_info header,
                     struct range_parameter evpr,
                     float **azimuth,
                     double dval,
                     int azi,
                     int bin,
                     int ns,
                     double **rhi,
                     int **rhicount);
void assignValue2RHIvar(struct volume_info header,
                        struct range_parameter evpr,
                        float **azimuth,
                        double dval,
                        int azi,
                        int bin,
                        int ns,
                        double **rhi,
                        int **rhicount,
                        float **varrhi);

/* melting layer */
void getMeltinglayer(struct volume_info header,
                     char *nwppath,
                     float **azimuth,
                     short ***melting,
                     float ***T_top,
                     float ***T_bottom,
                     int *flag_nwp);
void computeTemperature(float angle,
                        int *lookup_ray,
                        float **T,
                        float **h,
                        float *T_ray);




/* DP variable processing */
void getZvolume(Radar *radar,
                int index,
                struct volume_info* vol_info,
                float **azimuth_super,
                short ***polargrid_Z_3d_super);
void getDPvolume(Radar *radar,
                 int index,
                 struct volume_info* vol_info,
                 int ***polargrid_DP_3d_super,
                 int scale);
void unfoldPhidp(struct volume_info header,
                 int ***Phidp_3d,
                 int ***Rho_3d,
                 int angle,
                 int scale);
int unfoldRay(float *Phidp,
              float *Rho,
              float Rho_thresh);
float Standard_deviation(int num,
                         float *data);
void averageZ(struct volume_info header,
              short ***z,
              int w,
              int angle);
void averageZ_ray(float *ray,
                 int w);
void averageDPvariables(struct volume_info header,
                        int ***DP,
                        int w,
                        int angle,
                        int scale);
void average_ray(float *ray,
                 int w);
void medianDPvariables(struct volume_info header,
                       int ***DP,
                       int w,
                       int angle,
                       int scale);
void medianfilter_ray(float *ray,
                      int w);
float medianfilter(float *arr,
                   int n);
void computeKdp(struct volume_info header,
                int ***phi,
                int ***rho,
                short ***z,
                int ***kdp,
                int w,
                int angle,
                int scale);
void adjustKdp(struct volume_info header,
               int angle,
               short ***z,
               int ***shortkdp,
               int ***longkdp);
void getkdp_ray(float *phi,
                float *rho,
                float *z,
                int w,
                float *kdp);
float Calculate_kdp(float *phi,
                    int m,
                    float g_size);
float Calculate_lls_kdp(float *phi,
                        int m,
                        float g_size);
float computeAlpha(struct volume_info header,
                   short ***z,
                   int ***zdr,
                   int ***rho,
                   short ***melting,
                   int scale);
void sortAscending(float *array,
                   int n);
float getLeastSquareSlope(float *x,
                          float *y,
                          int n);
float getWeightedLeastSquareSlope(float *x,
                                  float *y,
                                  int *nofpair,
                                  int n);

/* Rainfall estimation */
void estimateRZdefault(struct volume_info header,
                       short ***z,
                       int ***echo,
                       float ***rr,
                       int angle);
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
                int scale);
void estimateRR_RAT(struct volume_info header,
                    short ***z,
                    int ***phi,
                    int ***kdp,
                    int ***echo,
                    short **rtype,
                    float **azimuth,
                    float ***rr,
                    float ***airT,
                    float a,
                    int angle,
                    int scale);
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
                    float ***rr);
void getRayRA(float *z,
              float *phi,
              float *kdp,
              int *eclass,
              int *rtype,
              short *melt,
              float alpha,
              float *rate);
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
               int nazi);
void getRayRR_RAT(float *z,
                  float *phi,
                  float *kdp,
                  float *t,
                  int *eclass,
                  int *rtype,
                  float alpha,
                  float *rate);
float getIntercept(float t);
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
                  float *rate);
void estimateRZC(struct volume_info header,
                 short ***z,
                 int ***kdp,
                 int ***echo,
                 short ***rtype,
                 float ***rr,
                 int angle);
void estimateRZCS(struct volume_info header,
                  short ***z,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle);
void estimateRZDR(struct volume_info header,
                  short ***z,
                  int ***zdr,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle,
                  int scale);
void estimateRKDP(struct volume_info header,
                  short ***z,
                  int ***kdp,
                  int ***echo,
                  short ***rtype,
                  float ***rr,
                  int angle,
                  int scale);
void estimateRZCwVPR(struct volume_info header,
                     struct range_parameter evpr,
                     short ***z,
                     int ***echo,
                     short ***rtype,
                     float **vpr,
                     float *eangle,
                     int angle,
                     float ***rr);
void estimateRZCSwVPR(struct volume_info header,
                      struct range_parameter evpr,
                      short ***z,
                      int ***echo,
                      short ***rtype,
                      float **vpr,
                      float *eangle,
                      int angle,
                      float ***rr);
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
                      float ***rr);
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
                      float ***rr);

void get2DFixedGrid(float *razimuth,
                    float BWidth,
                    float **data);
void copy3Darray(int ***in,
                 int ***out);
float int_dBZ2float_dBZ(short int_dBZ);
short float_dBZ2int_dBZ(float float_dBZ);
float getSystemPhidp(int ***phi,
                     int scale,
                     int ***flag);
void averageProduct(float **rr,
                    int dw);
float getAverage(float **rr,
                 int azi,
                 int rng,
                 int max_azi,
                 int max_rng,
                 int dw);
int getEpochTime(struct volume_info header);






/* averaging Z product */
void averageZProduct(float (*product)[MAXIMUM_RANGE_SUPER],
                     float **new,
                     int max_azi,
                     int max_rng,
                     int dxy);
float getZmean(float (*product)[MAXIMUM_RANGE_SUPER],
               int azi,
               int rng,
               int max_azi,
               int max_rng,
               int dxy);
