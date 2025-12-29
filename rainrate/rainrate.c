/******************************************************************
--> rainrate_sa.c

    : Estimate rain rate using the specific attenuation method



*******************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "libphydro.h"
#include "rsl.h"

char *infile, *nwppath, *output_path, *huc;
int box_azi1, box_azi2, box_bin1, box_bin2;
int product_type, norain, ap_flag, range_flag;
float cappi_height=0., tower_height=0., gaussian_width;
int qpe_method, range_start, range_end, max_height, width_sector, move_sector;
double dheight;

int  getopt(argc, argv, opts);
void getArgs(int argc, char **argv);

int main(int argc, char *argv[])
{
    FILE *fout;
    struct volume_info header;
    struct bounding box;
    struct range_parameter evpr;
    char huc_string[9]="";
	int i, rn, bn, digit_huc, i_hcappi, dw, *status, ini_status=FALSE;
    float **rainrate;
    char date_string[12], time_string[8], year_string[5], day_string[3],
         hour_string[7], min_string[3], sec_string[3], zero_string[2]="0",
         bar_string[2]="_", tem_string[3], irc_string[2], iap_string[2],
         out1_string[14]=".out", norain_string[3], hcappi_string[3] ;
    char fname[200], temppath[200];
    char *chmon[12]={"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                     "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

    /* get arguments */
    getArgs(argc,argv);

    /* parameter setting */
    box.azi1=box_azi1;
    box.azi2=box_azi2;
    box.bin1=box_bin1;
    box.bin2=box_bin2;
    evpr.start_range=range_start;
    evpr.end_range=range_end;
    evpr.sector_width=width_sector;
    evpr.sector_move=move_sector;
    evpr.hmax=max_height;
    evpr.dh=dheight;
    status=&ini_status;

    /* Check HUC */
    if (atof(huc)==0) {
       for (i=0; i<8; i++){
            huc_string[i]='0';
       }
    }
    else {
        /* if HUC is smaller than 8 digits, fill "0" like "00001234" */
        digit_huc=strlen(huc);
        if (digit_huc!=8){
           for (i=0; i<8; i++){
               if (i<8-digit_huc)
                   huc_string[i]='0';
               else
                   huc_string[i]=huc[i-(8-digit_huc)];
           }
        }
        else {
           for (i=0; i<8; i++){
                huc_string[i]=huc[i];
           }
        }
    }

    /* define array */
    rainrate=get2dFloatArray(MAXIMUM_AZIMUTH_SUPER,MAXIMUM_RANGE_SUPER);

    /* Product selection */
    if (product_type>=1 && product_type<=10) {
        getRain(infile,nwppath,product_type,cappi_height,tower_height,ap_flag,range_flag,
                norain,&header,box,evpr,rainrate,qpe_method,status);
    }
    else {


    }

    /* average product */
    dw=1; // averaging window size
    averageProduct(rainrate,dw);

    if (*status == TRUE)
        fprintf(stderr, "Preprocessing: %s\n", infile);
    else {
        fprintf(stderr, "Program will be terminated at %s\n",infile);
        exit(-1);
        return -1;
    }

    /* Output file generation ---------------------------------------------- */
    strcpy(fname,header.radar_id);
    strcat(fname,huc_string);
    strcat(fname,bar_string);

    /* date & time */
    sprintf(year_string, "%d", header.year); // year
    sprintf(day_string, "%d", header.day); // day
    if (header.day<10){
        strcpy(tem_string,day_string);
        strcpy(day_string,zero_string);
        strcat(day_string,tem_string);
    }
    strcpy(date_string, day_string);
    strcat(date_string, chmon[header.mon-1]);
    strcat(date_string, year_string);

    sprintf(hour_string, "%d", header.hour); // hour
    if (header.hour<10){
        strcpy(tem_string,hour_string);
        strcpy(hour_string,zero_string);
        strcat(hour_string,tem_string);
    }

    sprintf(min_string, "%d", header.min); // minute
    if (header.min<10){
        strcpy(tem_string,min_string);
        strcpy(min_string,zero_string);
        strcat(min_string,tem_string);
    }

    sprintf(sec_string, "%d", header.sec); // second
    if (header.sec<10){
        strcpy(tem_string,sec_string);
        strcpy(sec_string,zero_string);
        strcat(sec_string,tem_string);
    }

    strcpy(time_string,bar_string);
    strcat(time_string,hour_string);
    strcat(time_string,min_string);
    strcat(time_string,sec_string);

    /* product type */
    if (product_type==0) // Volume scan
        strcat(fname, "V");
    else if (product_type==1) // Base scan
        strcat(fname, "B");
    else if (product_type>=2 && product_type<=9) // Elevation angle
        strcat(fname, "A");
    else if (product_type==10) // CAPPI
        strcat(fname, "C");
    else if (product_type==11) // PPS hybrid scan
        strcat(fname, "P");

    /* bin size (default: 1km) */
    strcat(fname, "1");

    /* threshold for no rain */
    if (norain<10)
        strcat(fname, "0");
    sprintf(norain_string, "%d", norain);
    strcat(fname, norain_string);

    /* CAPPI */
    if (product_type==10) {
        i_hcappi=(int)(cappi_height*10);
        sprintf(hcappi_string, "%d", i_hcappi);
        strcat(fname, hcappi_string);
    }
    else {
        strcat(fname,"0");
        sprintf(iap_string, "%d", ap_flag);
        strcat(fname, iap_string);
    }

    /* rain rate conversion: spectific attenuation */
    strcat(fname, "S");

    /* Range correction */
    sprintf(irc_string, "%d", range_flag);
    strcat(fname, irc_string);

    /* Advection correction: None */
    strcat(fname, "0");

    /* Coordinate system: Polar */
    strcat(fname, "P");
    strcat(fname, bar_string);
    strcat(fname, "00"); // accumulation interval
    strcat(fname, bar_string);
    strcat(fname, date_string);
    strcat(fname, time_string);
    strcat(fname, out1_string);
    strcpy(temppath,output_path);
    strcat(temppath,fname);
    fout=fopen(temppath,"w");

    if (product_type==11) {

    }
    else {
        for (rn=0; rn<MAXIMUM_AZIMUTH_SUPER; rn++) {
             for (bn=0; bn<MAXIMUM_RANGE_SUPER; bn++) {
                  if (bn/4<box.bin1 || bn/4>box.bin2)
                      continue;
                  else {
                     if (box.azi1<box.azi2) {
                         if (rn/2<box.azi1 || rn/2>box.azi2)
                             continue;
                         else
                             fprintf(fout, "%d ",(int)(rainrate[rn][bn]*RATE_SCALING+0.5));
                     }
                     else if (box.azi1>box.azi2){
                         if (rn/2>box.azi2 && rn/2<box.azi1)
                             continue;
                         else
                             fprintf(fout, "%d ",(int)(rainrate[rn][bn]*RATE_SCALING+0.5));
                     }
                     if ((int)(bn*0.25+0.25)==box.bin2+1)
                         fprintf(fout,"\n");
                  }
             }
        }
    }
    fclose(fout);
    fprintf(stdout,temppath);
    return 0;
}
void showUsage()
{
  fprintf(stderr,"Usage: ascindex [options] file\n\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr," -i str            Input file\n");
  fprintf(stderr," -p str            Output path\n");
  fprintf(stderr," -d str            NWP path\n");
  fprintf(stderr," -u str            HUC (Hydrologic Unit Code)\n");
  fprintf(stderr," -x int1, 2, 3, 4  Bounding Box.\n");
  fprintf(stderr," -a int            QPE method (estimator)\n");
  fprintf(stderr,"                     1. R(Z): NEXRAD Z-R\n");
  fprintf(stderr,"                     2. R(Z): Based on classification (NEXRAD and M-P)\n");
  fprintf(stderr,"                     3. R(Z,Zdr)\n");
  fprintf(stderr,"                     4. R(Kdp)\n");
  fprintf(stderr,"                     5. R(A)\n");
  fprintf(stderr,"                     6. R(A)with optimized alpha\n");
  fprintf(stderr," -y int            Product type\n");
  fprintf(stderr," -n int            Threshold for no rain (dBZ)\n");
  fprintf(stderr," -q                AP algorithm \n");
  fprintf(stderr," -r                Range correction\n");
  fprintf(stderr," -h flt            CAPPI height (km)\n");
  fprintf(stderr," -t flt            Radar tower height (m)\n");
  fprintf(stderr," -w flt            Smoothing width (km)\n");
  fprintf(stderr," -v flt            Vertical resolution for VPR (km)\n");
  fprintf(stderr," -g int            Maximum height for VPR (km)\n");
  fprintf(stderr," -e int 1, 2       Analysis range for VPR (from, to)\n");
  fprintf(stderr," -k int            Sector width for VPR (degree)\n");
  fprintf(stderr," -m int            Azimuthal movement for smoothing (degree)\n");
  exit(-1);
}
void getArgs(int argc, char * argv[])
{
  int   ch, errflg, n1, n2, n3, n4, n_sum;
  int n_i, n_p, n_u, n_x, n_y, n_n, n_q, n_r, n_h, n_t,
      n_w, n_v, n_g, n_e, n_k, n_m;
  int huc_flag, i_type, i_rc;
  float tmpf;
  extern char *optarg;
  extern int optind,optopt;

  errflg = 0;
  if (argc <2) errflg++;

  /* for duplicate input argument */
  n_i=0;  n_p=0;  n_u=0; n_x=0;
  n_y=0;  n_n=0;  n_q=1; n_r=1;
  n_h=0;  n_t=0;  n_w=0; n_v=0;
  n_g=0;  n_e=0;  n_k=0; n_m=0;

  huc="0";
  ap_flag = FALSE;
  range_flag = FALSE;
  huc_flag = FALSE;
  box_azi1 = 0;
  box_azi2 = 359;
  box_bin1 = 0;
  box_bin2 = 229;
  i_type = FALSE;
  i_rc = FALSE;

  /* Default behaviour */
  while((ch=getopt(argc,argv,"i:p:d:u:x:a:y:n:h:t:w:v:g:e:k:m:qr"))!=-1){
	  switch (ch) {

	  case 'i':
		  infile=optarg;
		  n_i=n_i+1;
		  break;

	  case 'p':
		  output_path=optarg;
		  n_p=n_p+1;
		  break;

	  case 'd':
		  nwppath=optarg;
		  break;

	  case 'u':                  /* HUC (Hydrologic Unit Code) */
		  huc = optarg;
		  huc_flag = TRUE;
		  n_u=n_u+1;
		  break;

	  case 'x':                   /* Bounding Box*/
		  sscanf(optarg,"%d%d%d%d",&n1,&n2,&n3,&n4);
		  box_azi1 = n1;
		  box_azi2 = n2;
		  box_bin1 = n3;
		  box_bin2 = n4;
		  if (box_bin2 > 229)
			  box_bin2 = 229;
		  n_x=n_x+1;
		  break;

	  case 'a':                  /* QPE method (rainfall estimator) */
	      tmpf=(float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -y requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  qpe_method=(int)tmpf;

		  break;


	  case 'y':                  /* Product type */
		  tmpf = (float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -y requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  product_type = (int)tmpf;
		  if (product_type == 10)
			  i_type = TRUE;
		  n_y=n_y+1;
		  break;

	  case 'n':                  /* Threshold for no rain (dBZ) */
		  tmpf = (float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -n requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  norain = (int)tmpf;
		  n_n=n_n+1;
		  break;

	  case 'q':                  /* AP algorithm */
		  ap_flag = TRUE;
		  break;

	  case 'r':                  /* Range correction */
		  range_flag = TRUE;
		  i_rc = TRUE;
		  break;

	  case 'h':                   /* CAPPI height */
		  tmpf = (float)atof(optarg);
		  if (!isdigit((int)*optarg) || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -h requires float argument > 0\n");
			  errflg++;
			  break;
		  }
		  cappi_height = tmpf;
		  n_h=n_h+1;
		  break;

	  case 't':                   /* Radar tower height */
		  tmpf = (float)atof(optarg);
		  if (!isdigit((int)*optarg) || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -t requires float argument > 0\n");
			  errflg++;
			  break;
		  }
		  tower_height = tmpf;
		  n_t=n_t+1;
		  break;

	  case 'w':                   /* Smoothing width (km) */
		  tmpf = (float)atof(optarg);
		  if (!isdigit((int)*optarg) || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -w requires float argument > 0\n");
			  errflg++;
			  break;
		  }
		  gaussian_width = tmpf;
		  n_w=n_w+1;
		  break;

	  case 'v':                   /* Vertical resolution for VPR (km) */
		  tmpf = (float)atof(optarg);
		  if (!isdigit((int)*optarg) || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -v requires float argument > 0\n");
			  errflg++;
			  break;
		  }
		  dheight = tmpf;
		  n_v=n_v+1;
		  break;

	  case 'g':                  /* Maximum height for VPR (km) */
		  tmpf = (float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -g requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  max_height = (int)tmpf;
		  n_g=n_g+1;
		  break;

	  case 'e':                   /* Analysis range for VPR */
		  sscanf(optarg,"%d%d",&n1,&n2);
		  range_start = n1;
		  range_end = n2;
		  n_e=n_e+1;
		  break;

	  case 'k':                  /* Sector width for VPR (degree) */
		  tmpf = (float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -k requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  width_sector = (int)tmpf;
		  n_k=n_k+1;
		  break;

	  case 'm':                  /* Azimuthal movement for smoothing (degree) */
		  tmpf = (float)atof(optarg);
		  if ((tmpf-(int)tmpf)>0 || tmpf < 0.0) {
			  fprintf(stderr,"ERROR: -m requires interger argument > 0\n");
			  errflg++;
			  break;
		  }
		  move_sector = (int)tmpf;
		  n_m=n_m+1;
		  break;

	  case ':':
		  fprintf(stderr,"Option -%c requires and argument\n",optopt);
		  errflg++;
		  break;

	  case '?':
		  fprintf(stderr,"Unrecognized option: -%c\n",optopt);
		  errflg++;
	  }
  }

  /* Add additional checks here */

  if (errflg > 0)
	  showUsage();

  if (huc_flag == TRUE) {
	  if (i_type == TRUE && range_flag == TRUE) {
		  n_sum = n_i+n_p+n_u+n_x+n_y+n_n+n_q+n_r+n_h+n_t+n_w+n_v+n_g+n_e+n_k+n_m;
		  if (n_sum < 16) {
			  if ((16- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 16-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 16) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else if (i_type == TRUE && range_flag == FALSE) {
		  n_sum = n_i+n_p+n_u+n_x+n_y+n_n+n_q+n_r+n_h+n_t+n_w;
		  if (n_sum < 11) {
			  if ((11- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 11-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 11) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else if (i_type == FALSE && range_flag == TRUE) {
		  n_sum = n_i+n_p+n_u+n_x+n_y+n_n+n_q+n_r+n_v+n_g+n_e+n_k+n_m;
		  if (n_sum < 13) {
			  if ((13- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 13-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 13) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else {
		  n_sum = n_i+n_p+n_u+n_x+n_y+n_n+n_q+n_r;
		  if (n_sum < 8) {
			  if ((8- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 8-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 8) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
  }

  else {
	  if (i_type == TRUE && range_flag == TRUE) {
		  n_sum = n_i+n_p+n_y+n_n+n_q+n_r+n_h+n_t+n_w+n_v+n_g+n_e+n_k+n_m;
		  if (n_sum < 14) {
			  if ((14- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 14-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 14) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else if (i_type == TRUE && range_flag == FALSE) {
		  n_sum = n_i+n_p+n_y+n_n+n_q+n_r+n_h+n_t+n_w;
		  if (n_sum < 9) {
			  if ((9- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 9-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 9) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else if (i_type == FALSE && range_flag == TRUE) {
		  n_sum = n_i+n_p+n_y+n_n+n_q+n_r+n_v+n_g+n_e+n_k+n_m;
		  if (n_sum < 11) {
			  if ((11- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 11-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 11) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
	  else {
		  n_sum = n_i+n_p+n_y+n_n+n_q+n_r;
		  if (n_sum < 6) {
			  if ((6- n_sum) ==1) {
				  fprintf(stderr, "ERROR: Needs one more argument\n");
				  exit(-1);
			  }
			  else {
				  fprintf(stderr, "ERROR: Needs %d more arguments\n", 6-n_sum);
				  exit(-1);
			  }
		  }
		  if (n_sum > 6) {
			  fprintf(stderr, "ERROR: Duplicate arguments\n");
			  exit(-1);
		  }
	  }
  }
}
