# Radar data processing libraries (/libphydro_v2/)

NEXRAD Level-II radar data procesing library for polarimetric rainfall estimation

- compile: run compile.sh and create a library file (libphydro.a)


# Rain rate estimation using polarimetrci algorithms (/rainrate/)

- dependency: NASA RSL (https://github.com/adokter/rsl)
- compile: run compile.sh and create an executable file (e.g., rainrate);
- example data (/rainrate/SampleData/): Level II (KEAX, Kansas City, MO) and NWP temperature and look-up data
- arguments:    
   -i str &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp; Input file  
   -p str &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Output path  
   -d str &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; NWP path  
   -u str &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; HUC (Hydrologic Unit Code)  
   -x int1, 2, 3, 4 &emsp;&emsp;&ensp; Bounding Box  
   -a int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; QPE method (estimator)  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 1. R(Z): NEXRAD Z-R  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 2. R(Z): Based on classification (NEXRAD and M-P)  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 3. R(Z,Zdr)  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 4. R(Kdp)  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 5. R(A)  
   &emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 6. R(A)with optimized alpha  
   -y int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Product type (PPI: 1-9, CAPPI: 10)  
   -n int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Threshold for no rain (dBZ)  
   -q &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; AP algorithm  
   -r &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; Range correction  
   -h flt &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; CAPPI height (km)  
   -t flt &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; Radar tower height (m)  
   -w flt &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Smoothing width (km)  
   -v flt &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Vertical resolution for VPR (km)  
   -g int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp; Maximum height for VPR (km)  
   -e int 1, 2 &emsp;&emsp;&emsp;&emsp;&ensp; Analysis range for VPR (from, to)  
   -k int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp; Sector width for VPR (degree)  
   -m int &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Azimuthal movement for smoothing (degree)  
