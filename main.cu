// TWiCE; Tephra fall simulator for Windy Condition Eruption
// Developed by K. Mannen in 2022
// Modified in 20230714 TO ASSIGN kw and ks from configfile
// Modified on 20230909 TO CALCULATE PLUME TRAJECTORY USING INDEPENDENT ds
// Modified on 20240728 OUTPUT FALL DATA FOR EACH 0.1 PHI INTERVAL; CORRECTION OF CALCULATION; PARAMETER NAME CHANGE, INTERVAL_PHI -> INTERVAL_DECIMAL_PHI;
// 						Output decimal_falldriftX_#.txt, decimal_falldriftY_#.txt, decimal_falltime_#.txt
// 						Output time after vent in plume.txt file
//						U = V after reaching Hb
// Modified on 20240729 OUTPUT particle_segregation_##.txt, which depict amount of segregation from each position along plume axis s and phidec sizeclass
// 20241119 Typo correction 255 -> 225 in L1518
// 20260426 Starting Add Comments
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

/* CUDA Check*/
#include <cuda_runtime.h>

#define CUDA_CHECK(call) do {                                      \
    cudaError_t err = (call);                                      \
    if (err != cudaSuccess) {                                      \
        fprintf(stderr,                                            \
                "CUDA error at %s:%d: %s\n",                       \
                __FILE__, __LINE__, cudaGetErrorString(err));      \
        exit(EXIT_FAILURE);                                        \
    }                                                              \
} while (0)

#define CUDA_KERNEL_CHECK() do {                                   \
    CUDA_CHECK(cudaGetLastError());                                \
    CUDA_CHECK(cudaDeviceSynchronize());                           \
} while (0)
// end of CUDA Check

//#define TEPHRA2
#define CUDA
//#define TEST	//OUTPUT lspml.txt which is massloading contribution for each Location, particle Source, Phi in decimal

double  GRAVITY = 9.81;

/* The following Global Variables are assigned some default values */

double DIFFUSION_COEFFICIENT = 200.0;
double FALL_TIME_THRESHOLD = 3600.0;
double EDDY_CONST = .04;
double PUMICE_DENSITY = 1000.0;

double WIND_INTERVAL;
double PLUME_HEIGHT = -9999;
double ERUPTION_MASS = 1e10;

double MAX_GRAINSIZE = -7.0;
double MIN_GRAINSIZE = 7.0;
double INTERVAL_DECIMAL_PHI = 0.1;

int COLLAPSE_THEN_OFF = 1;	// If eruption plume collapse, tephra dispersal won't be calculated. 20250105
int PHIDECDIM;
int WRITE_DEPCENT_TRAJECTORY = 0;
int WRITE_FALL_INFO_FILES = 0;
int WRITE_COLUMN_FILES = 0;
int WRITE_CONF = 0;
int WRITE_SDIMCUTOFF = 0;
int WRITE_MASSLOADING = 1;
int WRITE_DECIMAL_MASSLOADING = 0;
int WRITE_DECIMAL_FALL_TRAJ = 0;

double MEDIAN_GRAINSIZE = -1;
double STD_GRAINSIZE = 2.0;

double VENT_EASTING = 0.0;
double VENT_NORTHING = 0.0;
double VENT_ELEVATION = 0.0;

double INITIAL_WATER_CONTENT = -9999;
double MAGMA_DISCHARGE_RATE = -9999;
double MAGMA_TEMPERATURE = 1200;
double MESH_SIZE_IN_KM = 1.0;
double INITIAL_PLUME_VELOCITY = -9999;
double VENT_RADIUS = -9999;
double MINIMUM_DEPOSIT_FOR_MD_CALC = 1;
double PLUME_THICKNESS = 800;
double PLUME_RADIUS_CORRECTION = 1.0;

double MINIMUM_CONTRIBUTION = -9999;	// Minimum contribution of particle segregation to mass loading (massloading on land per particle segregation kg/sq-m / m)

double Ht;

double S_MAX= 100000;
int SDIM_FOR_PLUME_CALC = -9999;
int SDIM_FOR_FALL_CALC = -9999;
int SDIMCUTOFF;
double S_DELTA_FOR_PLUME_CALC = 10;
double S_DELTA_FOR_FALL_CALC = 100;
double Z_DELTA = 100;
int ZDIM;   // number of h interval

int LOCDIM; // number of locations to calc


// ORIGINALLY IN WINDY.C
double ENTRAIN_COEFF_KS = 0.09;  // another k should be introduced for gas thurst region but uniform value in this code
double ENTRAIN_COEFF_KW = 0.9;

// MAP BOUNDARY
double MAPENDE = -9999999;
double MAPENDW =  9999999;
double MAPENDS =  9999999;
double MAPENDN = -9999999;

typedef struct {
int j;
double x;
double y;
double z;
double dist;
double ttlmassloading;
double smallerthan1mm;
double meandiameter;
double dep[20]; // mass of each phi size class but for l2 dep[0] and dep[1] indicate area of isopach and its square root
} DEP;	//locdatastruct or l[j]

typedef struct {	// Total particle segregation during the eruption
double phi;
double theoretical;	// Deduced from particle distribution function
double actual;		// Can be lower than theoretical because of S_MAX shortage especially for small particles
} RELEASE;

typedef struct {	// Amount of particle segrigation from an interval on the plume axis
	double mass_from_ds[20]; //[] means number of size classes in phi scale
} SEG;

typedef struct {
  int day;
  int hour;
  double wind_height; /* height a.s.l. in km */
  double wind_speed; 	/* the average windspeed in m/s */
  double wind_dir;  	/* average wind direction in +/- degrees from north */
  double t_atm;
  double p_atm;
} WIND;
static WIND *W1;

#define MAX_LINE 200
#define M_2PI 2*M_PI

//prototypes
int init_globals(char *config_file);
void printparticlereleased(RELEASE *r);
void printsegregation_per_ds(SEG *r);
void printxyz(FILE *in, const char *type, int i, double *srcX, double *srcY, double *srcZ);
void printxyzq(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q, double *t);
void printxyze(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q);
void printfallsummary(const char *c, double *h, double *ttlfalltimesummary, double base, int stepnumber, double stepdelta); // modified on 2024.07.28
void printsegsummary(const char *c, double *massreleased_per_ds_and_phidec, double base, int stepnumber, double stepdelta);


void printttlml(int phiint, double *ttlml);
void printlspml(int phiint, double *lspml);
void printS4eachphidec(DEP *l, int phiint, double *lspml);

int get_line_number(FILE *in_wind);
int get_wind_line_number(FILE *in_wind);
void get_sdimcutoff(double *sigma_square, SEG *massreleased_per_ds_and_phidec, int phiint);

void atmosphere(int windlinenum, double *h, double *atmT, double *atmP, double *windX, double *windY, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);
void interval_fall_calc(int zmax, int phidecimal, double grainsize, double *h, double *atmP, double *atmT, double *windX, double *windY, double *driftX, double *driftY, double *ttlfalltime);
void writettlfallsummary(int phiint, double *ttlfalltime, double *ttlfalltimesummary);
void drift_from_a_certain_source(double *source_x, double *source_y, double *source_height, double *sourceRadius, double *TotalFallTime, double *driftX, double *driftY, double *driftX_s, double *driftY_s, double *sigma_square);
void writetrajectory(int phiint, double *driftX_s, double *driftY_s, double *sigma_squre, double *massreleased, SEG *seg);
double calc_sigma_square(double source_radius, double falltime);
void mass_release_calc(int zmax, int phidecimal, double phi, double *h, double *atmP, double *atmT, double *windX, double *windY, double *massreleased);
double confirm_released_mass(double *massreleased);
void locwrite(DEP *l, double *locX, double *locY, double *locZ);
void depwrite(int size, DEP *l, double *massloading);
void sumwrite(DEP *l, double *ttlmassloading, double *cummassphi);
void printdeposit(DEP *locdatastruct);
void clearary(int dim, double *ary);
int compare_ttlmassloading(const void *a, const void *b);
int compare_Md(const void * a, const void * b);
void create_source_array(double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *sourceT, double *plume_trajX, double *plume_trajY, double *plume_traj_Z, double *plume_trajR, double *plume_trajT);

void createisopachdata(DEP *l);
void extractisopachdata(DEP *l);
void countmeandiameter(DEP *l);
double calc_dir(double x, double y);

void obtaintheoreticallyreleased(RELEASE *r);

//External in windy.c
void read_wind(FILE *f, double *h, double *v, double *d, double *t, double *p);
void read_loc(FILE *f, double *x, double *y, double *z);

double windy(int linenum, double *sourceX, double *sourceY, double *sourceZ, double *sourceRadius, double *timeaftervent, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);
void rk(int, double);
void makewindstruct(int imax, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);

double func12(double, double, double, double);
double func13(double, double, double, double, double);
double func14(double, double, double, double, double);
double func15(double, double, double, double, double, double, double);
double func16(double, double, double, double);

double func17(double, double, double, double);				// calc plume density
double func18(double);				// particle content
double func19(double);				// Rg calc
double func20(double);				// Cp calc

double calcCp0(void);			// Cp0 calc (eq. 20.5; in text between eq20 and 21)
double calc_Tatm(double, int);
double calc_Patm(double, int);
double func22(double, double);	// pressure profile
double func23(double, double);	// atmospheric density
double get_V(double, int);
double get_dir(double, int);

/* Non-Cuda Functions (start)*/
void calc_mass_loading_element(int phisize, double *sourceZ, double *driftX_s, double *driftY_s, double *sigma_square, double *locX, double *locY, double *locZ, double *mlj, double *massreleased);
void calc_mass_loading_location(int phiint, double *mlj, double *massloading, double *ttl, double *cummassphi);
/* Non-Cuda Functions (end)*/

// CUDA function
#ifdef CUDA
void calc_mass_loading(double *sourceZ, double *driftX_s, double *driftY_s, double *sigma_square, double *locX, double *locY, double *locZ, double *lspml, double *ttlml, double *massreleased);
void cumulative_mass_x_phi(int phiint, double *lspml, double *massloading, double *ttl, double *cummassphi);
__global__ void funcD01a(int, int, int, int, float, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *);
__global__ void funcD01b(int N, int LOCDIM, float *lspmlD, float *ttlmlD);
#endif

//External in grain.c
double calc_particle_terminal_velocity(double h, double ashdiam, double part_density, double p, double t);
void phigenerator();
void phiconvert();
double calc_pdf_fraction(double phi);

int main(int argc, char *argv[]) {
	/* read config file */
	char string[30];
	init_globals(argv[1]);		// to set WRITE_CONF which shows parameter in conf before run
	init_globals(argv[1]);
	PHIDECDIM = (int)(1 / INTERVAL_DECIMAL_PHI);
	phiconvert();   // make max and minimum grain sizes ordered

	/* read wind file */
	FILE *in_wind;
	int windlinenum = 0;
	in_wind = fopen(argv[2], "r");
	windlinenum = get_wind_line_number(in_wind);
	//printf("\n\nthe total line number of wind file is %d\n", windlinenum);
	rewind(in_wind);

	double *wind_alt, *wind_v, *wind_dir, *wind_tmp, *wind_pres;
	wind_alt = (double *)malloc(windlinenum * sizeof(double));
	wind_v = (double *)malloc(windlinenum * sizeof(double));
	wind_dir = (double *)malloc(windlinenum * sizeof(double));
	wind_tmp = (double *)malloc(windlinenum * sizeof(double));
	wind_pres = (double *)malloc(windlinenum * sizeof(double));

	read_wind(in_wind, wind_alt, wind_v, wind_dir, wind_tmp, wind_pres);
	fclose(in_wind);

	/* read loc file */
	FILE *in_loc;
	in_loc = fopen(argv[3], "r");
	LOCDIM = get_line_number(in_loc);
	//printf("\n\nthe total line number of location file is %d\n", LOCDIM);
	rewind(in_loc);

	double *locX, *locY, *locZ;
	locX = (double *)malloc(LOCDIM * sizeof(double));
	locY = (double *)malloc(LOCDIM * sizeof(double));
	locZ = (double *)malloc(LOCDIM * sizeof(double));

    read_loc(in_loc, locX, locY, locZ);
	fclose(in_loc);
	/* end of read loc file */

	SDIM_FOR_PLUME_CALC = (S_MAX + S_DELTA_FOR_PLUME_CALC - 1) / S_DELTA_FOR_PLUME_CALC; // dimension of plume calculation; calculated plume should be longer than source plume
	SDIM_FOR_FALL_CALC = S_MAX / S_DELTA_FOR_FALL_CALC; // dimension of source
	
	//printf("SDIM_P=%d\tSDIM_F=%d\n", SDIM_FOR_PLUME_CALC, SDIM_FOR_FALL_CALC);
	
	SDIMCUTOFF = SDIM_FOR_FALL_CALC;

	// Trajectory of plume center for particle calculation
	double *plume_trajX, *plume_trajY, *plume_trajZ, *plume_trajR, *plume_trajT;
	plume_trajX = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
	plume_trajY = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
	plume_trajZ = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
	plume_trajR = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
	plume_trajT = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double)); //20240828 time from vent to a certain s
	Ht = windy(windlinenum, plume_trajX, plume_trajY, plume_trajZ, plume_trajR, plume_trajT, wind_alt, wind_v, wind_dir, wind_tmp, wind_pres);
	
	// Trajectory of plume center for FALL CALCULATION (array of point data of particle source along the plume)
	double *sourceX, *sourceY, *sourceZ, *sourceRadius, *sourceT;
	sourceX = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
	sourceY = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
	sourceZ = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
	sourceRadius = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
	sourceT = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
	
	create_source_array(sourceX, sourceY, sourceZ, sourceRadius, sourceT, plume_trajX, plume_trajY, plume_trajZ, plume_trajR, plume_trajT);
	
	// Write original trajectory file
	if(WRITE_COLUMN_FILES){
		FILE *outfile;
		outfile = fopen("plumetraj.txt", "w");
    const char *header2 = "calc_step\tx\ty\tz\tR\ttime\n";
		printxyzq(outfile, header2, SDIM_FOR_PLUME_CALC, plume_trajX, plume_trajY, plume_trajZ, plume_trajR, plume_trajT);
		fclose(outfile);
	}
	
	// Write original trajectory file
	if(WRITE_COLUMN_FILES){
		FILE *outfile;
		outfile = fopen("plumesourceposition.txt", "w");
    const char *header2 = "source\tx\ty\tz\tR\ttime\n";
		printxyzq(outfile, header2, SDIM_FOR_FALL_CALC, sourceX, sourceY, sourceZ, sourceRadius, sourceT);
		fclose(outfile);
	}

	// Obtain atmospheric condition at each interval
	double *h, *atmT, *atmP, *windX, *windY;
	int zmax = ceil(Ht / Z_DELTA);
	ZDIM = zmax + 1;
	h = (double *)malloc(ZDIM * sizeof(double));
	atmT = (double *)malloc(ZDIM * sizeof(double)); atmP = (double *)malloc(ZDIM * sizeof(double));
	windX = (double *)malloc(ZDIM * sizeof(double)); windY = (double *)malloc(ZDIM * sizeof(double));
	for(int z = 0; z < zmax; z++){
		h[z] = z * Z_DELTA;
	}
	h[zmax] = Ht;
	atmosphere(windlinenum, h, atmT, atmP, windX, windY, wind_v, wind_dir, wind_tmp, wind_pres);

	// Obtain fall time at each 
	/*The Loop 158; main loop for mass loading calculation for each location on the ground*/
	// Grain size loop consists of two loops
	//
	int phisize;  // phi number of the fraction 
	double grainsize, phi;
	//double released_mass_of_fraction = 0.0;

	double *driftX, *driftY;
	double *ttlfalltime;
	ttlfalltime = (double*)calloc(ZDIM, sizeof(double));
	driftX = (double*)calloc(ZDIM, sizeof(double));
	driftY = (double*)calloc(ZDIM, sizeof(double));

	double *ttlfalltimesummary, *ttldriftXsummary, *ttldriftYsummary;
	ttlfalltimesummary = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
	ttldriftXsummary = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
	ttldriftYsummary = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
	
	double *ttlfalltimephidec, *ttldriftXphidec, *ttldriftYphidec;
	ttlfalltimephidec = (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));
	ttldriftXphidec =   (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));
	ttldriftYphidec =   (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));

	double *massreleased_per_ds_and_phidec;
	massreleased_per_ds_and_phidec = (double*)calloc(SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
	
	SEG *massreleased_per_ds;
	massreleased_per_ds = (SEG*)calloc(SDIM_FOR_FALL_CALC, sizeof(SEG));

	double *driftX_s, *driftY_s, *sigma_square;
	driftX_s = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
	driftY_s = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
	sigma_square = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));

	double *tmpmassloading, *ttlmassloading, *cummassphi;
	tmpmassloading = (double*)calloc(LOCDIM, sizeof(double));
	ttlmassloading = (double*)calloc(LOCDIM, sizeof(double));
	cummassphi = (double*)calloc(LOCDIM, sizeof(double));

	DEP *locdatastruct;
	locdatastruct = (DEP *)calloc(LOCDIM, sizeof(DEP));
	locwrite(locdatastruct, locX, locY, locZ);

	RELEASE *r;
	r = (RELEASE *)calloc(MIN_GRAINSIZE - MAX_GRAINSIZE, sizeof(RELEASE));
	obtaintheoreticallyreleased(r);

	clearary(LOCDIM, ttlmassloading);
	clearary(LOCDIM, cummassphi);

	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		for(int phidecimal = 0; phidecimal < PHIDECDIM; phidecimal++){

			phi = phiint + MAX_GRAINSIZE + 1 - phidecimal * INTERVAL_DECIMAL_PHI;
			//printf("INTERVAL_DECIMAL_PHI=%1.1f\tPHIDECDIM = %d\tphi = %1.1f\n", INTERVAL_DECIMAL_PHI, PHIDECDIM, phi);
			grainsize = pow(2, -phi) * 0.001;	// grain size in mm

			/* F10 Calculate drift center and falltime for each height interval */
			interval_fall_calc(zmax, phidecimal, grainsize, h, atmP, atmT, windX, windY, driftX, driftY, ttlfalltime);
			
			//	20240728 store 0.1 phi interval fallout
			writettlfallsummary(phidecimal, ttlfalltime, ttlfalltimephidec); // 20240728
			writettlfallsummary(phidecimal, driftX, ttldriftXphidec); // 20240728
			writettlfallsummary(phidecimal, driftY, ttldriftYphidec); // 20240728
			
			/*printf("grainsize = %1.4e m\n", grainsize);
			for(int t = 0; t < ZDIM; t++){
				printf("%1.4f,%1.4f,%1.4f\n", ttlfalltime[t], driftX[t], driftY[t]);
			}*/
			// end of 20240728
			
			// output 1 phi interval fallout
			if(phidecimal == 0){
				writettlfallsummary(phiint, ttlfalltime, ttlfalltimesummary);
				writettlfallsummary(phiint, driftX, ttldriftXsummary);
				writettlfallsummary(phiint, driftY, ttldriftYsummary);
			}
			/* Calculate particle segregation from each plume interval */
			mass_release_calc(zmax, phidecimal, phi, h, atmP, atmT, windX, windY, massreleased_per_ds_and_phidec);
		}// END OF DECIMAL PHI LOOP
		//released_mass_of_fraction = confirm_released_mass(massreleased_per_ds_and_phidec);
		drift_from_a_certain_source(sourceX, sourceY, sourceZ, sourceRadius, ttlfalltimephidec, ttldriftXphidec, ttldriftYphidec, driftX_s, driftY_s, sigma_square);
		writetrajectory(phiint, driftX_s, driftY_s, sigma_square, massreleased_per_ds_and_phidec, massreleased_per_ds);  // mass released for 1phi interval is also calculated from 0.1 phi interval data
		get_sdimcutoff(sigma_square, massreleased_per_ds, phiint);	// obtain SDIMCUTOFF
		
		if(WRITE_DECIMAL_FALL_TRAJ){// 20240728 output 0.1 phi interval fallout to each 1 phi interval file
			sprintf(string, "decimal_falltime_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			printfallsummary(string, h, ttlfalltimephidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
			sprintf(string, "decimal_falldriftX_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			printfallsummary(string, h, ttldriftXphidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
			sprintf(string, "decimal_falldriftY_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			printfallsummary(string, h, ttldriftYphidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);	
		}  // 20240728 END
		
		if(WRITE_FALL_INFO_FILES){
			sprintf(string, "particle_segregation_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			printsegsummary(string, massreleased_per_ds_and_phidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
		}

		r[phiint].phi = phiint + MAX_GRAINSIZE + 1;
		r[phiint].actual = confirm_released_mass(massreleased_per_ds_and_phidec);
		
		phisize = phiint + MAX_GRAINSIZE + 1;
		if(WRITE_SDIMCUTOFF) printf("PHI = %d\tSDIMCUTOFF = %d\tMASS_RELEASED = %1.4e\n",  phisize, SDIMCUTOFF, r[phiint].actual);
		
		//if(SDIMCUTOFF < 0){SDIMCUTOFF = 100;}
		
		if(SDIMCUTOFF > 0){
			double *lspml; // lspml stands for "local-s-phi mass loading"
						   // =  massloading data for each combination of location, s (position in plume) and phi (grain size, decimal phi)
			lspml = (double*)calloc(PHIDECDIM * SDIMCUTOFF * LOCDIM, sizeof(double));
		
#ifdef CUDA
			calc_mass_loading(sourceZ, driftX_s, driftY_s, sigma_square, locX, locY, locZ, lspml, tmpmassloading, massreleased_per_ds_and_phidec);
			cumulative_mass_x_phi(phiint, lspml, tmpmassloading, ttlmassloading, cummassphi);
#else
			calc_mass_loading_element(phisize, sourceZ, driftX_s, driftY_s, sigma_square, locX, locY, locZ, lspml, massreleased_per_ds_and_phidec);
			calc_mass_loading_location(phiint, lspml, tmpmassloading, ttlmassloading, cummassphi);
#endif
			
#ifdef TEST
			printlspml(phiint, lspml);
#endif
			if(WRITE_DECIMAL_MASSLOADING){printS4eachphidec(locdatastruct, phiint, lspml);}
			free(lspml);
		}
		
		depwrite(phiint, locdatastruct, tmpmassloading);
		clearary(LOCDIM, tmpmassloading);
		
	}// END OF INTEGER PHI LOOP

	if(WRITE_FALL_INFO_FILES){
	const char *name1 = "falltime.txt";
	printfallsummary(name1, h, ttlfalltimesummary, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	const char *name2 = "falldriftX.txt";
	printfallsummary(name2, h, ttldriftXsummary, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	const char *name3 = "falldriftY.txt";
	printfallsummary(name3, h, ttldriftYsummary, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	printparticlereleased(r);			//particle_released.txt
	printsegregation_per_ds(massreleased_per_ds); //segregation_per_ds.txt
	}

	/* write massloading.txt */
	sumwrite(locdatastruct, ttlmassloading, cummassphi);
	if(WRITE_MASSLOADING){printdeposit(locdatastruct);		//massloading.txt
	createisopachdata(locdatastruct);}						//S_vs_Area.txt
	
	free(wind_alt); free(wind_v); free(wind_dir); free(wind_tmp); free(wind_pres); 
	free(locX); free(locY); free(locZ);
	free(sourceX); free(sourceY); free(sourceZ); free(sourceRadius);
	free(plume_trajX); free(plume_trajY); free(plume_trajZ); free(plume_trajR);	// 20230808
	free(h); free(atmT); free(atmP); free(windX); free(windY);
	free(W1); //20180218
}	// End of main

void create_source_array(double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *sourceT, double *plume_trajX, double *plume_trajY, double *plume_trajZ, double *plume_trajR, double *plume_trajT){
	double r = 0.0;
	for(int j = 0; j < SDIM_FOR_FALL_CALC; j++){
		for(int i = 0; i < SDIM_FOR_PLUME_CALC + 1; i++){
			if((double)(j + 1) * S_DELTA_FOR_FALL_CALC == (double)i * S_DELTA_FOR_PLUME_CALC){
				sourceX[j] = plume_trajX[i-1]; sourceY[j] = plume_trajY[i-1]; sourceZ[j] = plume_trajZ[i-1]; sourceR[j] = plume_trajR[i-1]; sourceT[j] = plume_trajT[i-1];
				//printf("LINE481 i=%d\tj=%d\tr=0.0\n", i, j);
				break;
			}else if((double)(j + 1) * S_DELTA_FOR_FALL_CALC < (double)i * S_DELTA_FOR_PLUME_CALC && (double)(j + 1) * S_DELTA_FOR_FALL_CALC > (double)(i - 1) * S_DELTA_FOR_PLUME_CALC){
				r = (double)(j + 1) * S_DELTA_FOR_FALL_CALC - (double)(i - 1) * S_DELTA_FOR_PLUME_CALC;
				r = r / ((double)i * S_DELTA_FOR_PLUME_CALC - (double)(i - 1) * S_DELTA_FOR_PLUME_CALC);
				i--;
				sourceX[j] = plume_trajX[i - 1] + r * (plume_trajX[i] - plume_trajX[i - 1]);
				sourceY[j] = plume_trajY[i - 1] + r * (plume_trajY[i] - plume_trajY[i - 1]);
				sourceZ[j] = plume_trajZ[i - 1] + r * (plume_trajZ[i] - plume_trajZ[i - 1]);
				sourceR[j] = plume_trajR[i - 1] + r * (plume_trajR[i] - plume_trajR[i - 1]);
				sourceT[j] = plume_trajT[i - 1] + r * (plume_trajT[i] - plume_trajT[i - 1]);
				i++;
				//printf("LINE491 i=%d\tj=%d\tr=%1.4f\n", i, j, r);
				break;
			}
		}
	}
}

void writettlfallsummary(int phiint, double *ttlfalltime, double *ttlfalltimesummary){
	for(int z = 0; z < ZDIM; z++){
		ttlfalltimesummary[z + phiint * ZDIM] = ttlfalltime[z];
	}
}

void printfallsummary(const char *c, double *h, double *ttlfalltimesummary, double base, int stepnum, double stepdelta){
	FILE *outfile;
	outfile = fopen(c, "w");

	fprintf(outfile, "z\th(m)");
	for(int phiint = 0; phiint < stepnum; phiint++){
		fprintf(outfile, "\tphi%1.1f", base - phiint * stepdelta);
	}
	fprintf(outfile, "\n");
	for(int z = ZDIM - 1; z >= 0; z--){
		fprintf(outfile, "%d\t%1.1f", z, h[z]);
		for(int phiint = 0; phiint < stepnum; phiint++){
			fprintf(outfile, "\t%1.4f", ttlfalltimesummary[z + phiint * ZDIM]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}

void printsegsummary(const char *c, double *massreleased_perds_and_phidec, double base, int stepnum, double stepdelta){
	FILE *outfile;
	outfile = fopen(c, "w");

	fprintf(outfile, "i");
	for(int phiint = 0; phiint < stepnum; phiint++){
		fprintf(outfile, "\tphi%1.1f", base - phiint * stepdelta);
	}
	fprintf(outfile, "\n");
	for(int i = 0; i < SDIM_FOR_FALL_CALC; i++){
		fprintf(outfile, "%d", i);
		for(int phiint = 0; phiint < stepnum; phiint++){
			fprintf(outfile, "\t%1.4f", massreleased_perds_and_phidec[i + phiint * SDIM_FOR_FALL_CALC]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}

double confirm_released_mass(double *massreleased_per_ds_and_phidec){
	double totalofthefraction = 0.0;
	int s;
	
	for(int i = 0; i < SDIM_FOR_FALL_CALC * PHIDECDIM; i++){
		s = i % SDIM_FOR_FALL_CALC;
		if(s < SDIMCUTOFF){
			totalofthefraction += massreleased_per_ds_and_phidec[i];
		}else{
			totalofthefraction += 0;
		}
		
	}
	return(totalofthefraction);
}

void clearary(int dim, double *ary){
	for(int i = 0; i < dim; i++){
		ary[i] = 0.00;
	}
}

// F20
void drift_from_a_certain_source(double *source_x, double *source_y, double *source_height, double *sourceRadius, double *TotalFallTime, double *driftX, double *driftY, double *driftX_s, double *driftY_s, double *sigma_square){
	int s, z, phidec, idz, idz_s;
	int z_source;		// z_source means interval count of z axis of just above the source height
	double residue_up, fall_time_residue_up;
	double falltime, ttldriftX, ttldriftY;

	for(int idx = 0; idx < ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM; idx++){		// idx is count for driftXY_s and sigma_square
		z = idx % ZDIM;
		s = (idx / ZDIM) % SDIM_FOR_FALL_CALC;
		phidec = idx / (ZDIM * SDIM_FOR_FALL_CALC);
		z_source = ceil(source_height[s] / Z_DELTA); // source_height means source height
		if(z < z_source){
			idz = (phidec * ZDIM) + z; idz_s = (phidec * ZDIM) + z_source; // idz is count for driftXY and TotalFallTime
			residue_up = source_height[s] - (z_source - 1) * Z_DELTA;
			fall_time_residue_up =  (TotalFallTime[idz_s - 1] - TotalFallTime[idz_s]) * residue_up / Z_DELTA;

			falltime = TotalFallTime[idz] - TotalFallTime[idz_s - 1] + fall_time_residue_up;

			ttldriftX = (driftX[idz] - driftX[idz_s - 1]) + (driftX[idz_s - 1] - driftX[idz_s]) * residue_up / Z_DELTA;
			ttldriftY = (driftY[idz] - driftY[idz_s - 1]) + (driftY[idz_s - 1] - driftY[idz_s]) * residue_up / Z_DELTA;

			//printf("phidec=%d\ts=%d\tz=%d\tdrftX=%1.4f\tttldrftX = %1.4f\n", phidec, s, z, driftX[idz], ttldriftX);

			driftX_s[idx] = source_x[s] + ttldriftX;
			driftY_s[idx] = source_y[s] + ttldriftY;
			sigma_square[idx] = calc_sigma_square(sourceRadius[s] * PLUME_RADIUS_CORRECTION, falltime); // F21
		}
	}
} // End of the function (F20)

#ifdef CUDA

// Structures
typedef struct {
    float *lspmlD;
    float *ttlmlD;

    float *sourceZD;
    float *centXD;
    float *centYD;
    float *sigsqD;

    float *locXD;
    float *locYD;
    float *locZD;

    float *massreleasedD;
} DeviceBuffers;


typedef struct {
    float *ttlmlF;

    float *sourceZF;
    float *centXF;
    float *centYF;
    float *sigsqF;

    float *locXF;
    float *locYF;
    float *locZF;

    float *massreleasedF;
} HostBuffers;


typedef struct {
    size_t PSZC;
    size_t LSP;

    DeviceBuffers device;
    HostBuffers host;
} Buffers;
//

// FUNCTIONS FOR CUDA

/* 3. Device Buffer Allocation */
static void allocate_device_buffers_struct(Buffers *b)
{
    CUDA_CHECK(cudaMalloc((void**)&b->device.lspmlD,
        b->LSP * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.ttlmlD,
        LOCDIM * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.sourceZD,
        SDIMCUTOFF * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.centXD,
        b->PSZC * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.centYD,
        b->PSZC * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.sigsqD,
        b->PSZC * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locXD,
        LOCDIM * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locYD,
        LOCDIM * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locZD,
        LOCDIM * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.massreleasedD,
        SDIMCUTOFF * PHIDECDIM * sizeof(float)));
}
/* end of #3*/

/* 4. Host Data Packing: Fixed Data */
static void pack_fixed_data_buffers(
    Buffers *b,
    double *sourceZ,
    double *driftcentXs,
    double *driftcentYs,
    double *sigma_square,
    double *massreleased
){
	int phidec, s, z;  // indices for phi subdivision, source, and height level

	/*
	* CUDA device memory is handled as linear memory.
	* Therefore, multi-dimensional indices are flattened into
	* one-dimensional array indices before data transfer to GPU.
	*/
	int ipsz;  // flattened index for (phidec, s, z)
	int ips;   // flattened index for (phidec, s)

    for(int i = 0; i < SDIMCUTOFF; i++){
        b->host.sourceZF[i] = (float)sourceZ[i];
    }

    for(size_t ipszc = 0; ipszc < b->PSZC; ipszc++){
        phidec = (int)(ipszc / (SDIMCUTOFF * ZDIM));
        s      = (int)((ipszc / ZDIM) % SDIMCUTOFF);
        z      = (int)(ipszc % ZDIM);

        ipsz = (phidec * SDIM_FOR_FALL_CALC * ZDIM) + (s * ZDIM) + z;

        b->host.centXF[ipszc] = (float)driftcentXs[ipsz];
        b->host.centYF[ipszc] = (float)driftcentYs[ipsz];
        b->host.sigsqF[ipszc] = (float)sigma_square[ipsz];
    }

    for(int ipsc = 0; ipsc < SDIMCUTOFF * PHIDECDIM; ipsc++){
        s = ipsc % SDIMCUTOFF;
        phidec = ipsc / SDIMCUTOFF;

        ips = phidec * SDIM_FOR_FALL_CALC + s;

        b->host.massreleasedF[ipsc] = (float)massreleased[ips];
    }
}
/* End of 4*/

/* 5. Host Data Packing: Location Data */
static void pack_location_data_buffers(
    Buffers *b,
    double *locX,
    double *locY,
    double *locZ
){
    for(int j = 0; j < LOCDIM; j++){
        b->host.locXF[j] = (float)locX[j];
        b->host.locYF[j] = (float)locY[j];
        b->host.locZF[j] = (float)locZ[j];
    }
}
/* End of 5 */

/* 6. Copy Fixed Data to Device */
static void copy_fixed_data_to_device_buffers(Buffers *b)
{
	CUDA_CHECK(cudaMemcpy(
		b->device.sourceZD,
		b->host.sourceZF,
		SDIMCUTOFF * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.centXD,
		b->host.centXF,
		b->PSZC * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.centYD,
		b->host.centYF,
		b->PSZC * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.sigsqD,
		b->host.sigsqF,
		b->PSZC * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.massreleasedD,
		b->host.massreleasedF,
		SDIMCUTOFF * PHIDECDIM * sizeof(float),
		cudaMemcpyHostToDevice
	));
}
/* end of 6 */

/* 7. Copy Location Data to Device */
static void copy_location_data_to_device_buffers(Buffers *b)
{
    CUDA_CHECK(cudaMemcpy(
        b->device.locXD,
        b->host.locXF,
        LOCDIM * sizeof(float),
        cudaMemcpyHostToDevice
    ));

    CUDA_CHECK(cudaMemcpy(
        b->device.locYD,
        b->host.locYF,
        LOCDIM * sizeof(float),
        cudaMemcpyHostToDevice
    ));

    CUDA_CHECK(cudaMemcpy(
        b->device.locZD,
        b->host.locZF,
        LOCDIM * sizeof(float),
        cudaMemcpyHostToDevice
    ));
}
/* end of #7*/

/* 8. Kernel Launch */
static void launch_mass_loading_kernels_buffers(Buffers *b)
{
    if (b->LSP > INT_MAX) {
        fprintf(stderr, "Error: LSP too large: %zu\n", b->LSP);
        exit(EXIT_FAILURE);
    }

    int N = (int)b->LSP;

    int blocksize = 128;
    dim3 block(blocksize, 1, 1);
    dim3 grid((N + block.x - 1) / block.x, 1, 1);

    funcD01a<<<grid, block>>>(
        N, ZDIM, SDIMCUTOFF, PHIDECDIM, (float)Z_DELTA,
        b->device.lspmlD,
        b->device.ttlmlD,
        b->device.sourceZD,
        b->device.centXD,
        b->device.centYD,
        b->device.sigsqD,
        b->device.locXD,
        b->device.locYD,
        b->device.locZD,
        b->device.massreleasedD
    );
    CUDA_KERNEL_CHECK();

    dim3 grid2((LOCDIM + block.x - 1) / block.x, 1, 1);

    funcD01b<<<grid2, block>>>(
        N,
        LOCDIM,
        b->device.lspmlD,
        b->device.ttlmlD
    );
    CUDA_KERNEL_CHECK();

    CUDA_CHECK(cudaDeviceSynchronize());
}
/* end of #8*/

/* 9. Copy Result to Host */
static void copy_result_to_host_buffers(Buffers *b, double *ttlml)
{
    CUDA_CHECK(cudaMemcpy(
        b->host.ttlmlF,
        b->device.ttlmlD,
        LOCDIM * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    for(int i = 0; i < LOCDIM; i++){
        ttlml[i] = b->host.ttlmlF[i];
    }
}
/* end of #9 */

/* 10 Cleanup new */
static void cleanup_buffers(Buffers *b)
{
	free(b->host.locXF);
	free(b->host.locYF);
	free(b->host.locZF);
	free(b->host.ttlmlF);

	free(b->host.sourceZF);
	free(b->host.massreleasedF);
	free(b->host.centXF);
	free(b->host.centYF);
	free(b->host.sigsqF);

	CUDA_CHECK(cudaFree(b->device.locXD));
	CUDA_CHECK(cudaFree(b->device.locYD));
	CUDA_CHECK(cudaFree(b->device.locZD));

	CUDA_CHECK(cudaFree(b->device.sourceZD));
	CUDA_CHECK(cudaFree(b->device.centXD));
	CUDA_CHECK(cudaFree(b->device.centYD));
	CUDA_CHECK(cudaFree(b->device.sigsqD));
	CUDA_CHECK(cudaFree(b->device.massreleasedD));

	CUDA_CHECK(cudaFree(b->device.lspmlD));
	CUDA_CHECK(cudaFree(b->device.ttlmlD));
}
/* end of 10b */

/* New functions inserted on May 2, 2026*/
static void prepare_mass_loading(
    Buffers *b,
    double *sourceZ,
    double *driftcentXs,
    double *driftcentYs,
    double *sigma_square,
    double *massreleased,
    double *locX,
    double *locY,
    double *locZ
){
    /* device allocation */
    allocate_device_buffers_struct(b);

    /* fixed packing */
    pack_fixed_data_buffers(
        b,
        sourceZ,
        driftcentXs,
        driftcentYs,
        sigma_square,
        massreleased
    );

    /* location packing */
    pack_location_data_buffers(
        b,
        locX,
        locY,
        locZ
    );

    /* fixed copy */
    copy_fixed_data_to_device_buffers(b);
}

static void compute_mass_loading(
    Buffers *b,
    double *ttlml
)
{
    copy_location_data_to_device_buffers(b);
    launch_mass_loading_kernels_buffers(b);
    copy_result_to_host_buffers(b, ttlml);
}//end of the new functions



// D01a and D01b FOR GPU PROCESSING
// Calculate mass loading of a certain grain size on a certain point on the ground (Sloc) from a certain source: Sloc(phi, s)
void calc_mass_loading(double *sourceZ, double *driftcentXs, double *driftcentYs, double *sigma_square, double *locX, double *locY, double *locZ, double *lspml, double *ttlml, double *massreleased){
	/* This function prepares input arrays for GPU processing,
	* transfers them to GPU memory, launches CUDA kernels,
	* and retrieves the reduced total mass loading per location.
	*
	* On the GPU:
	*   funcD01a evaluates contributions per (location, phidec, source).
	*   funcD01b performs a reduction over (phidec, source)
	*           to obtain total mass loading per location.
	*
	* D in the name of parameter (e.g. ttlmlD) comes from "Device", which means such parameters
	* are used in GPU calculation
	* F in the name of parameter (e.g. ttlmlF) means such parameters are temporaly ones in CPU
	*/

	/* 1. Define size of arrays used in GPU */
	
	/* ---- dimensions and sizes --------------------------------------
	 * PSZC : size of arrays indexed by (phidec, source, height_interval)
	 *        e.g., driftcentXs, driftcentYs, sigma_square
	 *        PSZC stands Particle size, Souce and Z of a particle Cloud
	 *
	 * LSP  : size of arrays indexed by (location, phidec, source)
	 *        e.g., lspml
	 */
	size_t LSP;
	size_t PSZC; 
	
	/* ---- index definitions --------------------------------------
	 * phidec : phi (grain size) subdivision index
 	 * s      : source index along plume axis
 	 * z      : vertical layer index
	 */
	//int phidec, s, z;

	/* ---- flattened indices -----------------------------------------
	* Multi-dimensional indices (phidec, source, z) are mapped to
 	* 1D arrays for GPU memory access (CUDA global memory is linear).
	*
	* ipsz : index for (phidec, source, height_interval)
	* ips  : index for (phidec, source)
	*/
	//int ipsz;
	//int ips;

	PSZC = PHIDECDIM * SDIMCUTOFF * ZDIM;	//PSZC = PHIDECDIM * SDIM_FOR_FALL_CALC* ZDIM;
	LSP = (size_t)LOCDIM * SDIMCUTOFF * PHIDECDIM;	//LSP = LOCDIM * SDIM_FOR_FALL_CALC* PHIDECDIM;

	/* end of 1.*/

	/* 2. Host Buffer Allocation */

	float *ttlmlF;
	ttlmlF = (float *)malloc(LOCDIM * sizeof(float));
	if (!ttlmlF) {
    fprintf(stderr, "Error: malloc failed for ttlmlF\n");
    exit(EXIT_FAILURE);
	}

	// Allocate host buffers used for double-to-float conversion
	float *sourceZF, *centXF, *centYF, *sigsqF, *locXF, *locYF, *locZF, *massreleasedF;
	sourceZF = (float *)malloc(SDIMCUTOFF * sizeof(float));
	centXF = (float *)malloc(PSZC * sizeof(float));
	centYF = (float *)malloc(PSZC * sizeof(float));
	sigsqF = (float *)malloc(PSZC * sizeof(float));

	locXF = (float *)malloc(LOCDIM * sizeof(float));
	locYF = (float *)malloc(LOCDIM * sizeof(float));
	locZF = (float *)malloc(LOCDIM * sizeof(float));

	massreleasedF = (float *)malloc(PHIDECDIM * SDIMCUTOFF * sizeof(float));
	
	if (!ttlmlF || !sourceZF || !centXF || !centYF || !sigsqF ||
    !locXF || !locYF || !locZF || !massreleasedF) {
    fprintf(stderr, "Error: malloc failed for host buffers\n");
    exit(EXIT_FAILURE);
	}	

	/* end of 2.*/

	/* 3. Device Buffer Allocation */
	Buffers b;

	b.PSZC = PSZC;
	b.LSP  = LSP;

	/* host bridge: required before prepare */
	b.host.ttlmlF = ttlmlF;

	b.host.sourceZF = sourceZF;
	b.host.centXF   = centXF;
	b.host.centYF   = centYF;
	b.host.sigsqF   = sigsqF;

	b.host.locXF = locXF;
	b.host.locYF = locYF;
	b.host.locZF = locZF;

	b.host.massreleasedF = massreleasedF;

	prepare_mass_loading(
		&b,
		sourceZ,
		driftcentXs,
		driftcentYs,
		sigma_square,
		massreleased,
		locX,
		locY,
		locZ
	);
	compute_mass_loading(&b, ttlml);
	cleanup_buffers(&b);
	/* end of host bridge */

} // End of the function

__global__ void funcD01a(int N, int zdim, int sdim, int phidecdim, float zdelta, float *lspmlD, float *ttlmlD, float *sourceZD, float *centX, float *centY, float *sigma_square, float *locX, float *locY, float *locZ, float *massreleased){
		int j, s, z, phidec, ips, ipsz;
		float depcentX, depcentY, sigma2, square_distance;
		
		unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
		if(tid < N){
			j = tid / (phidecdim * sdim);
			s = (tid / phidecdim) % sdim;
			phidec = tid % phidecdim;
			z = locZ[j] / zdelta;

			ipsz = (phidec * sdim * zdim) + (s * zdim) + z;
			ips = s + phidec * sdim;

			depcentX = centX[ipsz+1] + (centX[ipsz] - centX[ipsz+1]) * (zdelta * (z + 1) - locZ[j]) / zdelta;
			//printf("depcentX = %1.1f\n", centX[0]);
			depcentY = centY[ipsz+1] + (centY[ipsz] - centY[ipsz+1]) * (zdelta * (z + 1) - locZ[j]) / zdelta;
			sigma2 = sigma_square[ipsz+1] + (sigma_square[ipsz] - sigma_square[ipsz+1]) * (zdelta * (z + 1) - locZ[j]) / zdelta;
			square_distance = pow((depcentX - locX[j]), 2) + pow((depcentY - locY[j]), 2);

			if(locZ[j] < sourceZD[s]){	//20241224
				lspmlD[tid] = 1 / (M_2PI * sigma2) * exp(-square_distance / (2 * sigma2)) * massreleased[ips]; // Original formulation appears in Bonadonna+ (2005)
#ifdef TEPHRA2
				lspmlD[tid] = 1 / (M_PI * sigma2) * exp(-square_distance / (sigma2)) * massreleased[ips];    // Formulation used in Tephra2 and WT
#endif
			//printf("%d\t%d\t%1.4e\t%1.4e\t%1.4e\n", s, phidec, massreleased[ips], lspmlD[tid], depcentX); //see also Line535
			}
			tid += blockDim.x * gridDim.x;
		}
}


__global__ void funcD01b(int N, int locdim, float *lspmlD, float *ttlmlD){
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if(tid < locdim){
		ttlmlD[tid] = 0.0;
		for(int i = 0; i < (N / locdim); i++){
			ttlmlD[tid] += lspmlD[tid * (N / locdim) + i];
		}
	}
}

// D01b
// Calculate mass loading on a certain grain size on a certain point on the ground.
// All grainsizes and sources are summed up from the "elements", which is calculated by D01a.
void cumulative_mass_x_phi(int phiint, double *lspml, double *tmp, double *ttl, double *cummassphi){
	double phi;

	for(int j = 0; j < LOCDIM; j++){
			phi = phiint + MAX_GRAINSIZE + 1;
			ttl[j] += tmp[j];
			cummassphi[j] += tmp[j] * phi;
	}
}
#endif

/* NON CUDA FUNCTIONS (START)*/
// D01a
// Calculate mass loading of a certain grain size on a certain point on the ground (Sloc) from a certain source: Sloc(phi, s)
void calc_mass_loading_element(int phisize, double *sourceZ, double *driftcentXs, double *driftcentYs, double *sigma_square, double *locX, double *locY, double *locZ, double *lspml, double *massreleased){
	int j, s, phidec, z;
	int ips;    // counter for arrays having grainsize(phidec) - plumelength(s; non cut off) order such as massreleased
	int ipsz;	// counter for arrays having grainsize(phidec) - plumelength(s; non cut off) - height(z) order such as driftX_s 
	double depcentX, depcentY, sigma2, square_distance;
	//char string[30];
	//FILE *outfile;
	
	//sprintf(string, "plume_fall_%dphi.txt", phisize);
	//outfile = fopen(string, "w");
	//fprintf(outfile, "idx\ti\tj\tphisize\tdepcentX\tdepcentY\tsourceZ\tdep-locX\tdep-locY\tsigma2\tsquare_distance\tmassloading\tsourcemagnitude\n");
	
	for(int idx = 0; idx < PHIDECDIM * SDIMCUTOFF * LOCDIM; idx++){ // idx is count for location(j) - plumelength(s cutoff) - grainsize(phidec) order
		j = idx / (PHIDECDIM * SDIMCUTOFF);
		s = (idx / PHIDECDIM) % SDIMCUTOFF;
		phidec = idx % PHIDECDIM;

		z = (int)(locZ[j] / Z_DELTA);
		ipsz = (phidec * SDIM_FOR_FALL_CALC * ZDIM) + (s * ZDIM) + z;
		ips = s + phidec * SDIM_FOR_FALL_CALC;

		if(locZ[j] < sourceZ[s]){
			depcentX = driftcentXs[ipsz+1] + (driftcentXs[ipsz] - driftcentXs[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			depcentY = driftcentYs[ipsz+1] + (driftcentYs[ipsz] - driftcentYs[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			sigma2 = sigma_square[ipsz+1] + (sigma_square[ipsz] - sigma_square[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			square_distance = pow((depcentX - locX[j]), 2) + pow((depcentY - locY[j]), 2);
			
			lspml[idx] = 1 / (M_2PI * sigma2) * exp(-square_distance / (2 * sigma2)) * massreleased[ips];
#ifdef TEPHRA2
			lspml[idx] = 1 / (M_PI * sigma2) * exp(-square_distance / (sigma2)) * massreleased[ips];	// Formulation used in Tephra2 and WT
#endif
			//if(j == 0 && phidec == 0){fprintf(outfile, "%d\t%d\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.6e\t%1.6e\t%1.6e\n", idx, s, j, phisize - phidec * 0.1, depcentX, depcentY, sourceZ[s], depcentX - locX[j], depcentY - locY[j], sigma2, square_distance, lspml[idx], massreleased[ips]);}
		}
	}
} // End of the function

// D01b
// Calculate mass loading on a certain grain size on a certain point on the ground.
// All grainsizes and sources are summed up from the "elements", which is calculated by D01a.
void calc_mass_loading_location(int phiint, double *lspml, double *massloading, double *ttl, double *cummassphi){
	int idx;
	double phi;

	for(int j = 0; j < LOCDIM; j++){
		for(idx = PHIDECDIM * SDIMCUTOFF * j; idx < PHIDECDIM * SDIMCUTOFF * (j + 1); idx++){
			phi = phiint + MAX_GRAINSIZE + 1; // - phidec * INTERVAL_DECIMAL_PHI;
			massloading[j] += lspml[idx];
			ttl[j] += lspml[idx];
			cummassphi[j] += lspml[idx] * phi;
		}
	}
	//printf("LINE302 j = %d\n", j);
}

/* NON CUDA FUNCTIONS (END)*/

void interval_fall_calc(int zmax, int phidecimal, double grainsize, double *h, double *atmP, double *atmT, double *windX, double *windY, double *driftX, double *driftY, double *ttlfalltime){
	// F10
	
#ifdef TEPHRA2
	double v0, falltime;
#else
	double v1, v0, falltime;
#endif

	//printf("\n\nh\tp\tt\tair_density\tair_viscosity\tRe\ttermfallv\n");
				// zmax is count for Ht
	v0 = calc_particle_terminal_velocity(h[zmax-1], grainsize, PUMICE_DENSITY, atmP[zmax-1], atmT[zmax-1]);
	
#ifdef TEPHRA2
	falltime = (h[zmax]-h[zmax-1]) / (v0);
	driftX[zmax-1] = falltime * (windX[zmax]);
	driftY[zmax-1] = falltime * (windY[zmax]);
#else
	v1 = calc_particle_terminal_velocity(h[zmax], grainsize, PUMICE_DENSITY, atmP[zmax], atmT[zmax]);	
	falltime = (h[zmax]-h[zmax-1]) / ((v1 + v0) / 2);
	driftX[zmax-1] = falltime * (windX[zmax-1] + windX[zmax]) / 2;
	driftY[zmax-1] = falltime * (windY[zmax-1] + windY[zmax]) / 2;
#endif
		
	ttlfalltime[zmax-1] = falltime;

	for(int z = zmax - 2; z >= 0; z--){ // The Loop 166
		v0 = calc_particle_terminal_velocity(h[z], grainsize, PUMICE_DENSITY, atmP[z], atmT[z]);
#ifdef TEPHRA2
		falltime = Z_DELTA / (v0);
		driftX[z] = driftX[z+1] + falltime * (windX[z+1]);
		driftY[z] = driftY[z+1] + falltime * (windY[z+1]);		
#else
		v1 = calc_particle_terminal_velocity(h[z+1], grainsize, PUMICE_DENSITY, atmP[z+1], atmT[z+1]);
		falltime = Z_DELTA / ((v1 + v0) / 2);
		driftX[z] = driftX[z+1] + falltime * (windX[z+1] + windX[z]) / 2;
		driftY[z] = driftY[z+1] + falltime * (windY[z+1] + windY[z]) / 2;
#endif
		ttlfalltime[z] = ttlfalltime[z+1] + falltime;
	}
}

// Calculate amount of particle segregation from plume
void mass_release_calc(int zmax, int phidecimal, double phi, double *h, double *atmP, double *atmT, double *windX, double *windY, double *massreleased){
	                               //no need of phidecimal: phi already has decimal number.
	double vphi, vw;
	double beta;
	double demon1, demon2;
	double pdf_fraction;
	double grainsize;

	grainsize = pow(2, -phi) * 0.001;

	pdf_fraction = calc_pdf_fraction(phi);

	vw = pow(windX[zmax], 2) + pow(windY[zmax], 2);
	vw = sqrt(vw);
	vphi = calc_particle_terminal_velocity(h[zmax], grainsize, PUMICE_DENSITY, atmP[zmax], atmT[zmax]);

	beta = vphi / (vw * PLUME_THICKNESS);
	
	//printf("%1.1f\t%1.4e\n", phi, vphi);

	for(int s = 0; s < SDIM_FOR_FALL_CALC; s++){
		if(s==0){
			demon1 = 0;
		}else{
			demon1 = -1 * beta * (s) * S_DELTA_FOR_FALL_CALC;
		}
		
		demon2 = -1 * beta * (s + 1) * S_DELTA_FOR_FALL_CALC;
		massreleased[s + phidecimal * SDIM_FOR_FALL_CALC] = ERUPTION_MASS * (exp(demon1) - exp(demon2)) * pdf_fraction;
		//printf("%1.4e\t%1.4e\n", demon1, demon2);
		//printf("%1.1f\t%d\t%1.4e\n", phi, s, massreleased[s + phidecimal * SDIM_FOR_FALL_CALC]);
	}
}

void locwrite(DEP *mll, double *locX, double *locY, double *locZ){
  for(int j = 0; j < LOCDIM ; j++){
    mll[j].x = locX[j];
		mll[j].y = locY[j];
		mll[j].z = locZ[j];
		mll[j].dist = sqrt(pow(locX[j], 2) + pow(locY[j], 2));
  }
}

void depwrite(int size, DEP *mll, double *loading){
	int phi;
	phi = size + MAX_GRAINSIZE + 1;
	
	for(int j = 0; j < LOCDIM; j++){
		mll[j].dep[size] = loading[j];
		if(phi > 0){mll[j].smallerthan1mm += loading[j];}	
	}
}


void sumwrite(DEP *mll, double *ttlmassloading, double *cummassphi){
	for(int j = 0; j < LOCDIM; j++){
		mll[j].ttlmassloading = ttlmassloading[j];
		if(ttlmassloading[j] > 0){
			mll[j].meandiameter = cummassphi[j] / ttlmassloading[j];
		}else{
			mll[j].meandiameter = -9999;
		}
	}
}

void printdeposit(DEP *mll){
	FILE *outfile;
	outfile = fopen("massloading.txt", "w");

	/*write header*/
	fprintf(outfile, "x(m)\ty(m)\tz(m)\tdistfromvent(m)");
	fprintf(outfile, "\tttlmassloading(kg/sq-m)\tF(percent)\tMean(phi)");
	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		fprintf(outfile, "\t%1.0fphi", MAX_GRAINSIZE + phiint + 1);
	}
	fprintf(outfile, "\n");

	/*write data*/
	for(int j = 0; j < LOCDIM; j++){
		fprintf(outfile, "%1.0f\t%1.0f\t%1.0f\t%1.1f", mll[j].x + VENT_EASTING, mll[j].y + VENT_NORTHING, mll[j].z, mll[j].dist);
		fprintf(outfile, "\t%1.4e\t%1.4f\t%1.4f", mll[j].ttlmassloading, mll[j].smallerthan1mm/mll[j].ttlmassloading*100, mll[j].meandiameter);
		for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
			fprintf(outfile, "\t%1.4e", mll[j].dep[phiint]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}

void atmosphere(int windlinenum, double *h, double *atmT, double *atmP, double *windX, double *windY, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres){
	double dir, v;
	double *vary, *dirary;
	vary = (double *)malloc(ZDIM * sizeof(double)), dirary = (double *)malloc(ZDIM * sizeof(double));
	
	//printf("windlinenum = %d\n", windlinenum);
	
	
	FILE *outfile;
	if(WRITE_COLUMN_FILES){
		outfile = fopen("atmosphere_used.txt", "w");
		fprintf(outfile, "z\th(m)\twind_dir\twind_v(m/s)\twindX(m/s)\twindY(m/s)\ttemp(K)\tpres(Pa)\n");		
	}

	for(int z = 0; z < ZDIM; z++){
		if (z == 0){
		atmT[z] = wind_tmp[0];
		atmP[z] = wind_pres[0] * 100; // hPa -> Pa
		v = wind_v[0];
		dir = wind_dir[0];
		vary[z] = v;
		dirary[z] = dir;
		windY[z] = v * cos(dir / 360 * 2 * M_PI);
		windX[z] = v * sin(dir / 360 * 2 * M_PI);
		}else{
		atmT[z] = calc_Tatm(h[z], windlinenum);
		atmP[z] = calc_Patm(h[z], windlinenum);
		v = get_V(h[z], windlinenum);
		dir = get_dir(h[z], windlinenum);
		vary[z] = v;
		dirary[z] = dir;
		windY[z] = v * cos(dir / 360 * 2 * M_PI);
		windX[z] = v * sin(dir / 360 * 2 * M_PI);
		}
		if(WRITE_COLUMN_FILES){fprintf(outfile, "%d\t%1.0f\t%1.0f\t%1.1f\t%1.1f\t%1.1f\t%1.1f\t%1.1f\n", z, h[z], dirary[z], vary[z], windX[z], windY[z], atmT[z], atmP[z]);}
	}
	free(vary); free(dirary);
	if(WRITE_COLUMN_FILES){fclose(outfile);}
}

double calc_sigma_square(double source_radius, double falltime) {
	// time needed for point source to diffuse until sigma equals to the plume radius
	double virtual_falltime = 0.0;
	double sigma_square = 0.0;

	if (falltime < FALL_TIME_THRESHOLD){
		// coarse particle	Bonadonna+(2005) Eq.6
		virtual_falltime = source_radius * source_radius / (4 * DIFFUSION_COEFFICIENT);
		sigma_square = 4 * DIFFUSION_COEFFICIENT * (falltime + virtual_falltime);
	}else{
		// fine particle	Bonadonna+(2005) Eq.8
		virtual_falltime = pow((5 * source_radius * source_radius) / (8 * EDDY_CONST), 0.4);
		sigma_square = 8 * EDDY_CONST / 5 * pow(falltime + virtual_falltime, 2.5);
	}
	return(sigma_square);
} // End of the function



int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;

  in_config = fopen(config_file, "r");

  if (in_config == NULL) {
    fprintf(stderr,
	    "Cannot open configuration file=[%s]:[%s]. Exiting.\n", config_file, strerror(errno));
    return 1;
  }

  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, in_config) != NULL) {
    /*fprintf(stderr, "%s\n", line); */
    if (line[0] == '#' || line[0] == '\n') continue;

    token = strtok_r(line, space, ptr1);
    if (!strncmp(token, "DIFFUSION_COEFFICIENT", strlen("DIFFUSION_COEFFICIENT"))) {
      token = strtok_r(NULL,space,ptr1);
      DIFFUSION_COEFFICIENT = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "DIFFUSION_COEFFICIENT=%.1f\n", DIFFUSION_COEFFICIENT);
	}
	else if (!strncmp(token, "EDDY_CONST", strlen("EDDY_CONST"))) {
      token = strtok_r(NULL,space,ptr1);
      EDDY_CONST = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "EDDY_CONST=%g\n", EDDY_CONST);
    }
	else if (!strncmp(token, "ENTRAIN_COEFF_KS", strlen("ENTRAIN_COEFF_KS"))) {
      token = strtok_r(NULL,space,ptr1);
      ENTRAIN_COEFF_KS = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "ENTRAIN_COEFF_KS=%g\n", ENTRAIN_COEFF_KS);
    }
	else if (!strncmp(token, "ENTRAIN_COEFF_KW", strlen("ENTRAIN_COEFF_KW"))) {
      token = strtok_r(NULL,space,ptr1);
      ENTRAIN_COEFF_KW = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "ENTRAIN_COEFF_KW=%g\n", ENTRAIN_COEFF_KW);
    }
    else if (!strncmp(token, "FALL_TIME_THRESHOLD", strlen("FALL_TIME_THRESHOLD"))) {
      token = strtok_r(NULL,space,ptr1);
      FALL_TIME_THRESHOLD = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "FALL_TIME_THRESHOLD=%.1f\n", FALL_TIME_THRESHOLD);
    }
    /*else if (!strncmp(token, "LITHIC_DENSITY", strlen("LITHIC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      LITHIC_DENSITY = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "LITHIC_DENSITY=%.1f\n", LITHIC_DENSITY);
    }*/
    else if (!strncmp(token, "PUMICE_DENSITY", strlen("PUMICE_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      PUMICE_DENSITY = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "PUMICE_DENSITY=%.1f\n", PUMICE_DENSITY);
    }
    else if (!strncmp(token, "Z_DELTA", strlen("Z_DELTA"))) {
      token = strtok_r(NULL, space, ptr1);
      Z_DELTA = (int)atoi(token);
      if(WRITE_CONF) fprintf(stderr, "Z_DELTA = %1.1f\n", Z_DELTA);
    }
    else if (!strncmp(token, "MINIMUM_CONTRIBUTION", strlen("MINIMUM_CONTRIBUTION"))) {
      token = strtok_r(NULL, space, ptr1);
      MINIMUM_CONTRIBUTION = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MINIMUM_CONTRIBUTION = %g\n", MINIMUM_CONTRIBUTION);
    }
    else if (!strncmp(token, "S_DELTA_FOR_PLUME_CALC", strlen("S_DELTA_FOR_PLUME_CALC"))) {
      token = strtok_r(NULL, space, ptr1);
      S_DELTA_FOR_PLUME_CALC = (int)atoi(token);
      if(WRITE_CONF) fprintf(stderr, "S_DELTA_FOR_PLUME_CALC = %1.1f\n", S_DELTA_FOR_PLUME_CALC);
    }
    else if (!strncmp(token, "S_DELTA_FOR_FALL_CALC", strlen("S_DELTA_FOR_FALL_CALC"))) {
      token = strtok_r(NULL, space, ptr1);
      S_DELTA_FOR_FALL_CALC = (int)atoi(token);
      if(WRITE_CONF) fprintf(stderr, "S_DELTA_FOR_FALL_CALC = %1.1f\n", S_DELTA_FOR_FALL_CALC);
    }
    else if (!strncmp(token, "ERUPTION_MASS", strlen("ERUPTION_MASS"))) {
      token = strtok_r(NULL, space, ptr1);
      ERUPTION_MASS = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "ERUPTION_MASS = %g\n", ERUPTION_MASS);
    }
    else if (!strncmp(token, "MAX_GRAINSIZE", strlen("MAX_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAX_GRAINSIZE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MAX_GRAINSIZE = %.0f\n", MAX_GRAINSIZE);
    }
    else if (!strncmp(token, "MIN_GRAINSIZE", strlen("MIN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MIN_GRAINSIZE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MIN_GRAINSIZE = %.0f\n", MIN_GRAINSIZE);
    }
    else if (!strncmp(token, "INTERVAL_DECIMAL_PHI", strlen("INTERVAL_DECIMAL_PHI"))) {
      token = strtok_r(NULL, space, ptr1);
      INTERVAL_DECIMAL_PHI = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "INTERVAL_DECIMAL_PHI = %.2f\n", INTERVAL_DECIMAL_PHI);
    }
    else if (!strncmp(token, "COLLAPSE_THEN_OFF", strlen("COLLAPSE_THEN_OFF"))) {
      token = strtok_r(NULL, space, ptr1);
      COLLAPSE_THEN_OFF = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "COLLAPSE_THEN_OFF = %d\n", COLLAPSE_THEN_OFF);
    }
    else if (!strncmp(token, "WRITE_DEPCENT_TRAJECTORY", strlen("WRITE_DEPCENT_TRAJECTORY"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_DEPCENT_TRAJECTORY = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_DEPCENT_TRAJECTORY = %d\n", WRITE_DEPCENT_TRAJECTORY);
    }
    else if (!strncmp(token, "WRITE_COLUMN_FILES", strlen("WRITE_COLUMN_FILES"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_COLUMN_FILES = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_COLUMN_FILES = %d\n", WRITE_COLUMN_FILES);
    }
    else if (!strncmp(token, "WRITE_DECIMAL_MASSLOADING", strlen("WRITE_COLUMN_FILES"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_DECIMAL_MASSLOADING = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_DECIMAL_MASSLOADING = %d\n", WRITE_DECIMAL_MASSLOADING);
    }
    else if (!strncmp(token, "WRITE_DECIMAL_FALL_TRAJ", strlen("WRITE_DECIMAL_FALL_TRAJ"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_DECIMAL_FALL_TRAJ = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_DECIMAL_FALL_TRAJ = %d\n", WRITE_DECIMAL_FALL_TRAJ);
    }
    else if (!strncmp(token, "WRITE_FALL_INFO_FILES", strlen("WRITE_FALL_INFO_FILES"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_FALL_INFO_FILES = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_FALL_INFO_FILES = %d\n", WRITE_FALL_INFO_FILES);
    }
    else if (!strncmp(token, "WRITE_CONF", strlen("WRITE_CONF"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_CONF = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_CONF = %d\n", WRITE_CONF);
    }
    else if (!strncmp(token, "WRITE_MASSLOADING", strlen("WRITE_MASSLOADING"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_MASSLOADING = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_MASSLOADING = %d\n", WRITE_MASSLOADING);
    }
    else if (!strncmp(token, "WRITE_SDIMCUTOFF", strlen("WRITE_SDIMCUTOFF"))) {
      token = strtok_r(NULL, space, ptr1);
      WRITE_SDIMCUTOFF = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "WRITE_SDIMCUTOFF = %d\n", WRITE_SDIMCUTOFF);
    }
    else if (!strncmp(token, "MEDIAN_GRAINSIZE", strlen("MEDIAN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MEDIAN_GRAINSIZE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MEDIAN_GRAINSIZE = %.2f\n", MEDIAN_GRAINSIZE);
    }
    else if (!strncmp(token, "MINIMUM_DEPOSIT_FOR_MD_CALC", strlen("MINIMUM_DEPOSIT_FOR_MD_CALC"))) {
      token = strtok_r(NULL, space, ptr1);
      MINIMUM_DEPOSIT_FOR_MD_CALC = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MINIMUM_DEPOSIT_FOR_MD_CALC = %.4f\n", MINIMUM_DEPOSIT_FOR_MD_CALC);
    }
    else if (!strncmp(token, "STD_GRAINSIZE", strlen("STD_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      STD_GRAINSIZE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "STD_GRAINSIZE = %.2f\n", STD_GRAINSIZE);
    }
    else if (!strncmp(token, "VENT_EASTING", strlen("VENT_EASTING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_EASTING = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "VENT_EASTING = %.1f\n", VENT_EASTING);
    }
    else if (!strncmp(token, "VENT_NORTHING", strlen("VENT_NORTHING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_NORTHING = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "VENT_NORTHING = %.1f\n", VENT_NORTHING);
    }
    else if (!strncmp(token, "VENT_ELEVATION", strlen("VENT_ELEVATION"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_ELEVATION = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "VENT_ELEVATION = %.1f\n", VENT_ELEVATION);
    }
    else if (!strncmp(token, "INITIAL_WATER_CONTENT", strlen("INITIAL_WATER_CONTENT"))) {
      token = strtok_r(NULL, space, ptr1);
      INITIAL_WATER_CONTENT = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "INITIAL_WATER_CONTENT = %.4f\n", INITIAL_WATER_CONTENT);
    }
    else if (!strncmp(token, "MAGMA_DISCHARGE_RATE", strlen("MAGMA_DISCHARGE_RATE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAGMA_DISCHARGE_RATE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MAGMA_DISCHARGE_RATE = %.1f\n", MAGMA_DISCHARGE_RATE);
    }
    else if (!strncmp(token, "MAGMA_TEMPERATURE", strlen("MAGMA_TEMPERATURE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAGMA_TEMPERATURE = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MAGMA_TEMPERATURE = %.1f\n", MAGMA_TEMPERATURE);
    }
    else if (!strncmp(token, "MESH_SIZE_IN_KM", strlen("MESH_SIZE_IN_KM"))) {
      token = strtok_r(NULL, space, ptr1);
      MESH_SIZE_IN_KM = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "MESH_SIZE_IN_KM = %.2f\n", MESH_SIZE_IN_KM);
    }
    else if (!strncmp(token, "INITIAL_PLUME_VELOCITY", strlen("INITIAL_PLUME_VELOCITY"))) {
      token = strtok_r(NULL, space, ptr1);
      INITIAL_PLUME_VELOCITY = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "INITIAL_PLUME_VELOCITY = %.1f\n", INITIAL_PLUME_VELOCITY);
    }
    else if (!strncmp(token, "VENT_RADIUS", strlen("VENT_RADIUS"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_RADIUS = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "VENT_RADIUS = %.1f\n", VENT_RADIUS);
    }
    else if (!strncmp(token, "S_MAX", strlen("S_MAX"))) {
      token = strtok_r(NULL, space, ptr1);
      S_MAX = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "S_MAX = %.1f\n", S_MAX);
    }
    else if (!strncmp(token, "PLUME_THICKNESS", strlen("PLUME_THICKNESS"))) {        /* added by Kaz 09-Mar-2020 */
      token = strtok_r(NULL, space, ptr1);
      PLUME_THICKNESS = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "PLUME_THICKNESS = %.1f\n", PLUME_THICKNESS);
    }
    else if (!strncmp(token, "PLUME_RADIUS_CORRECTION", strlen("PLUME_RADIUS_CORRECTION"))) {        /* added by Kaz 09-Mar-2020 */
      token = strtok_r(NULL, space, ptr1);
      PLUME_RADIUS_CORRECTION = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "PLUME_RADIUS_CORRECTION = %.1f\n", PLUME_RADIUS_CORRECTION);
    }
    else continue;
  }
  (void) fclose(in_config);
  return 0;
}

void printxyz(FILE *in, const char *header, int imax, double *sourceX, double *sourceY, double *sourceZ){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\n", i, sourceX[i], sourceY[i], sourceZ[i]);
  }
}

void printxyzq(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q, double *t){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", i, x[i], y[i], z[i], q[i], t[i]);
  }
}

void printxyze(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4e\n", i, x[i], y[i], z[i], q[i]);
  }
}


int get_line_number(FILE *f){
	char line[1000];
	int i = 0;
	while(NULL != fgets(line, 1000, f)){
		if(line[0] == '#')continue;
		i++;
	}
	return(i);
}

int get_wind_line_number(FILE *f){
	char line[1000];
	int i = 0;
	int additional = 1;
	int ret;
	double wind_height, wind_speed, wind_dir, wind_temp, wind_pres;
	while(NULL != fgets(line, 1000, f)){
		if(line[0] == '#')continue;
		else{
		while(ret=sscanf(line,
		"%lf %lf %lf %lf %lf",
		&wind_height,
		&wind_speed,
		&wind_dir,
		&wind_temp,
		&wind_pres), ret != 5){}
		}
		if(wind_height == 0 && i==0){
			additional = 0;
		}
		//printf("wind_height\t%1.1f\n", height[i]);
		i++;
	}
	return(i + additional);
}

void read_wind(FILE *f, double *height, double *speed, double *dir, double *atm_temp, double *atm_pres){
	char line[1000];
	int i = 0;
	int ret;
	double wind_height, wind_speed, wind_dir, wind_temp, wind_pres;
	while(NULL != fgets(line, 1000, f)){
		if(line[0] == '#')continue;
		else{
		while(ret=sscanf(line,
		"%lf %lf %lf %lf %lf",
		&wind_height,
		&wind_speed,
		&wind_dir,
		&wind_temp,
		&wind_pres), ret != 5){}
		}
		if(wind_height == 0 && i==0){
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}else if(wind_height > 0 && i==0){
			height[i] = 0;
			speed[i] = 0;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp + 0.0065 * wind_height;
			atm_pres[i] = wind_pres + 0.12 * wind_height;
			i++;
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}else{
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}
		//printf("wind_height\t%1.1f\n", height[i]);
		i++;
	}
}

void read_loc(FILE *f, double *x, double *y, double *z){
	char line[1000];
	int i = 0;
	double xtmp, ytmp, ztmp;
	int ret;
	//double wind_height, wind_speed, wind_dir, wind_temp, wind_pres;
	while(NULL != fgets(line, 1000, f)){
		if(line[0] == '#')continue;
		else{
		while(ret=sscanf(line,
		"%lf %lf %lf",
		&xtmp,
		&ytmp,
		&ztmp), ret != 3){}
		}
		x[i] = xtmp - VENT_EASTING; y[i] = ytmp - VENT_NORTHING; 
		if(ztmp < 0){z[i] = 0;}else{z[i] = ztmp;}
		if(xtmp-VENT_EASTING < MAPENDW){MAPENDW = xtmp-VENT_EASTING;}	//20241224
		if(xtmp-VENT_EASTING > MAPENDE){MAPENDE = xtmp-VENT_EASTING;}
		if(ytmp-VENT_NORTHING < MAPENDS){MAPENDS = ytmp-VENT_NORTHING;}
		if(ytmp-VENT_NORTHING > MAPENDN){MAPENDN = ytmp-VENT_NORTHING;}
		i++;
	}
	printf("MAP BOUNDARY W %1.0f\tE %1.0f\tS %1.0f\tN %1.0f\n", MAPENDW, MAPENDE, MAPENDS, MAPENDN);
}

void createisopachdata(DEP *l){
	qsort(l, LOCDIM, sizeof(DEP), compare_ttlmassloading);
	extractisopachdata(l);  // print out S-A relation

	//qsort(l, LOCDIM, sizeof(DEP), compare_Md);
	//countmeandiameter(l);
}

int compare_ttlmassloading(const void * a, const void * b){
    double z1 = ((DEP *)a)->ttlmassloading;
	double z2 = ((DEP *)b)->ttlmassloading;

	if (z1 < z2) {
        return 1;
    } else {
        return -1;
    }
}

int compare_Md(const void * a, const void * b){
    double z1 = ((DEP *)a)->meandiameter;
	double z2 = ((DEP *)b)->meandiameter;

	if (z1 > z2) {
        return 1;
    } else {
        return -1;
    }
}

void extractisopachdata(DEP *l){
	int jmax = 0;
	int i = 0;
	int imax;
	int j = 0;
	int sqrtA, sqrtAinit;
	
	double maxsqrtA;
	
	double distance = -9999;
	double dir, fvalue;
	double poweroftwo, ratio;
	double log2S;
	double demon1, demon2;
	double minS = 0.01;
	double maxS;

	double intS, intA, intSqrtA, intx, inty, intdist, intdistax, intmean, intF;
	double intxprevious, intyprevious, intdistaxprevious;

	FILE *outfile;

	for(j=0; j<LOCDIM; j++){
		if(l[j].x > MAPENDW && l[j].x < MAPENDE && l[j].y > MAPENDS && l[j].y < MAPENDN){
					if(l[j].dist > distance){
						distance = l[j].dist;
						jmax++;
					} 
		}
	}
	
	maxS = l[0].ttlmassloading;
	
	DEP *l2;			// l2 stores array of distribution axis, which is defined as most distal point above the certain thickness.
	l2 = (DEP *)calloc(jmax, sizeof(DEP));
	int j2 = 0;
	distance = -9999;
	outfile = fopen("S_vs_Area.txt", "w");
	fprintf(outfile, "massloading(kg/m2)\tarea(sqkm)\tsqrtA(km)\tdist_from_vent(m)\tx(m)\ty(m)\tdir(deg)\tMean(phi)\tF(percent)\n");

	for(j=0; j<LOCDIM; j++){
		if(l[j].x > MAPENDW && l[j].x < MAPENDE && l[j].y > MAPENDS && l[j].y < MAPENDN){
			if(l[j].dist > distance){
				l2[j2].j = j + 1; //20250112
				l2[j2].x = l[j].x;
				l2[j2].y = l[j].y;
				l2[j2].z = l[j].z;
				dir = calc_dir(l2[j2].x, l2[j2].y);
				l2[j2].dist = l[j].dist;
				l2[j2].ttlmassloading = l[j].ttlmassloading;
				l2[j2].smallerthan1mm = l[j].smallerthan1mm;
				l2[j2].meandiameter = l[j].meandiameter;
				l2[j2].dep[0] = l2[j2].j * MESH_SIZE_IN_KM * MESH_SIZE_IN_KM; // area of isopach 20250111
				l2[j2].dep[1] = sqrt(l2[j2].dep[0]); // square root area
				fvalue = l2[j2].smallerthan1mm / l2[j2].ttlmassloading * 100;
				
				// distance from the vent along the axis
				if(j2==0){
				l2[j2].dep[3] = sqrt(l2[j2].x * l2[j2].x + l2[j2].y * l2[j2].y);
				}else{
				demon1 = l2[j2-1].x - l2[j2].x; demon2 = l2[j2-1].y - l2[j2].y;
				l2[j2].dep[3] = sqrt(demon1 * demon1 + demon2 * demon2) + l2[j2-1].dep[3];
				}
				fprintf(outfile, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.4f\n", l2[j2].ttlmassloading, l2[j2].dep[0], l2[j2].dep[1], l2[j2].dist, l2[j2].x, l2[j2].y, dir, l2[j2].meandiameter, fvalue);
				distance = l[j].dist;
				maxsqrtA = l2[j2].dep[1];
				j2++;
			} 
		}else{break;}
	}
	fclose(outfile);
	
	
  outfile = fopen("S_vs_Area_summary.txt", "w");
  log2S = floor(log2(l[0].ttlmassloading)); //20250128
  if(l2[j2-1].ttlmassloading > minS){minS = l2[j2-1].ttlmassloading;} //20250128
  imax = (int)(log2S) - (int)ceil(log2(minS)) + 1; //20250105

  fprintf(outfile, "massloading(kg/m2)\tarea(sqkm)\tsqrtA(km)\tdist_from_vent(m)\tx(m)\ty(m)\tdir(deg)\tMean(phi)\tF(percent)\n");
  for(i=0; i < imax; i++){
  	for(j=0; j<jmax; j++){
        //poweroftwo = pow(2, log2S - float(i)); //20250105
  		poweroftwo = pow(2, log2S - (double)i); //20250105
  		if(l2[j].ttlmassloading >= poweroftwo && l2[j+1].ttlmassloading < poweroftwo){
  			ratio = (poweroftwo - l2[j+1].ttlmassloading) / (l2[j].ttlmassloading - l2[j+1].ttlmassloading);
  			intA = (l2[j].dep[0] - l2[j + 1].dep[0]) * ratio + l2[j + 1].dep[0];
  			intSqrtA = (l2[j].dep[1] - l2[j + 1].dep[1]) * ratio + l2[j + 1].dep[1];
  			intx = (l2[j].x - l2[j + 1].x) * ratio + l2[j + 1].x;
  			inty = (l2[j].y - l2[j + 1].y) * ratio + l2[j + 1].y;
  			dir = calc_dir(l2[j].x, l2[j].y);
  			intdist  = sqrt(intx * intx + inty * inty);
  			 
  			if(j==0){
  			intdistax = sqrt(intx * intx + inty * inty);
  			intdistaxprevious = intdistax; intxprevious = intx; intyprevious = inty;
  			}else{
  			demon1 = intx - intxprevious;
  			demon2 = inty - intyprevious;
  			intdistax = sqrt(demon1 * demon1 + demon2 * demon2) + intdistaxprevious;
  			intdistaxprevious = intdistax; intxprevious = intx; intyprevious = inty;
  			} 			
  			
  			intmean = (l2[j].meandiameter - l2[j + 1].meandiameter) * ratio + l2[j + 1].meandiameter;
  			intF =  (l2[j].smallerthan1mm / l2[j].ttlmassloading - l2[j+1].smallerthan1mm / l2[j+1].ttlmassloading) * ratio + l2[j+1].smallerthan1mm / l2[j+1].ttlmassloading;
  			intF = 100 * intF;
  			//printf("%1.1f\t%1.1f\t%1.1f\n", l2[j].ttlmassloading, poweroftwo, l2[j+1].ttlmassloading);
  			fprintf(outfile, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.4f\n", poweroftwo, intA, intSqrtA, intdist, intx, inty, dir, intmean, intF);
  			break;
  		}
  	}
  	if(j==jmax-1){break;}
  }
  fclose(outfile);
  
  // 20250129
  outfile = fopen("DF.txt", "w");

  fprintf(outfile, "MaxS(kg/m2)\t0.1MaxS\t0.01MaxS\tD\tF\n");
  fprintf(outfile, "%1.4f\t%1.4f\t%1.4f\t", maxS, 0.1*maxS, 0.01*maxS);
  	for(j=0; j<jmax; j++){
  		if(l2[j].ttlmassloading >= 0.1*maxS && l2[j+1].ttlmassloading < 0.1*maxS){
  		ratio = (0.1*maxS - l2[j+1].ttlmassloading) / (l2[j].ttlmassloading - l2[j+1].ttlmassloading);
  		intF =  (l2[j].smallerthan1mm / l2[j].ttlmassloading - l2[j+1].smallerthan1mm / l2[j+1].ttlmassloading) * ratio + l2[j+1].smallerthan1mm / l2[j+1].ttlmassloading;
  		}
  		
  		if(l2[j].ttlmassloading >= 0.01*maxS && l2[j+1].ttlmassloading < 0.01*maxS){
  		ratio = (0.01*maxS - l2[j+1].ttlmassloading) / (l2[j].ttlmassloading - l2[j+1].ttlmassloading);
  		intA = (l2[j].dep[0] - l2[j + 1].dep[0]) * ratio + l2[j + 1].dep[0];
  		break;
  		}
  	}
  fprintf(outfile, "%1.4f\t%1.4f\n", intA, intF);
  fclose(outfile);
  
  
  // 20250129
  outfile = fopen("Area_vs_S_summary.txt", "w");
  sqrtAinit = (int)ceil(l2[0].dep[1]);
  if(l2[j2-1].ttlmassloading > minS){minS = l2[j2-1].ttlmassloading;}
  imax = (int)(log2S) - (int)ceil(log2(minS)) + 1;
  maxsqrtA = ceil(maxsqrtA);

  fprintf(outfile, "massloading(kg/m2)\tarea(sqkm)\tsqrtA(km)\tdist_from_vent(m)\tx(m)\ty(m)\tdir(deg)\tMean(phi)\tF(percent)\n");
  
  for(sqrtA=sqrtAinit; sqrtA<int(maxsqrtA); sqrtA++){
	for(j=0; j<j2; j++){
		if(l2[j].dep[1] <= double(sqrtA) && l2[j+1].dep[1] > double(sqrtA)){
			ratio = ((double)(sqrtA) - l2[j].dep[1]) / (l2[j+1].dep[1] - l2[j].dep[1]);
			intA = (l2[j+1].dep[0] - l2[j].dep[0]) * ratio + l2[j].dep[0];
			intS = (l2[j+1].ttlmassloading - l2[j].ttlmassloading) * ratio + l2[j].ttlmassloading;
			intx = (l2[j + 1].x - l2[j].x) * ratio + l2[j].x;
			inty = (l2[j + 1].y - l2[j].y) * ratio + l2[j].y;
			dir = calc_dir(l2[j].x, l2[j].y);
			intdist  = sqrt(intx * intx + inty * inty);
			 
			if(j==0){
			intdistax = sqrt(intx * intx + inty * inty);
			intdistaxprevious = intdistax; intxprevious = intx; intyprevious = inty;
			}else{
			demon1 = intx - intxprevious;
			demon2 = inty - intyprevious;
			intdistax = sqrt(demon1 * demon1 + demon2 * demon2) + intdistaxprevious;
			intdistaxprevious = intdistax; intxprevious = intx; intyprevious = inty;
			} 			
			
			intmean = (l2[j + 1].meandiameter - l2[j].meandiameter) * ratio + l2[j].meandiameter;
			intF =  (l2[j + 1].smallerthan1mm / l2[j + 1].ttlmassloading - l2[j].smallerthan1mm / l2[j].ttlmassloading) * ratio + l2[j].smallerthan1mm / l2[j].ttlmassloading;
			intF = 100 * intF;
			//printf("%1.1f\t%1.1f\t%1.1f\n", l2[j].ttlmassloading, poweroftwo, l2[j+1].ttlmassloading);
			fprintf(outfile, "%1.4e\t%1.1f\t%1.1f\t%1.4f\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.4f\n", intS, intA, double(sqrtA), intdist, intx, inty, dir, intmean, intF);
			break;
		}
	}	// end of j loop
  } 		// end of sqrtA loop
  fclose(outfile);



}

void countmeandiameter(DEP *l){
	int count=0;
	int flag=0;
	double distance = 0.0;
	double phi;
	double x, y;
	double dir;

	FILE *outfile, *outfile2;
	outfile = fopen("Mean_vs_Area_summary.txt", "w");
	outfile2 = fopen("Mean_vs_Area_all.txt", "w");

	fprintf(outfile, "#Mean(phi)\tarea(sqkm)\tmax_dist(m)\tx_of_max_dist(m)\ty_of_max_dist(m)\tdir_of_max_dist(deg)\n");
	fprintf(outfile2, "#Mean(phi)\tarea(sqkm)\tmax_dist(m)\n");

	for(int j=0; j<LOCDIM; j++){
		if(l[j].ttlmassloading >= MINIMUM_DEPOSIT_FOR_MD_CALC){
			if(flag == 0){
				phi = ceil(l[j].meandiameter * 10) * 0.1;
				distance = l[j].dist;
				flag = 1;
			}
			if(distance < l[j].dist){
				distance = l[j].dist; x = l[j].x; y = l[j].y;
				dir = calc_dir(x, y);
			}
			if(phi > 5){break;}
			if(l[j].meandiameter > phi){
				fprintf(outfile, "%1.4f\t%1.4f\t%1.1f\t%1.1f\t%1.1f\t%1.1f\n", phi, (count-1) * MESH_SIZE_IN_KM * MESH_SIZE_IN_KM, distance, x, y, dir);
				//phi++
				phi = phi + 0.1;
			}
			count++;
			fprintf(outfile2, "%1.4f\t%1.4f\t%1.1f\n", l[j].meandiameter, count * MESH_SIZE_IN_KM * MESH_SIZE_IN_KM, distance);
		}	
	}
	fclose(outfile);
	fclose(outfile2);
}

double calc_dir(double x, double y){	// calculate direction of local point from vent
	double dir;
	double plusdir = 0.0;
	if(y < 0){
		plusdir = 180;
	}else if(x < 0){
		plusdir = 360;
	}else{
		plusdir = 0;
	}
	//printf("x = %1.1f\ty = %1.1f\tplusdir = %1.1f\n", x, y, plusdir);
 	dir = atan(x / y) / (M_2PI) * 360 + plusdir;
 	return(dir);
}

void obtaintheoreticallyreleased(RELEASE *r){
	double phi;
	double pdf_fraction = 0.0;

	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		for(int phidecimal = 0; phidecimal < PHIDECDIM; phidecimal++){
			phi = phiint + MAX_GRAINSIZE + 1 - phidecimal * INTERVAL_DECIMAL_PHI;
			pdf_fraction += calc_pdf_fraction(phi);
		}// END OF DECIMAL PHI LOOP
		r[phiint].theoretical = pdf_fraction * ERUPTION_MASS;
		pdf_fraction=0.0;
	}
}

void printparticlereleased(RELEASE *r){
	FILE *outfile;

	outfile = fopen("particle_released.txt", "w");
	fprintf(outfile, "#i\tFraction(phi)\tTheoretical(kg)\tWtPercent\tActual(kg)\tActual/Theoretical\n");

	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		fprintf(outfile, "%d\t%1.1f\t%1.4e\t%1.4f\t%1.4e\t%1.4f\n", phiint, r[phiint].phi, r[phiint].theoretical, r[phiint].theoretical/ERUPTION_MASS*100, r[phiint].actual, r[phiint].actual/r[phiint].theoretical);
	}

	fclose(outfile);
}

void printsegregation_per_ds(SEG *r){
	FILE *outfile;

	outfile = fopen("segregation_per_ds.txt", "w");
	
	fprintf(outfile, "#source");
	for(int i = 0; i < MIN_GRAINSIZE - MAX_GRAINSIZE; i++){
		fprintf(outfile, "\t%1.0f", i + MAX_GRAINSIZE + 1);
	}
	fprintf(outfile, "\n");
	
	for(int s = 0; s < SDIM_FOR_FALL_CALC; s++){
		fprintf(outfile, "%d", s);
		for(int phiint = 0; phiint < MIN_GRAINSIZE - MAX_GRAINSIZE; phiint++){
			fprintf(outfile, "\t%1.4e", r[s].mass_from_ds[phiint]);
		}
		fprintf(outfile, "\n");
	}

	fclose(outfile);
}

//// ORIGINALLY IN GRAIN.C

// F11
double calc_particle_terminal_velocity(double h, double ashdiam, double part_density, double p, double t) {
	// double h is almost useless. Just shows h of calculated P & T condition
  // Modified from function “particle_fall_time” in tephra_calc.c
	double air_density, air_viscosity, temp;
 	double vtl, vti, vtt;
 	double reynolds_number;
 	double particle_terminal_velocity;
	double gravity = 9.81;

	// p in Pa here. Change to hPa
  air_density = p * 0.0034837 / t;  // US STD Atomosphere 1976 P 15. Eq. 42
  air_viscosity = 1.458e-6 * pow(t, 1.5) / (t + 110.4);	// US STD Atomosphere 1976 P 19. Eq. 51

#ifdef TEPHRA2
	// air density and viscosity in the previous versions such as Tephra2 and WT
  	air_density = 1.293 * exp(-h / 8200);
  	air_viscosity = 0.000018325;
#endif
	/*  Based on Bonadonna and Phillips (2003) JGR 108, 2034. Eq. A4
    	vtl is terminal velocity (m/s) in laminar regime Re < 6
    	vti is terminal velocity (m/s) in intermediate regime 6 <Re <500
    	vtt is terminal velocity (m/s) in turbulent regime Re > 500*/
  	vtl = gravity * ashdiam * ashdiam * (part_density - air_density) / (18 * air_viscosity);
  	reynolds_number = ashdiam * air_density * vtl / air_viscosity;
  	particle_terminal_velocity = vtl;

	if (reynolds_number >= 6.0) {
    		temp = 4 * gravity * gravity * (part_density - air_density) * (part_density - air_density) / (225 * air_density * air_viscosity);
    		vti = ashdiam * pow(temp, 1.0 / 3.0);
    		reynolds_number = ashdiam * air_density * vti / air_viscosity;
    		particle_terminal_velocity = vti;

  		if (reynolds_number >= 500.0) {
    			vtt = sqrt( 3.1 * gravity * ashdiam * (part_density - air_density) / air_density);
    			reynolds_number =  ashdiam * air_density * vtt / air_viscosity;
    			particle_terminal_velocity = vtt;
  		}
  	}
  return particle_terminal_velocity;
}

void phiconvert(){
	if(MIN_GRAINSIZE < MAX_GRAINSIZE){
		int tmp = MIN_GRAINSIZE;
		MIN_GRAINSIZE = MAX_GRAINSIZE;
		MAX_GRAINSIZE = tmp;
	}
}


double calc_pdf_fraction(double phi){				// calculate fraction of the particle size phi
																						// based on particle distribution function
	double demon1, demon2, demon3;
	double frac;
	//                 1                 -(x-myu)^2
	// f(x) = ------------------   exp ---------------
	//        sqrt(2pi * sigma^2)         2*sigma^2

	demon1 = sqrt(2 * M_PI * pow(STD_GRAINSIZE, 2));
	demon2 = pow((phi - INTERVAL_DECIMAL_PHI / 2 - MEDIAN_GRAINSIZE), 2);
	demon3 = 2 * pow(STD_GRAINSIZE, 2);

	frac = 1 / demon1 * exp(-1 * demon2 / demon3) * INTERVAL_DECIMAL_PHI;
	return(frac);
}

//// ORIGINALLY IN WINDY.C
// Global

// atmosphere structure

double g_dir;

double Ra =  285;
double Rg0 =  462;
//double g =  9.81;


double t0 =  293;	//293;
double x, north, east;	// horizontal position; x means max length

double p = 100000;	//101325.0;

double ds;
double dz;
double U = -9999;
double R = -9999;
double n0;

double rho_s = 1200;
double rho_w = 1000;

double Ca = 998; //713;
double Cs = 1617; //1100;
double Cv = 1850; //1850

double theta = M_PI / 2; // PI / 2

// param in func11
double E, Cp, Cp0;

// param in func12 - 15
double M, Q, rho_c, rho_a, Ue, V;

// param in func16
int flag=0; //	0, gas-thrust;   1, buoyant;     2, umbrella
			//  rho_a > rho_c    rho_a < rho_c   rho_a > rho_c

// param in func17-19
double n, Q0, Rg;


double gz, gs;
double ta, dp_over_dz;
double smax;

//////////////////

double windy(int imax, double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *taftervent, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres){	// The main routine in this file
	int i = 0;
	int total = imax; // total line number of wind file

	double Hg = -9999, Hb = -9999, Ht = -9999;
	double U0, R0;

	double z0;
	double T;
	double time_after_vent = 0.0;
	double gz_previous, x_previous, g_dir_previous, north_previous;
	double east_previous, ta_previous, p_previous, rho_a_previous, rho_c_previous, n_previous;
	double Q_previous, Cp_previous, Rg_previous, V_previous, M_previous; //theta_previous;
	double U_previous, R_previous, T_previous;

	n0 = INITIAL_WATER_CONTENT;
	z0 = VENT_ELEVATION;
    ds = S_DELTA_FOR_PLUME_CALC;

	gz= z0;
	gs = 0;
	FILE *f, *f2;

	T = MAGMA_TEMPERATURE;

	makewindstruct(imax, wind_alt, wind_v, wind_dir, wind_tmp, wind_pres);

	if(WRITE_COLUMN_FILES) f = fopen("plume.txt", "w");

	// initialize

	ta = calc_Tatm(gz, total);
	p = calc_Patm(gz, total);
	Cp0 = calcCp0();	// func 20.5

	rho_a = func23(p, ta);
	n=n0;

	Rg=Rg0;
	rho_c=func17(n0, p, Rg, T); // get rho_c

	if(n0 < 0 || n0 > 1){
      fprintf(stderr,
  	      "ERROR\nYou need proper INITIAL_WATER_CONTENT in config file\nPROGRAM HAS BEEN HALTED\n\n");
      exit(1);
	}

	if(MAGMA_DISCHARGE_RATE < 0 || INITIAL_PLUME_VELOCITY < 0 || VENT_RADIUS < 0){
		if(MAGMA_DISCHARGE_RATE < 0 && INITIAL_PLUME_VELOCITY > 0 && VENT_RADIUS > 0){
			Q = rho_c * INITIAL_PLUME_VELOCITY * VENT_RADIUS * VENT_RADIUS;
			U = INITIAL_PLUME_VELOCITY;
			R = VENT_RADIUS;
		}else if(MAGMA_DISCHARGE_RATE > 0 && INITIAL_PLUME_VELOCITY < 0 && VENT_RADIUS > 0){
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2013)
			U = Q / (rho_c * VENT_RADIUS * VENT_RADIUS);
			R = VENT_RADIUS;
		}else if(MAGMA_DISCHARGE_RATE > 0 && INITIAL_PLUME_VELOCITY > 0 && VENT_RADIUS < 0){
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2013)
			U = INITIAL_PLUME_VELOCITY;
			R = sqrt(Q / (rho_c * INITIAL_PLUME_VELOCITY));
		}else{
	        fprintf(stderr,
	    	      "ERROR\nYou need to assign at least two parameters properly from U, R and Q in the config file\nPROGRAM HAS BEEN HALTED 179\n\n");
	        exit(1);
		}
	}else{
        fprintf(stderr,
    	      "ERROR\nYou need to assign at least two parameters properly from U, R and Q in the config file\nPROGRAM HAS BEEN HALTED 184\n\n");
        exit(1);
	}


	// initialize (func 11)
	// Q = rho_c * U * R * R;
	M = rho_c * U * U * R * R;
	E = Q * Cp0 * T;

	Q0 = Q;
	U0 = U;
	R0 = R;

	Cp = Cp0;

	V = get_V(gz, total);		// wind velocity
	g_dir = get_dir(gz, total);	// wind direction
	x = 0.0;
	north = 0.0;
	east = 0.0;
    //sourceX[i] = east; sourceY[i] = north; sourceZ[i] = gz, sourceR[i] = R, taftervent[i] = time_after_vent;

  // Calculate plume parameters until reaching Hb: See while loop after Line 259
  if(WRITE_COLUMN_FILES){
  	fprintf(f, "#z\ts\tx\tdir\tnorthing\teasting\tTa\tP\tatm_dens\tcol_dens\tn\t");
		fprintf(f, "Q\tCp\tRg\tV\tUe\tM\ttheta\t");
		fprintf(f, "U\tR\tTm\tTime\n");
  }

	if(S_MAX < 0){smax = 99999;}else{smax = S_MAX;}

	//while(i < 100000 && M > 0.0 && theta > 0.0 && gs <= smax){	// Till 2023.08.08
	while(i < SDIM_FOR_PLUME_CALC && M > 0.0 && theta > 0.0 && gs <= smax){
		//Ht = gz; // when M < 0 (static) or theta < 0 (windy), z just before the height is considered as Ht
		// top of the gas thrust region
		if(flag==0 && rho_a - rho_c > 0){flag=1; Hg = gz;}	// top of the gas-thrust region
		// top of the convective region
		if(flag==1 && rho_a - rho_c < 0){flag=2; Hb = gz;}	// top of the convective region
		if(WRITE_COLUMN_FILES){
			T = E / Q / Cp;
			fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", gz, gs, x, g_dir, north, east, ta, p, rho_a, rho_c, n);
			fprintf(f, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Q, Cp, Rg, V, Ue, M, theta);
			fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4f\n", U, R, T, time_after_vent);
		}

		gz_previous = gz; x_previous = x; g_dir_previous = g_dir; north_previous = north;
		east_previous = east; ta_previous = ta; p_previous = p; rho_a_previous = rho_a; rho_c_previous = rho_c; n_previous = n;
		Q_previous = Q; Cp_previous = Cp; Rg_previous = Rg; V_previous = V; M_previous = M; //theta_previous = theta;
		U_previous = U; R_previous = R; T_previous = T;
		
		rk(total, T);
		time_after_vent += ds / ((U_previous + U) / 2);
		taftervent[i] = time_after_vent;
		
		if(theta > 0.0){ //No i increment before here means i = 0 is not at crater but at next step after the crater
			sourceX[i] = east; sourceY[i] = north; sourceZ[i] = gz, sourceR[i] = R;
		}else{
			gz = gz_previous; 
			Ht = gz;

			x = x_previous; g_dir = g_dir_previous; north = north_previous;
			east = east_previous; ta = ta_previous; p = p_previous; rho_a = rho_a_previous; rho_c = rho_c_previous; n = n_previous;
			Q = Q_previous; Cp = Cp_previous; Rg = Rg_previous; V = V_previous; Ue = 0; M = M_previous; theta = 0;
			U = U_previous; R = R_previous; T = T_previous;
			i--;
			//sourceX[i] = east; sourceY[i] = north; sourceZ[i] = Ht, sourceR[i] = R;	
		}
		//time_after_vent += ds / ((U_previous + U) / 2);
		//taftervent[i] = time_after_vent;
		//printf("L1652i = %d\n", i);
		i++;
	}
	
	// Processing when plume reached Ht
	if(WRITE_COLUMN_FILES) {
		f2 = fopen("plume_parameters.txt", "w");
		fprintf(f2, "#Q0\tU0\tR0\tHg\tHb\tHt\tR@Ht\tColumnT\tAtmT\n");
		fprintf(f2, "%1.4e\t%1.4e\t%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", Q0 * M_PI, U0, R0, Hg, Hb, Ht, R, T, ta);
		fclose(f2);
	}

	if(Hb == -9999 && COLLAPSE_THEN_OFF){
		if(WRITE_COLUMN_FILES) fclose(f);
		printf("Plume collapsed. No tephra dispersal calculated.\n");
		exit(1);
	}

	PLUME_HEIGHT = Ht;
    gs = gs - ds;

	// Print out plume parameters after reaching Hb
	V = get_V(Ht, total);
	while (i < SDIM_FOR_PLUME_CALC){
		gs = gs + ds;
		x += ds * cos(theta);
		north = north + ds * cos(g_dir / 360 * 2 * M_PI);
		east = east + ds * sin(g_dir / 360 * 2 * M_PI);
		if(WRITE_COLUMN_FILES){
			fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Ht, gs, x, g_dir, north, east, ta, p, rho_a, rho_c, n);
			fprintf(f, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Q, Cp, Rg, V, Ue, M, theta);
			fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4f\n", V, R, T, time_after_vent);
		}
    	//printf("i = %d\n", i);
		sourceX[i] = east; sourceY[i] = north; sourceZ[i] = Ht, sourceR[i] = R;
		time_after_vent += ds / V;
		taftervent[i] = time_after_vent;
		i++;
	}

	if(WRITE_COLUMN_FILES) fclose(f);
	//printf("kokodayo %1.4f\t%d\n", S_MAX, SDIM_FOR_FALL_CALC);
	//read_plume_file(SDIM_FOR_FALL_CALC);

	return(Ht);
}

void rk(int total, double T){
	//double dp_over_ds, dQ_over_ds, dM_over_ds, dtheta_over_ds, dE_over_ds;
	//double dp_over_ds1, dQ_over_ds1, dM_over_ds1, dtheta_over_ds1, dE_over_ds1;
	//double dp_over_ds2, dQ_over_ds2, dM_over_ds2, dtheta_over_ds2, dE_over_ds2;
	//double dp_over_ds3, dQ_over_ds3, dM_over_ds3, dtheta_over_ds3, dE_over_ds3;

	double dQ_over_ds, dM_over_ds, dtheta_over_ds, dE_over_ds;
	double dQ_over_ds1, dM_over_ds1, dtheta_over_ds1, dE_over_ds1;
	double dQ_over_ds2, dM_over_ds2, dtheta_over_ds2, dE_over_ds2;
	double dQ_over_ds3, dM_over_ds3, dtheta_over_ds3, dE_over_ds3;
	double dQ_over_ds4, dM_over_ds4, dtheta_over_ds4, dE_over_ds4;
	double dx;

	double E_tmp, M_tmp, rho_a_tmp, rho_c_tmp, p_tmp, theta_tmp, Q_tmp;

	E_tmp = E; M_tmp = M; rho_a_tmp=rho_a; rho_c_tmp=rho_c; theta_tmp = theta; Q_tmp = Q;


	/////////////////////////////////
	// k1     ///////////////////////
	//dp_over_ds1 = func22(p, ta);
	dQ_over_ds1 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds1 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds1 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds1 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds1);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds1 * ds * 0.5;
	M_tmp = M + dM_over_ds1 * ds * 0.5;
	theta_tmp = theta + dtheta_over_ds1 * ds * 0.5;
	E_tmp = E + dE_over_ds1 * ds * 0.5;

	// centre position of the next step
	dz = ds * sin(theta_tmp) * 0.5;


	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(gz+ dz, total);
	ta = calc_Tatm(gz+ dz, total);
	rho_a_tmp = func23(p_tmp, ta);

	Cp = func20(n);					// calc Cp
	//printf("k1\n");
	V  = get_V(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k2     ///////////////////////
	//dp_over_ds2 = func22(p_tmp, ta);
	dQ_over_ds2 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds2 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds2 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds2 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds1);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds2 * ds * 0.5;
	M_tmp = M + dM_over_ds2 * ds * 0.5;
	theta_tmp = theta + dtheta_over_ds2 * ds * 0.5;
	E_tmp = E + dE_over_ds2 * ds * 0.5;

	// centre position of the next step
	dz = ds * 0.5 * sin(theta_tmp);


	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(gz+ dz, total);
	ta = calc_Tatm(gz+ dz, total);
	rho_a_tmp = func23(p_tmp, ta);

	Cp = func20(n);					// calc Cp
	//printf("k2\n");
	V = get_V(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k3     ///////////////////////
	//dp_over_ds3 = func22(p_tmp, ta);
	dQ_over_ds3 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds3 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds3 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds3 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds2);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds3 * ds;
	M_tmp = M + dM_over_ds3 * ds;
	theta_tmp = theta + dtheta_over_ds3 * ds;
	E_tmp = E + dE_over_ds3 * ds;

	// centre position of the next step
	dz = ds * sin(theta_tmp);


	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(gz+ dz, total);
	ta = calc_Tatm(gz+ dz, total);
	rho_a_tmp = func23(p_tmp, ta);

	Cp = func20(n);					// calc Cp
	//printf("k3\n");
	V = get_V(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k4     ///////////////////////
	//dp_over_ds4 = func22(p_tmp, ta);
	dQ_over_ds4 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds4 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds4 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds4 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds3);

	////////////////////////////////////
	////////////////////////////////////
	//// set new step value   //////////

	//dp_over_ds = 	 (dp_over_ds1 + 2 * dp_over_ds2 + 2 * dp_over_ds3 + dp_over_ds4)/6;
	dQ_over_ds = 	 (dQ_over_ds1 + 2 * dQ_over_ds2 + 2 * dQ_over_ds3 + dQ_over_ds4)/6;
	dM_over_ds = 	 (dM_over_ds1 + 2 * dM_over_ds2 + 2 * dM_over_ds3 + dM_over_ds4)/6;
	dtheta_over_ds = (dtheta_over_ds1 + 2 * dtheta_over_ds2 + 2 * dtheta_over_ds3 + dtheta_over_ds4)/6;
	dE_over_ds =     (dE_over_ds1 + 2 * dE_over_ds2 + 2 * dE_over_ds3 + dE_over_ds4)/6;


	//printf("U = %1.1f\tdQ_over_ds=%1.4f\tdM_over_ds=%1.4f\n", U, dQ_over_ds, dM_over_ds);

	// generate next step parameters
	Q = Q + dQ_over_ds * ds;
	M = M + dM_over_ds * ds;
	if(M<0){M=0;}
	theta = theta + dtheta_over_ds * ds;
	E = E + dE_over_ds * ds;

	// centre position of the next step
	dx = ds * cos(theta);
	dz = ds * sin(theta);
	north = north + dx * cos(g_dir / 360 * 2 * M_PI);
	east = east + dx * sin(g_dir / 360 * 2 * M_PI);
	x = x + dx;
	gs= gs+ ds;
	gz= gz+ dz;

	// atmosphreic content of the next step
	n = func18(Q);				// calc n
	p = calc_Patm(gz, total);
	ta = calc_Tatm(gz, total);
	rho_a = func23(p, ta);

	Cp = func20(n);					// calc Cp
	//printf("k4\n");
	V = get_V(gz, total);
	g_dir = get_dir(gz, total);



	T = E / Q / Cp;

	Rg = func19(n);					// calc Rg
	rho_c = func17(n, p, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M, Q, theta_tmp, V);
	U = M / Q;
	R = sqrt(Q / (U * rho_c));
}

double func12(double M_tmp, double rho_a_tmp, double rho_c_tmp, double Q_tmp){			// plume mass flux
	double dQ_over_ds;

	dQ_over_ds = 2 * rho_a_tmp * Ue * Q_tmp / sqrt(rho_c_tmp * M_tmp);

	//printf("%1.4f\n", dQ_over_ds);
	return dQ_over_ds;
}

double func13(double rho_a_tmp, double rho_c_tmp, double M_tmp, double theta_tmp, double Q_tmp){
	double dM_over_ds;

	dM_over_ds = GRAVITY * (rho_a_tmp - rho_c) * Q_tmp * Q_tmp / (rho_c_tmp * M_tmp) * sin(theta_tmp);
	dM_over_ds = dM_over_ds + 2 * rho_a_tmp * Q_tmp / sqrt(rho_c_tmp * M_tmp) * Ue * V * cos(theta_tmp);

	return dM_over_ds;
}

double func14(double M_tmp, double Q_tmp, double rho_a_tmp, double rho_c_tmp, double theta_tmp){
	double dtheta_over_ds;

	dtheta_over_ds = GRAVITY * (rho_a_tmp - rho_c_tmp) * Q_tmp * Q_tmp * cos(theta_tmp) / (rho_c_tmp * M_tmp * M_tmp);
	dtheta_over_ds = dtheta_over_ds - 2 * rho_a_tmp * Q_tmp * Ue * V * sin(theta_tmp) / (M * sqrt(rho_c_tmp * M_tmp));

	//printf("dtheta = %1.4f\n", dtheta_over_ds);
	return dtheta_over_ds;
}

double func15(double M_tmp, double Q_tmp, double rho_a_tmp, double rho_c_tmp, double Ta, double theta_tmp, double dQ_over_ds){
	double dE_over_ds;
	double term1, term2, term3, term4;

	term1 = (Ca * Ta + Ue * Ue / 2) * dQ_over_ds;
	term2 = M_tmp * M_tmp / (2 * Q_tmp * Q_tmp) * dQ_over_ds;
	term3 = rho_a_tmp / rho_c_tmp * Q_tmp * GRAVITY * sin(theta_tmp);
	term4 = 2 * rho_a_tmp * Ue * V * cos(theta_tmp) * sqrt(M_tmp / rho_c_tmp);
	dE_over_ds = term1 + term2 - term3 - term4;

	return dE_over_ds;
}

double func16(double M_tmp, double Q_tmp, double theta_tmp, double V_tmp){
	double ue_tmp;
	double ks_tmp;

	//if(flag==0){ks_tmp=sqrt(rho_a/rho_c)/16;}
	//else{ks_tmp=ks;}	// use these lines when you use ks for gas thrust region; include rho_a and rho_c as local

	ks_tmp=ENTRAIN_COEFF_KS;	// gas thrust region also assumes 0.09
				// remove this line when you take
				// ks = f(rho_a. rho_c)

	//printf("flag=%d\tks=%1.4f\n", flag, ks_tmp);

	ue_tmp = ks_tmp * fabs(M_tmp/Q_tmp - V_tmp * cos(theta_tmp)) + ENTRAIN_COEFF_KW * fabs(V_tmp * sin(theta_tmp));

	return ue_tmp;
}

double func17(double n_tmp, double p_tmp, double Rg_tmp, double T_tmp){			// plume density
	double rho_tmp;

	rho_tmp = (1 - n_tmp) / rho_s + n_tmp * Rg_tmp * T_tmp / p_tmp;
	rho_c = 1 / rho_tmp;
	//printf("rho     = %1.4f\n", rho);
	return rho_c;
}

double func18(double Q_tmp){			// solid content in the plume
	double n_tmp;
	n_tmp = 1 - (1 - n0) * Q0 / Q_tmp;

	return n_tmp;
}

double func19(double n_tmp){
	double Rg_tmp;
	Rg_tmp = Ra + (Rg0 - Ra) * n0 * (1 - n_tmp) / (n_tmp * (1 - n0));

	return Rg_tmp;
}

double func20(double n_tmp){
	double Cp_tmp;
	Cp_tmp = Ca + (Cp0 - Ca) * (1 - n_tmp) / (1 - n0);

	return Cp_tmp;
}

double calcCp0(){
	double Cp0_tmp;
	Cp0_tmp = n0 * Cv + (1 - n0) * Cs;
	return Cp0_tmp;
}


double calc_Tatm(double h, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double t_atm=9999.9;

	t_atm=W1[total-1].t_atm;

	while(i<total){
		if(h < W1[i].wind_height){
			t_atm = W1[i-1].t_atm + (W1[i].t_atm - W1[i-1].t_atm) * (h - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			break;
		}
		i++;
	}
	//printf("%d\t%1.4f\t%1.4f\n", i, h, v);
	return t_atm;
}

double calc_Patm(double h, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double t_atm0=9999.9, t_atm=9999.9;
	double p_atm0=9999.9, p_atm=9999.9;
	double a;

	t_atm=W1[total-1].t_atm;
	t_atm0=W1[total-1].t_atm;
	p_atm0=W1[total-1].p_atm;

	while(i<total){
		if(h < W1[i].wind_height){
			t_atm = W1[i-1].t_atm + (W1[i].t_atm - W1[i-1].t_atm) * (h - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			t_atm0 = W1[i-1].t_atm;
			p_atm0 = W1[i-1].p_atm;
			break;
		}
		i++;
	}

	a = (h - W1[i-1].wind_height) * GRAVITY * 2 / (Ra * (t_atm + t_atm0));
	p_atm = p_atm0 / exp(a) * 100; // hPa -> Pa
	//printf("h = %1.4f\ta_0 = %1.4f\tta_0 = %1.4f\tta_1 = %1.4f\te = %1.4f\tp = %1.4f\tp0 = %1.4f\n", z, a, t_atm0, t_atm, exp(a), p_atm, p_atm0);
	return p_atm;
}

double func22(double p_tmp, double t_tmp){	// atmospheric pressure
	double dp_over_ds;

	dp_over_ds = -1 * (GRAVITY * p_tmp) / (Ra * t_tmp);
	return dp_over_ds;
}

double func23(double p_tmp, double t_tmp){	// atmospheric density
	double rho_tmp;

	rho_tmp = p_tmp / (Ra * t_tmp);
	return rho_tmp;
}


double get_V(double h, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double v=9999.9;

	v=W1[total-1].wind_speed;

	while(i<total){
		if(h < W1[i].wind_height){
			v = W1[i-1].wind_speed + (W1[i].wind_speed - W1[i-1].wind_speed)*(h - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			break;
		}
		i++;
	}
	//printf("%d\t%1.4f\t%1.4f\n", i, h, v);
	return v;
}

double get_dir(double h, int total){	//return wind direction based on discrete wind data
	int i=1;
	double dir=9999.9;
	double ratio, wind1, wind2;

	dir=W1[total-1].wind_dir;

	while(i<total){
		wind1 = W1[i-1].wind_dir;
		wind2 = W1[i].wind_dir;
		//printf("h = %1.4f\t, w = %1.4f\n", h, W1[i].wind_height);
		if(h < W1[i].wind_height){
			ratio = (h - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			if(wind2 - wind1 > 180){
				wind1 = wind1 + 360;
			}else if(wind1 - wind2 > 180){
				wind2 = wind2 + 360;
			}

			dir = ratio * (wind2 - wind1) + wind1;

			if(dir>360){dir=dir-360;}
			if(dir<0){dir=dir+360;}

			break;
		}
		i++;
	}
	//printf("dir\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", i, z, ratio, dir, wind1, wind2);
	return dir;
}

void makewindstruct(int imax, double *alt, double *v, double *dir, double *temp, double *pres){
	W1 = (WIND *)malloc((imax) * sizeof(WIND));

	for(int i = 0; i < imax; i++){
    W1[i].day=0;
    W1[i].hour=0;
    W1[i].wind_height=alt[i];
    W1[i].wind_speed=v[i];
    W1[i].wind_dir=dir[i];
		W1[i].t_atm=temp[i];
		W1[i].p_atm=pres[i];
	}

	/*for(int i = 0; i < imax; i++){
		printf("%d\t%1.4f\t%1.4f\n", i, W1[i].wind_dir, W1[i].wind_speed);
	}*/
}


void writetrajectory(int phiint, double *x, double *y, double *sig, double *r, SEG *seg_ds){
		int idz; //(phidec * sdim * zdim) + (s * zdim) + z
		int phi;
		double released;	// mass released in the s interval of the size fraction
		char string[20];
		FILE *outfile;

		phi = phiint + MAX_GRAINSIZE + 1; 
		//printf("maxgrainsize = %1.0f phiint = %d\n", MAX_GRAINSIZE, phiint);
		if(WRITE_DEPCENT_TRAJECTORY){
			sprintf(string, "depcenttraj%d.txt", phi);
			outfile = fopen(string, "w");
			fprintf(outfile, "i\tx_from_vent\ty_from_vent\tx_coord\ty_coord\tsigma_square\tmass_released\n");
		}

		for(int i = 0; i < SDIM_FOR_FALL_CALC; i++){
			released = 0.0;
			for(int gsize = 0; gsize < PHIDECDIM; gsize++){
				released += r[i + gsize * SDIM_FOR_FALL_CALC];
			}
			
			seg_ds[i].mass_from_ds[phiint] = released; 
			
			idz = i * ZDIM;
			if(WRITE_DEPCENT_TRAJECTORY) fprintf(outfile, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4e\n", i, x[idz], y[idz], x[idz] + VENT_EASTING, y[idz] + VENT_NORTHING, sig[idz], released);
		}
		if(WRITE_DEPCENT_TRAJECTORY) fclose(outfile);
}

void get_sdimcutoff(double *sigma_square, SEG *massreleased_str, int phiint){
	double averagedmassloading;
	//int phisize;
	
	//phisize = phiint + MAX_GRAINSIZE + 1;
	//printf("s\tsegregate\tsigmasquare\n");
	for(int s = 0; s < SDIM_FOR_FALL_CALC; s++){
		averagedmassloading = massreleased_str[s].mass_from_ds[phiint] / sigma_square[s * ZDIM];
		//printf("%d\t%d\t%1.4e\t%1.4e\t%1.4e\n", phisize, s, massreleased_str[s].mass_from_ds[phiint], sigma_square[s * ZDIM], averagedmassloading);
		if(averagedmassloading < MINIMUM_CONTRIBUTION * S_DELTA_FOR_FALL_CALC){SDIMCUTOFF = s; break;}
		//if(massreleased_str[s].mass_from_ds[phiint] / S_DELTA_FOR_FALL_CALC < MINIMUM_CONTRIBUTION){SDIMCUTOFF = s; break;}
	}
}

void printttlml(int phiint, double *ttlml){
	int j, phi;
	char string[20];
	FILE *outfile;
	
	phi = phiint + MAX_GRAINSIZE + 1; 
	sprintf(string, "ttlml%d.txt", phi);
	outfile = fopen(string, "w");

	fprintf(outfile, "j\tttlml\n");
	for(int i = 0; i < LOCDIM; i++){
		fprintf(outfile, "%d\t%1.4e\n", j, ttlml[i]);
	}
	fclose(outfile);
}

void printlspml(int phiint, double *lspml){
	int j, s, phidec, phi;
	char string[20];
	FILE *outfile;
	
	phi = phiint + MAX_GRAINSIZE + 1; 
	sprintf(string, "lspmlD%d.txt", phi);
	outfile = fopen(string, "w");

	fprintf(outfile, "j\ts\tphi\tphidec\tlspml\n");
	for(int i = 0; i < LOCDIM * SDIM_FOR_FALL_CALC * PHIDECDIM; i++){
		phidec = i % PHIDECDIM;
		s = (i / PHIDECDIM) % SDIM_FOR_FALL_CALC;
		j = i / (PHIDECDIM * SDIM_FOR_FALL_CALC);
		fprintf(outfile, "%d\t%d\t%1.1f\t%d\t%1.4e\n", j, s, phi - phidec * 0.1, phidec, lspml[i]);
	}
	fclose(outfile);
}

void printS4eachphidec(DEP *l, int phiint, double *lspml){ // Output S for 0.1 phi interval at each location (2024/12/3)
	int idx, idxmax, j, jprevious, phidec, phi;
	double massloadingeachphidec[PHIDECDIM];
	char string[20];
	FILE *outfile;
	
	phi = phiint + MAX_GRAINSIZE + 1;
	idxmax = PHIDECDIM * SDIMCUTOFF * LOCDIM;
	idx = 0;
	j = 0;
	
	// create file and output header
	sprintf(string, "S_decimalphi_%d.txt", phi);
	outfile = fopen(string, "w");
	fprintf(outfile, "X\tY\tZ");
	
	while(idx < PHIDECDIM){
		fprintf(outfile, "\t%1.1f", phi - idx * 0.1);
		idx++;
	}
	fprintf(outfile, "\n");
	// end of creating file
	
	idx = 0;
	jprevious = 0;
	
	while(idx < idxmax + 1){
		j = idx / (PHIDECDIM * SDIMCUTOFF);
		phidec = idx % PHIDECDIM;
		
		if(j == jprevious){
			massloadingeachphidec[phidec] += lspml[idx];
		}else{
			phidec = 0;
			fprintf(outfile, "%1.1f\t%1.1f\t%1.1f", l[j].x + VENT_EASTING, l[j].y + VENT_NORTHING, l[j].z);
			while(phidec < PHIDECDIM){
				fprintf(outfile, "\t%1.4e", massloadingeachphidec[phidec]);
				massloadingeachphidec[phidec] = 0.0;
				phidec++;
			}
			fprintf(outfile, "\n");
		}
		jprevious = j;
		idx++;
	}
}
