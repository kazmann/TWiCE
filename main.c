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
    /*CUDA_CHECK(cudaDeviceSynchronize());*/                           \
} while (0)
// end of CUDA Check

#define TEPHRA2
//#define CUDA
//#define TEST	//OUTPUT massloading_loc_source_phi.txt which is massloading contribution for each Location, particle Source, Phi in decimal

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
int PRINT_PROGRESS = 0;
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

// Plume trajectory
double S_MAX= 100000;					// Length of Plume
int SDIM_FOR_PLUME_CALC = -9999;		// Number of segment along the plume for trajectory calculation
int SDIM_FOR_FALL_CALC = -9999;			// Number of particle source along the plume
int SDIMCUTOFF;							// Released mass below this threshold is treated as zero and not considered a source.
double S_DELTA_FOR_PLUME_CALC = 10;		// Step size for Runge-Kutta calculation for plume trajectory
double S_DELTA_FOR_FALL_CALC = 100;		// Length of a plume segment treated as a single source along the plume axis
double Z_DELTA = 100;					// Vertical step size for fall calculation
int ZDIM;   							// number of step (Z interval) for fall calculation

// Locations of the ground
int LOCDIM; 							// number of locations to calc




/*
 * Global state variables for plume calculation.
 *
 * These variables are mainly used by plume_calculation() and
 * advance_plume_state_rk4(). They represent the current plume state,
 * atmospheric conditions, and model constants during plume integration.
 *
 * Note:
 * These are kept global to avoid passing a large number of tightly coupled
 * plume-state variables through the RK4 integration functions.
 */

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

// Entrainment coefficient for plume calculation
double ENTRAIN_COEFF_KS = 0.09;  // another k should be introduced for gas thurst region but uniform value in this code
double ENTRAIN_COEFF_KW = 0.9;   // See Woodhouse et al. (2012) https://doi.org/10.1029/2012JB009592Digital Object Identifier (DOI)


// MAP BOUNDARY
double MAPENDE = -9999999;
double MAPENDW =  9999999;
double MAPENDS =  9999999;
double MAPENDN = -9999999;

/*
 * DEP: data structure for each ground location
 * used to store data for final result file named massloading.txt
 *
 * j                  : index of location
 * x, y, z            : coordinates [m]
 * dist               : distance from vent [m]
 * ttlmassloading     : total mass loading at this location [kg/m^2]
 * smallerthan1mm     : mass fraction of particles smaller than 1 mm
 * meandiameter       : mean grain size at this location
 * dep[phi]           : mass loading per phi class
 *                      (special use in some modes: area of isopach, etc.)
 */
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
} DEP;	//location_properties or l[j]

// Total released particle mass for a phi interval
typedef struct {
    double phi;
    double theoretical; // From particle size distribution function
    double actual;      // May be lower due to S_MAX limitation,
                        // especially for small particles
} RELEASE;

// Particle mass for sub-phi intervals (e.g., 0.1 phi),
// distributed along the plume axis
typedef struct {
    double mass_from_ds[20]; // Number of size classes in phi scale
} SEG;

//High altitude meteorology
typedef struct {
  int day;
  int hour;
  double wind_height;	// height a.s.l. in km
  double wind_speed;	// the average windspeed in m/s
  double wind_dir;		// average wind direction in +/- degrees from north
  double t_atm;	// atmospheric temperature
  double p_atm;	// atmospheric pressure
} WIND;
static WIND *W1;

#define MAX_LINE 200
#define M_2PI 2*M_PI

//prototypes
int init_globals(char *config_file);
void write_particle_release_theoretical_vs_actual(RELEASE *r);
void write_massrelease_per_source_phi(SEG *r);
void printxyz(FILE *in, const char *type, int i, double *srcX, double *srcY, double *srcZ);
void printxyzq(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q, double *t);
void printxyze(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q);
void write_vertical_profiles_for_phi(const char *c, double *h, double *ttlfalltime_phiint, double base, int stepnumber, double stepdelta); // modified on 2024.07.28
void write_phi_s_table(const char *c, double *massreleased_per_ds_and_phidec, double base, int stepnumber, double stepdelta);


void write_total_massloading(int phiint, double *ttlml);
void write_massloading_loc_source_phi(int phiint, double *massloading_loc_source_phi);
void write_massloading_per_phidec_at_locations(DEP *l, int phiint, double *massloading_loc_source_phi);

int get_line_number(FILE *in_wind);
int get_wind_line_number(FILE *in_wind);
void get_sdimcutoff(double *cloud_sigma2, SEG *massreleased_per_ds_and_phidec, int phiint);

void interpolate_atmosphere_and_wind(int windlinenum, double *h, double *atmT, double *atmP, double *windX, double *windY, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);
void compute_falltime_and_drift_profile(int zmax, int phidecimal, double grainsize, double *h, double *atmP, double *atmT, double *windX, double *windY, double *driftX, double *driftY, double *ttlfalltime);
void store_profile_for_phi(int phiint, double *ttlfalltime, double *ttlfalltime_phiint);
void drift_from_a_certain_source(double *source_x, double *source_y, double *source_height, double *sourceRadius, double *TotalFallTime, double *driftX, double *driftY, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2);
void write_cloud_trajectory_and_mass(int phiint, double *cloud_center_x, double *cloud_center_y, double *sigma_squre, double *massreleased, SEG *seg);
double calc_cloud_sigma2(double source_radius, double falltime);
void compute_mass_release_along_plume(int zmax, int phidecimal, double phi, double *h, double *atmP, double *atmT, double *windX, double *windY, double *massreleased);
double compute_total_released_mass(double *massreleased);
void set_coordinates_to_location_properties(DEP *l, double *locX, double *locY, double *locZ);
void store_massloading_for_phi(int size, DEP *l, double *massloading);
void store_total_massloading_and_mean_phi(DEP *l, double *ttlmassloading, double *cummassphi);
void write_deposit_summary(DEP *location_properties);	/* write massloading.txt */
void clear_array(int dim, double *ary);
int compare_ttlmassloading(const void *a, const void *b);
int compare_Md(const void * a, const void * b);
void set_source_points_on_plume(double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *sourceT, double *plume_trajX, double *plume_trajY, double *plume_traj_Z, double *plume_trajR, double *plume_trajT);

void createisopachdata(DEP *l);
void generate_isopach_analysis(DEP *l);
void countmeandiameter(DEP *l);
double compute_direction_from_vent(double x, double y);

void compute_theoretical_particle_release(RELEASE *r);


void read_wind_profile(FILE *f, double *h, double *v, double *d, double *t, double *p);
void read_locations_and_convert_to_vent_centered(FILE *f, double *x, double *y, double *z);

double plume_calculation(int linenum, double *sourceX, double *sourceY, double *sourceZ, double *sourceRadius, double *timeaftervent, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);
void advance_plume_state_rk4(int, double); // plume calculation using Runge-Kutta
void makewindstruct(int imax, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres);

double func12(double, double, double, double);
double func13(double, double, double, double, double);
double func14(double, double, double, double, double);
double func15(double, double, double, double, double, double, double);
double func16(double, double, double, double);

double func17(double, double, double, double);				// calc plume density
double func18(double);				// particle content
double func19(double);				// Rg calc
double calc_plume_heat_capacity(double);

double calc_Cp0(void);			// Cp0 calc (eq. 20.5; in text between eq20 and 21)
double calc_Tatm(double, int);
double calc_Patm(double, int);
double compute_pressure_gradient(double, double);	// pressure profile
double compute_air_density(double, double);	// atmospheric density
double interpolate_wind_speed(double, int);
double interpolate_wind_direction_across_360(double, int);

/* Non-Cuda Functions (start)*/
void calc_mass_loading_element(int phisize, double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *locX, double *locY, double *locZ, double *mlj, double *massreleased);
void calc_mass_loading_location(int phiint, double *mlj, double *massloading, double *ttl, double *cummassphi);
/* Non-Cuda Functions (end)*/

// CUDA function
#ifdef CUDA
void calc_mass_loading(double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *locX, double *locY, double *locZ, double *massloading_loc_source_phi, double *ttlml, double *massreleased);
void accumulate_massloading_for_phi(int phiint, double *massloading, double *ttl, double *cummassphi);
__global__ void funcD01_direct(int, int, int, int, float, float *, float *, float *, float *, float *, float *, float *, float *, float *);
__global__ void funcD01a(int, int, int, int, float, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *);
__global__ void funcD01b(int N, int LOCDIM, float *massloading_loc_source_phiD, float *ttlmlD);
#endif

//External in grain.c
double calc_particle_terminal_velocity(double h, double ashdiam, double part_density, double p, double t);
void phigenerator();
void phiconvert();
double calc_pdf_fraction(double phi);

/* Prototypes mainly used calc_mass_loading() */
/*
 * Device-side buffers used in mass loading calculation.
 *
 * This structure groups all GPU (device) memory pointers required
 * for particle transport and deposition calculations.
 *
 * Naming convention:
 *   *D suffix indicates device (GPU) memory.
 *
 * Members:
 *   massloading_loc_source_phiD :
 *       Mass loading for each (location, source, phidec)
 *
 *   ttlmlD :
 *       Total mass loading per location
 *
 *   sourceZD :
 *       Source height (z) along plume axis
 *
 *   centXD, centYD :
 *       Cloud center positions for each (s, z)
 *
 *   sigsqD :
 *       Variance (sigma^2) of particle cloud spread
 *
 *   locXD, locYD, locZD :
 *       Coordinates of observation locations
 *
 *   massreleasedD :
 *       Mass released per (s, phidec)
 */
typedef struct {
    float *massloading_loc_source_phiD;
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


/*
 * Host-side buffers corresponding to device buffers.
 *
 * This structure stores data on the CPU (host) that are:
 * - copied to the GPU
 * - or copied back from the GPU after computation
 *
 * Naming convention:
 *   *F suffix indicates host (CPU) memory.
 */
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

/*
 * Unified container for host and device buffers.
 *
 * This structure manages all memory used in the mass-loading calculation,
 * keeping host (CPU) and device (GPU) buffers together.
 *
 * Size definitions:
 *   PSZ :
 *     Total number of elements for arrays indexed by (phidec, s, z)
 *       = PHIDECDIM * SDIM * ZDIM
 *
 *   LSP :
 *     Total number of elements for arrays indexed by (location, s, phidec)
 *       = LOCDIM * SDIM * PHIDECDIM
 *
 * Members:
 *   device :
 *     Device (GPU) buffers
 *
 *   host :
 *     Host (CPU) buffers
 *
 * This structure centralizes memory management and simplifies
 * host–device data transfer and allocation/free operations.
 */
typedef struct {
    size_t PSZ;
    size_t LSP;

    DeviceBuffers device;
    HostBuffers host;
} Buffers;
//

static void copy_location_data_to_device_buffers(Buffers *b, int locN);
static void launch_mass_loading_kernels_buffers(Buffers *b, int locN);
static void copy_result_to_host_buffers(Buffers *b, double *ttlml, int loc0, int locN);
static void allocate_device_buffers_struct(Buffers *b, int chunk_locdim);
static void pack_fixed_data_buffers(Buffers *b,double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *massreleased);
static void copy_fixed_data_to_device_buffers(Buffers *b);
static void pack_location_data_buffers(Buffers *b, double *locX, double *locY, double *locZ, int loc0, int locN);
static void prepare_mass_loading(Buffers *b, double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *massreleased, int chunk_locdim);
static void compute_mass_loading(Buffers *b, double *locX, double *locY, double *locZ, double *ttlml, int loc0, int locN);
static void cleanup_buffers(Buffers *b);

void read_wind_file(const char *filename, int *windlinenum,
                    double **wind_alt, double **wind_v, double **wind_dir,
                    double **wind_tmp, double **wind_pres);

void read_loc_file(const char *filename, int *locdim,
                   double **locX, double **locY, double **locZ);

void build_plume_and_sources(int windlinenum,
                             double *wind_alt, double *wind_v, double *wind_dir,
                             double *wind_tmp, double *wind_pres,
                             double **plume_trajX, double **plume_trajY,
                             double **plume_trajZ, double **plume_trajR,
                             double **plume_trajT,
                             double **sourceX, double **sourceY,
                             double **sourceZ, double **sourceRadius,
                             double **sourceT);

void write_plume_files(int write_column_files,
                       int sdim_for_plume_calc,
                       int sdim_for_fall_calc,
                       double *plume_trajX, double *plume_trajY,
                       double *plume_trajZ, double *plume_trajR,
                       double *plume_trajT,
                       double *sourceX, double *sourceY,
                       double *sourceZ, double *sourceRadius,
                       double *sourceT);

void build_atmosphere_tables(int windlinenum, double Ht,
                             double *wind_v, double *wind_dir,
                             double *wind_tmp, double *wind_pres,
                             double **h, double **atmT, double **atmP,
                             double **windX, double **windY,
                             int *zmax);

void allocate_work_arrays(double **ttlfalltime, double **driftX, double **driftY,
                          double **ttlfalltime_phiint,
                          double **ttldriftX_phiint,
                          double **ttldriftY_phiint,
                          double **ttlfalltime_phidec,
                          double **ttldriftX_phidec,
                          double **ttldriftY_phidec,
                          double **massreleased_per_ds_and_phidec,
                          SEG **massreleased_per_ds,
                          double **cloud_center_x, double **cloud_center_y,
                          double **cloud_sigma2,
                          double **tmpmassloading,
                          double **ttlmassloading,
                          double **cummassphi,
                          DEP **location_properties,
                          RELEASE **r);

void initialize_simulation_state(DEP *location_properties,
                                 double *locX, double *locY, double *locZ,
                                 double *ttlmassloading,
                                 double *cummassphi,
                                 RELEASE *r);

void calculate_massloading(int zmax,
                           double *h, double *atmP, double *atmT,
                           double *windX, double *windY,
                           double *locX, double *locY, double *locZ,
                           double *sourceX, double *sourceY,
                           double *sourceZ, double *sourceRadius,
                           double *ttlfalltime,
                           double *driftX,
                           double *driftY,
                           double *ttlfalltime_phiint,
                           double *ttldriftX_phiint,
                           double *ttldriftY_phiint,
                           double *ttlfalltime_phidec,
                           double *ttldriftX_phidec,
                           double *ttldriftY_phidec,
                           double *massreleased_per_ds_and_phidec,
                           SEG *massreleased_per_ds,
                           double *cloud_center_x,
                           double *cloud_center_y,
                           double *cloud_sigma2,
                           double *tmpmassloading,
                           double *ttlmassloading,
                           double *cummassphi,
                           DEP *location_properties,
                           RELEASE *r);

void free_all(double *wind_alt, double *wind_v, double *wind_dir,
              double *wind_tmp, double *wind_pres,
              double *locX, double *locY, double *locZ,
              double *plume_trajX, double *plume_trajY,
              double *plume_trajZ, double *plume_trajR,
              double *plume_trajT,
              double *sourceX, double *sourceY, double *sourceZ,
              double *sourceRadius, double *sourceT,
              double *h, double *atmT, double *atmP,
              double *windX, double *windY,
              double *ttlfalltime, double *driftX, double *driftY,
              double *ttlfalltime_phiint,
              double *ttldriftX_phiint,
              double *ttldriftY_phiint,
              double *ttlfalltime_phidec,
              double *ttldriftX_phidec,
              double *ttldriftY_phidec,
              double *massreleased_per_ds_and_phidec,
              SEG *massreleased_per_ds,
              double *cloud_center_x, double *cloud_center_y,
              double *cloud_sigma2,
              double *tmpmassloading,
              double *ttlmassloading,
              double *cummassphi,
              DEP *location_properties,
              RELEASE *r);





/////////////////////////////////[PART 01]///////////////////////////////// 
//
//
// ==============================
// [1] Main
//
// simulation workflow:
// input → setup → compute → output → cleanup
// =========================

int main(int argc, char *argv[]) {
	/* 1. INPUT */

	if (argc != 4) {
		fprintf(stderr,
			"Usage: %s config.txt wind.txt loc.txt\n"
			"Example: %s conf.txt wind.dat loc.dat\n",
			argv[0], argv[0]);
		exit(EXIT_FAILURE);
	}

	/* 1.1. read config file */
	init_globals(argv[1]);		// to set WRITE_CONF which shows parameter in conf before run

	/* 1.1.1. set calculate grain size range*/
	PHIDECDIM = (int)round(1.0 / INTERVAL_DECIMAL_PHI);
	phiconvert();   // make max and minimum grain sizes ordered

	/* 1.2. read wind (atmospheric) file */
	int windlinenum = 0;
	double *wind_alt, *wind_v, *wind_dir, *wind_tmp, *wind_pres;
	read_wind_file(
		argv[2],
		&windlinenum,
		&wind_alt,
		&wind_v,
		&wind_dir,
		&wind_tmp,
		&wind_pres
	);

	/* 1.3. read loc file */
	double *locX, *locY, *locZ;
	read_loc_file(
		argv[3],
		&LOCDIM,
		&locX,
		&locY,
		&locZ
	);

	/* 2. PLUME and SOURCE SETUP */
	/* 2.1. Set plume intervals for plume calculation and source distribution */
	SDIM_FOR_PLUME_CALC = (S_MAX + S_DELTA_FOR_PLUME_CALC - 1) / S_DELTA_FOR_PLUME_CALC; // dimension of plume calculation; calculated plume should be longer than source plume
	SDIM_FOR_FALL_CALC = S_MAX / S_DELTA_FOR_FALL_CALC; // dimension of source
	SDIMCUTOFF = SDIM_FOR_FALL_CALC;

	/* 2.2. Obtain plume trajectory and set particle sources on it */
	double *plume_trajX, *plume_trajY, *plume_trajZ, *plume_trajR, *plume_trajT;
	double *sourceX, *sourceY, *sourceZ, *sourceRadius, *sourceT;
	build_plume_and_sources(
		windlinenum,
		wind_alt, wind_v, wind_dir, wind_tmp, wind_pres,
		&plume_trajX, &plume_trajY, &plume_trajZ, &plume_trajR, &plume_trajT,
		&sourceX, &sourceY, &sourceZ, &sourceRadius, &sourceT
	);
	
	/* 2.3. Write original trajectory file */
	write_plume_files(
		WRITE_COLUMN_FILES,
		SDIM_FOR_PLUME_CALC,
		SDIM_FOR_FALL_CALC,
		plume_trajX, plume_trajY, plume_trajZ, plume_trajR, plume_trajT,
		sourceX, sourceY, sourceZ, sourceRadius, sourceT
	);

	/* 3. Atmosphere and wind field (table preparation) */
	double *h, *atmT, *atmP, *windX, *windY;
	int zmax;
	build_atmosphere_tables(
		windlinenum,
		Ht,
		wind_v,
		wind_dir,
		wind_tmp,
		wind_pres,
		&h,
		&atmT,
		&atmP,
		&windX,
		&windY,
		&zmax
	);

	/*
	* 4. ALLOCATE WORKING ARRAYS FOR THE SIMULATION
	*
	* These arrays are used across multiple stages of the computation:
	* - fall time and drift calculations
	* - per-grainsize summaries
	* - per-source (s) distributions
	* - mass loading at ground locations
	*
	* All arrays are allocated here to centralize memory management
	* and make the data dependencies of the main computation explicit.
	*
	* Note:
	* The sizes depend on global dimensions such as ZDIM, PHIDECDIM,
	* SDIM_FOR_FALL_CALC, and LOCDIM.
	*/

	// TODO: These arrays can be grouped into a SimulationWorkspace struct

	// temporary arrays for fall time and drift during particle descent
	// used inside the phi loop in calculate_massloading()
	//
	// ttlfalltime[z] : elapsed time for a particle to reach height interval z [s]
	// driftX[z]      : horizontal drift distance in X-direction up to height interval z [m]
	// driftY[z]      : horizontal drift distance in Y-direction up to height interval z [m]
	//
	// dimensions: ZDIM
	double *ttlfalltime, *driftX, *driftY;

	// summary of fall time and drift for each integer phi class
	// ttlfalltime_phiint[phi][z] : elapsed time for particles of integer phi to reach height z [s]
	// ttldriftX_phiint[phi][z]   : horizontal drift distance in X-direction for integer phi [m]
	// ttldriftY_phiint[phi][z]   : horizontal drift distance in Y-direction for integer phi [m]
	// dimensions: (number of integer phi classes) × ZDIM
	double *ttlfalltime_phiint, *ttldriftX_phiint, *ttldriftY_phiint;

	// fall time and drift for each decimal phi (phidec) class
	// ttlfalltime_phidec[phidec][z] : elapsed time to reach height z [s]
	// ttldriftX_phidec[phidec][z]   : horizontal drift distance in X-direction [m]
	// ttldriftY_phidec[phidec][z]   : horizontal drift distance in Y-direction [m]
	// dimensions: PHIDECDIM × ZDIM
	double *ttlfalltime_phidec, *ttldriftX_phidec, *ttldriftY_phidec;

	
	double *massreleased_per_ds_and_phidec;
	SEG *massreleased_per_ds;

	// cloud center position and diffusion
	// (cloud: group of particles with same grain size [phi + phidec; e.g. fraction of 3.5 phi] 
	// released from a certain source [s; e.g. 53rd interval from the vent])
	double *cloud_center_x, *cloud_center_y, *cloud_sigma2;


	double *tmpmassloading, *ttlmassloading, *cummassphi;
	DEP *location_properties;
	RELEASE *r;

	allocate_work_arrays(
		&ttlfalltime, &driftX, &driftY,
		&ttlfalltime_phiint, &ttldriftX_phiint, &ttldriftY_phiint,
		&ttlfalltime_phidec, &ttldriftX_phidec, &ttldriftY_phidec,
		&massreleased_per_ds_and_phidec,
		&massreleased_per_ds,
		&cloud_center_x, &cloud_center_y, &cloud_sigma2,
		&tmpmassloading, &ttlmassloading, &cummassphi,
		&location_properties, &r
	);

	/* 5. INITIALIZATION*/
	initialize_simulation_state(
		location_properties,
		locX, locY, locZ,
		ttlmassloading,
		cummassphi,
		r
	);

	/* 6. MAIN COMPUTATION */
	calculate_massloading(
		zmax,
		h, atmP, atmT, windX, windY,
		locX, locY, locZ,
		sourceX, sourceY, sourceZ, sourceRadius,
		ttlfalltime,
		driftX,
		driftY,
		ttlfalltime_phiint,
		ttldriftX_phiint,
		ttldriftY_phiint,
		ttlfalltime_phidec,
		ttldriftX_phidec,
		ttldriftY_phidec,
		massreleased_per_ds_and_phidec,
		massreleased_per_ds,
		cloud_center_x,
		cloud_center_y,
		cloud_sigma2,
		tmpmassloading,
		ttlmassloading,
		cummassphi,
		location_properties,
		r
	);


	/* 7.1. massloading and isopach */
	store_total_massloading_and_mean_phi(location_properties, ttlmassloading, cummassphi);	//calculate mean diameter for each location
	if(WRITE_MASSLOADING){write_deposit_summary(location_properties);		//massloading.txt
	createisopachdata(location_properties);}						//S_vs_Area.txt


	/* 7. OUTPUT */
	/* 7.2. Particle relase and diffusion */
	if(WRITE_FALL_INFO_FILES){	// DEFINED IN CONFIG FILE
	const char *name1 = "falltime.txt";
	write_vertical_profiles_for_phi(name1, h, ttlfalltime_phiint, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	const char *name2 = "falldriftX.txt";
	write_vertical_profiles_for_phi(name2, h, ttldriftX_phiint, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	const char *name3 = "falldriftY.txt";
	write_vertical_profiles_for_phi(name3, h, ttldriftY_phiint, MAX_GRAINSIZE + 1, MIN_GRAINSIZE - MAX_GRAINSIZE, -1);
	write_particle_release_theoretical_vs_actual(r);			//particle_released.txt
	write_massrelease_per_source_phi(massreleased_per_ds); //segregation_per_ds.txt
	}
	
	/* 8. CLEAN UP */
	free_all(
		wind_alt, wind_v, wind_dir, wind_tmp, wind_pres,
		locX, locY, locZ,
		plume_trajX, plume_trajY, plume_trajZ, plume_trajR, plume_trajT,
		sourceX, sourceY, sourceZ, sourceRadius, sourceT,
		h, atmT, atmP, windX, windY,
		ttlfalltime, driftX, driftY,
		ttlfalltime_phiint, ttldriftX_phiint, ttldriftY_phiint,
		ttlfalltime_phidec, ttldriftX_phidec, ttldriftY_phidec,
		massreleased_per_ds_and_phidec, massreleased_per_ds,
		cloud_center_x, cloud_center_y, cloud_sigma2,
		tmpmassloading, ttlmassloading, cummassphi,
		location_properties, r
	);

	return 0;

}	// End of main 

/////////////////////////////////[PART 01]///////////////////////////////// 
//
//
// ==============================
// [1] Input
// ==============================
// [1.1.] read_wind_file
// [1.1.1.] read_wind
// [1.2] read_loc_file
// [1.2.1.] read_locations_and_convert_to_vent_centered
// [1.3] phiconvert
// 
/*[1.1.]
 * Read raw atmospheric (wind) data from the input file.
 *
 * This function only reads the original input profiles and stores them
 * as given in the file. Interpolation onto the simulation height grid
 * is done later by build_atmosphere_tables().
 *
 * Outputs:
 * - wind_alt[i]  : altitude of input data point [m]
 * - wind_v[i]    : wind speed [m/s]
 * - wind_dir[i]  : wind direction [deg]
 * - wind_tmp[i]  : temperature [K]
 * - wind_pres[i] : pressure [Pa]
 */
void read_wind_file(
    const char *filename,
    int *windlinenum,
    double **wind_alt,
    double **wind_v,
    double **wind_dir,
    double **wind_tmp,
    double **wind_pres
){
    FILE *in_wind = fopen(filename, "r");
	if (!in_wind) {
		fprintf(stderr, "Error: cannot open wind file: %s\n", filename);
		perror("fopen");
		exit(EXIT_FAILURE);
	}

    *windlinenum = get_wind_line_number(in_wind);
    rewind(in_wind);

	// allocate arrays to store input atmospheric (wind) data (argv[2])
	// wind_alt[i]  : altitude [m]
	// wind_v[i]    : wind speed [m/s]
	// wind_dir[i]  : wind direction [deg]
	// wind_tmp[i]  : temperature [K]
	// wind_pres[i] : pressure [Pa]
	// dimension: windlinenum (number of input data points)
    *wind_alt  = (double *)malloc(*windlinenum * sizeof(double));
    *wind_v    = (double *)malloc(*windlinenum * sizeof(double));
    *wind_dir  = (double *)malloc(*windlinenum * sizeof(double));
    *wind_tmp  = (double *)malloc(*windlinenum * sizeof(double));
    *wind_pres = (double *)malloc(*windlinenum * sizeof(double));

    read_wind_profile(in_wind, *wind_alt, *wind_v, *wind_dir, *wind_tmp, *wind_pres);

    fclose(in_wind);
}

/*
 * [1.1.1] Read wind and atmospheric profiles from input file.
 *
 * The input file must contain vertical profiles of:
 *
 *   height [m]
 *   wind speed
 *   wind direction
 *   atmospheric temperature
 *   atmospheric pressure
 *
 * Comment lines beginning with '#' are ignored.
 *
 * If the first valid level is located above z = 0 m,
 * an additional synthetic surface-level entry is inserted.
 *
 * The inserted surface values use:
 *   - zero wind speed
 *   - the same wind direction
 *   - simple lapse-rate extrapolation for temperature
 *   - simple hydrostatic extrapolation for pressure
 *
 * This guarantees that the atmospheric profile
 * always begins at ground level.
 *
 * Note:
 * The parsing logic must remain consistent with
 * get_wind_line_number().
 */
void read_wind_profile(FILE *f, double *height, double *speed, double *dir, double *atm_temp, double *atm_pres){
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

/* [1.2.]
 * Read ground location data from the input file.
 *
 * This function only reads the original location coordinates.
 * These coordinates are later stored in location_properties by
 * initialize_simulation_state().
 *
 * Outputs:
 * - locX[j] : X-coordinate of ground location j [m]
 * - locY[j] : Y-coordinate of ground location j [m]
 * - locZ[j] : elevation of ground location j [m]
 */
void read_loc_file(
    const char *filename,
    int *locdim,
    double **locX,
    double **locY,
    double **locZ
){
    FILE *in_loc = fopen(filename, "r");
	if (!in_loc) {
		fprintf(stderr, "Error: cannot open location file: %s\n", filename);
		perror("fopen");
		exit(EXIT_FAILURE);
	}

    *locdim = get_line_number(in_loc);
    rewind(in_loc);

    *locX = (double *)malloc(*locdim * sizeof(double));
    *locY = (double *)malloc(*locdim * sizeof(double));
    *locZ = (double *)malloc(*locdim * sizeof(double));

    read_locations_and_convert_to_vent_centered(in_loc, *locX, *locY, *locZ);

    fclose(in_loc);
}

/* [1.2.1.]
 * Read location coordinates and convert them to a vent-centered coordinate system.
 *
 * Input coordinates are given in map coordinates. This function converts them
 * to a vent-centered system (vent = 0):
 *   x = X_map - VENT_EASTING
 *   y = Y_map - VENT_NORTHING
 *
 * For each location:
 * - coordinates are stored relative to the vent
 * - negative elevation is clipped to z = 0
 * - map boundaries (MAPENDW/E/S/N) are updated
 *
 * Note:
 * The entire simulation uses a vent-centered coordinate system.
 */
void read_locations_and_convert_to_vent_centered(FILE *f, double *x, double *y, double *z){
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

/* [1.3.]
 *	The phi scale is inverted (smaller values mean larger particles),
 *	so any incorrect ordering is automatically corrected by this function
 */
void phiconvert(){
	if(MIN_GRAINSIZE < MAX_GRAINSIZE){
		int tmp = MIN_GRAINSIZE;
		MIN_GRAINSIZE = MAX_GRAINSIZE;
		MAX_GRAINSIZE = tmp;
	}
}
/////////////////////[END OF THE PART 01]/////////////////////

/////////////////////////////////[PART 02]///////////////////////////////// 
//
//
// ==============================
// [2] Initialization
// ==============================
// [2.1] init_globals
// [2.2] allocate_work_arrays
// [2.3] initialize_simulation_state
// 
/* [2.1.]
 * Initialize global parameters from configuration file.
 *
 * This function reads key-value pairs from the specified configuration file
 * and assigns values to global variables used throughout the simulation.
 *
 * Supported parameters include:
 * - physical constants (e.g., DIFFUSION_COEFFICIENT, EDDY_CONST)
 * - eruption conditions (e.g., ERUPTION_MASS, MAGMA_DISCHARGE_RATE)
 * - grain size distribution parameters (e.g., MAX_GRAINSIZE, STD_GRAINSIZE)
 * - numerical settings (e.g., Z_DELTA, S_DELTA_FOR_FALL_CALC)
 * - output control flags (e.g., WRITE_MASSLOADING, WRITE_COLUMN_FILES)
 *
 * Notes:
 * - Lines starting with '#' or empty lines are ignored.
 * - Values are parsed as either double or integer depending on parameter.
 * - Pressure input is assumed to be in hPa where applicable.
 * - If WRITE_CONF is enabled, all parsed values are printed to stderr.
 *
 * Return:
 *   0 on success, non-zero if file cannot be opened.
 */
int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;

  in_config = fopen(config_file, "r");
  	if (!config_file) {
		fprintf(stderr, "Error: cannot open config file: %s\n", config_file);
		perror("fopen");
		exit(EXIT_FAILURE);
	}

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
	  	if (Z_DELTA <= 0) {
			fprintf(stderr, "Invalid Z_DELTA: %s\n", token);
			exit(EXIT_FAILURE);
		}
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
	  if (S_DELTA_FOR_PLUME_CALC <= 0) {
			fprintf(stderr, "Invalid S_DELTA_FOR_PLUME_CALC: %s\n", token);
			exit(EXIT_FAILURE);
		}
      if(WRITE_CONF) fprintf(stderr, "S_DELTA_FOR_PLUME_CALC = %1.1f\n", S_DELTA_FOR_PLUME_CALC);
    }
    else if (!strncmp(token, "S_DELTA_FOR_FALL_CALC", strlen("S_DELTA_FOR_FALL_CALC"))) {
      token = strtok_r(NULL, space, ptr1);
      S_DELTA_FOR_FALL_CALC = (int)atoi(token);
		if (S_DELTA_FOR_FALL_CALC <= 0) {
			fprintf(stderr, "Invalid S_DELTA_FOR_FALL_CALC: %s\n", token);
			exit(EXIT_FAILURE);
		}
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
    else if (!strncmp(token, "PRINT_PROGRESS", strlen("PRINT_PROGRESS"))) {
      token = strtok_r(NULL, space, ptr1);
      PRINT_PROGRESS = strtod(token, NULL);
      if(WRITE_CONF) fprintf(stderr, "PRINT_PROGRESS = %d\n", PRINT_PROGRESS);
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


/* [2.2.]
 * Allocate working arrays used in the main mass-loading calculation.
 *
 * These arrays store temporary fall/drift profiles, per-phi summaries,
 * particle release distributions, cloud-center positions, cloud dispersion,
 * and accumulated deposit information at ground locations.
 *
 * Main array dimensions:
 * - ZDIM                                  : vertical grid
 * - PHIDECDIM                             : decimal phi classes
 * - SDIM_FOR_FALL_CALC                    : source points along plume axis
 * - LOCDIM                                : ground locations
 * - MIN_GRAINSIZE - MAX_GRAINSIZE         : integer phi classes
 *
 * Note:
 * Arrays are allocated here only. Their physical meanings are documented
 * at their declarations in main() and in the functions where they are used.
 */
void allocate_work_arrays(
    double **ttlfalltime,
    double **driftX,
    double **driftY,
    double **ttlfalltime_phiint,
    double **ttldriftX_phiint,
    double **ttldriftY_phiint,
    double **ttlfalltime_phidec,
    double **ttldriftX_phidec,
    double **ttldriftY_phidec,
    double **massreleased_per_ds_and_phidec,
    SEG **massreleased_per_ds,
    double **cloud_center_x,
    double **cloud_center_y,
    double **cloud_sigma2,
    double **tmpmassloading,
    double **ttlmassloading,
    double **cummassphi,
    DEP **location_properties,
    RELEASE **r
){
	// temporary fall-time and drift profiles
		// particle data, temporaly storage
    *ttlfalltime = (double*)calloc(ZDIM, sizeof(double));
    *driftX      = (double*)calloc(ZDIM, sizeof(double));
    *driftY      = (double*)calloc(ZDIM, sizeof(double));
		// particle data during fall for integer phi classes
    *ttlfalltime_phiint = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
    *ttldriftX_phiint   = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
    *ttldriftY_phiint   = (double*)calloc(ZDIM * (MIN_GRAINSIZE - MAX_GRAINSIZE), sizeof(double));
		// particle data during fall for integer decunak phi classes
    *ttlfalltime_phidec = (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));
    *ttldriftX_phidec   = (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));
    *ttldriftY_phidec   = (double*)calloc(ZDIM * PHIDECDIM, sizeof(double));

	// particle release distributions along plume axis
    *massreleased_per_ds_and_phidec = (double*)calloc(SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
    *massreleased_per_ds = (SEG*)calloc(SDIM_FOR_FALL_CALC, sizeof(SEG));

	// cloud center positions and horizontal dispersion
    *cloud_center_x = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
    *cloud_center_y = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));
    *cloud_sigma2 = (double*)calloc(ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM, sizeof(double));

	// mass loading and deposit accumulation at ground locations
    *tmpmassloading = (double*)calloc(LOCDIM, sizeof(double));
    *ttlmassloading = (double*)calloc(LOCDIM, sizeof(double));
    *cummassphi     = (double*)calloc(LOCDIM, sizeof(double));

	// output data structures
    *location_properties = (DEP *)calloc(LOCDIM, sizeof(DEP));
    *r = (RELEASE *)calloc(MIN_GRAINSIZE - MAX_GRAINSIZE, sizeof(RELEASE));
}

/* [2.3.]
 * Initialize simulation state before main computation.
 *
 * This function:
 * - assigns ground location coordinates to location_properties
 * - resets mass loading and cumulative mass arrays
 * - initializes theoretical released mass for each phi class
 */
void initialize_simulation_state(
    DEP *location_properties,
    double *locX,
    double *locY,
    double *locZ,
    double *ttlmassloading,
    double *cummassphi,
    RELEASE *r
){
	/* location_properties is a structure storing final deposition results for each location (massloading.txt) */
    set_coordinates_to_location_properties(location_properties, locX, locY, locZ);

    clear_array(LOCDIM, ttlmassloading);
    clear_array(LOCDIM, cummassphi);

    compute_theoretical_particle_release(r);
}

/* [2.3.1]
 * Compute theoretical particle release for each integer phi class.
 *
 * The release mass is obtained by integrating the grain-size PDF
 * over decimal phi bins within each integer phi interval, then
 * multiplying by the total eruption mass.
 */
void compute_theoretical_particle_release(RELEASE *r){
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
/////////////////////[END OF THE PART 02]/////////////////////


/////////////////////////////////[PART 03]///////////////////////////////// 
//
// ==============================
// [3] Plume & Source Construction
// ==============================
// [3.1] build_plume_and_sources
//   [3.1.1] plume_calculation
//       [3.1.1.1] advance_plume_state_rk4
//       [3.1.1.2] calc_plume_heat_capacity
//       [3.1.1.3] calc_Cp0
//       [3.1.1.4] func12
//       [3.1.1.5] func13
//       [3.1.1.6] func14
//       [3.1.1.7] func15
//       [3.1.1.8] func16
//       [3.1.1.9] func17
//       [3.1.1.10] func18
//       [3.1.1.11] func19
//   [3.1.2] set_source_points_on_plume
//   [3.1.3] makewindstruct

/* [3.1.]
 * Compute plume trajectory and define particle source points along it.
 *
 * First, solve the plume trajectory using input atmospheric data and store
 * the plume centerline coordinates, radius, and travel time.
 * Then, place discrete particle source points along the plume axis for
 * the fall and mass-loading calculations.
 *
 * Outputs:
 * - plume_trajX/Y/Z : plume centerline coordinates [m]
 * - plume_trajR     : plume radius [m]
 * - plume_trajT     : elapsed time from vent along plume axis [s]
 * - sourceX/Y/Z     : particle source coordinates [m]
 * - sourceRadius    : plume radius at each source point [m]
 * - sourceT         : elapsed time from vent to each source point [s]
 */
void build_plume_and_sources(
    int windlinenum,
    double *wind_alt,
    double *wind_v,
    double *wind_dir,
    double *wind_tmp,
    double *wind_pres,
    double **plume_trajX,
    double **plume_trajY,
    double **plume_trajZ,
    double **plume_trajR,
    double **plume_trajT,
    double **sourceX,
    double **sourceY,
    double **sourceZ,
    double **sourceRadius,
    double **sourceT
){
    // allocate arrays for plume trajectory
    *plume_trajX = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
    *plume_trajY = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
    *plume_trajZ = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
    *plume_trajR = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));
    *plume_trajT = (double *)malloc(SDIM_FOR_PLUME_CALC * sizeof(double));

	// compute plume centerline
    Ht = plume_calculation(
        windlinenum,
        *plume_trajX, *plume_trajY, *plume_trajZ,
        *plume_trajR, *plume_trajT,
        wind_alt, wind_v, wind_dir, wind_tmp, wind_pres
    );

    // allocate arrays for particle source points
    *sourceX = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
    *sourceY = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
    *sourceZ = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
    *sourceRadius = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));
    *sourceT = (double *)malloc(SDIM_FOR_FALL_CALC * sizeof(double));

	// interpolate source points along plume trajectory
    set_source_points_on_plume(
        *sourceX, *sourceY, *sourceZ, *sourceRadius, *sourceT,
        *plume_trajX, *plume_trajY, *plume_trajZ, *plume_trajR, *plume_trajT
    );
}

/* [3.1.1.]
 * Compute plume trajectory and source-point properties.
 *
 * This is the main plume calculation routine.
 * Starting from vent conditions, the plume state is advanced step by step
 * along the plume axis using advance_plume_state_rk4().
 *
 * At each step:
 * - plume state variables (Q, M, theta, E, etc.) are updated
 * - atmospheric conditions are interpolated from input data
 * - plume position and properties are stored as source points
 *
 * The integration continues until plume rise stops (theta <= 0 or M <= 0),
 * and the plume top height Ht is determined.
 *
 * Outputs:
 * - plume centerline trajectory
 * - source-point positions and properties for particle calculations
 * - plume height Ht
 */
double plume_calculation(int imax, double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *taftervent, double *wind_alt, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres){	// The main routine in this file
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
	Cp0 = calc_Cp0();

	rho_a = compute_air_density(p, ta);
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
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2012)
			U = Q / (rho_c * VENT_RADIUS * VENT_RADIUS);
			R = VENT_RADIUS;
		}else if(MAGMA_DISCHARGE_RATE > 0 && INITIAL_PLUME_VELOCITY > 0 && VENT_RADIUS < 0){
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2012)
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

	V = interpolate_wind_speed(gz, total);		// wind velocity
	g_dir = interpolate_wind_direction_across_360(gz, total);	// wind direction
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
		
		advance_plume_state_rk4(total, T);
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
	V = interpolate_wind_speed(Ht, total);
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


/* [3.1.1.1]
 * Advance plume state by one step using 4th-order Runge-Kutta (RK4).
 *
 * This function integrates the governing plume equations along the
 * trajectory coordinate s, updating:
 *   Q      : mass flux
 *   M      : momentum flux
 *   theta  : plume angle
 *   E      : energy flux
 *
 * At each RK stage, atmospheric conditions (pressure, temperature,
 * wind speed/direction) are interpolated based on height.
 *
 * The plume properties (density, velocity, radius, etc.) are updated
 * consistently with the current state.
 */
void advance_plume_state_rk4(int total, double T){
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
	//dp_over_ds1 = compute_pressure_gradient(p, ta);
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
	rho_a_tmp = compute_air_density(p_tmp, ta);

	Cp = calc_plume_heat_capacity(n);					// calc Cp
	//printf("k1\n");
	V  = interpolate_wind_speed(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k2     ///////////////////////
	//dp_over_ds2 = compute_pressure_gradient(p_tmp, ta);
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
	rho_a_tmp = compute_air_density(p_tmp, ta);

	Cp = calc_plume_heat_capacity(n);					// calc Cp
	//printf("k2\n");
	V = interpolate_wind_speed(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k3     ///////////////////////
	//dp_over_ds3 = compute_pressure_gradient(p_tmp, ta);
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
	rho_a_tmp = compute_air_density(p_tmp, ta);

	Cp = calc_plume_heat_capacity(n);					// calc Cp
	//printf("k3\n");
	V = interpolate_wind_speed(gz+ dz, total);

	T = E_tmp / Q_tmp / Cp;

	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));

	/////////////////////////////////
	// k4     ///////////////////////
	//dp_over_ds4 = compute_pressure_gradient(p_tmp, ta);
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
	rho_a = compute_air_density(p, ta);

	Cp = calc_plume_heat_capacity(n);					// calc Cp
	//printf("k4\n");
	V = interpolate_wind_speed(gz, total);
	g_dir = interpolate_wind_direction_across_360(gz, total);



	T = E / Q / Cp;

	Rg = func19(n);					// calc Rg
	rho_c = func17(n, p, Rg, T);	// calc plume density (rho_c)

	Ue = func16(M, Q, theta_tmp, V);
	U = M / Q;
	R = sqrt(Q / (U * rho_c));
}

/* [3.1.1.2.]
 * Compute mixture heat capacity Cp for a given gas fraction n.
 * (Equation 20 in Woodhouse et al.)
 *
 * Cp is linearly interpolated between:
 * - Ca : heat capacity of air
 * - Cp0: initial mixture heat capacity
 *
 * using gas fraction n:
 *
 *   Cp = Ca + (Cp0 - Ca) * (1 - n) / (1 - n0)
 */
double calc_plume_heat_capacity(double n_tmp){
	double Cp_tmp;
	Cp_tmp = Ca + (Cp0 - Ca) * (1 - n_tmp) / (1 - n0);

	return Cp_tmp;
}

/* [3.1.1.3.]
 * Compute mixture initial specific heat capacity of plume (Cp0).
 *
 * n0 = initial gas fraction
 * Cv = specific heat capacity of water vapor
 * Cs = specific heat capacity of solid pyroclast
 */
double calc_Cp0(){
	return n0 * Cv + (1 - n0) * Cs;
}


/* */
/* Function numbers correspond that in Woodhouse et al (2012) */
/* [3.1.1.4] */
/*  Flux change along the plume axis*/
double func12(double M_tmp, double rho_a_tmp, double rho_c_tmp, double Q_tmp){			// plume mass flux
	double dQ_over_ds;

	dQ_over_ds = 2 * rho_a_tmp * Ue * Q_tmp / sqrt(rho_c_tmp * M_tmp);

	//printf("%1.4f\n", dQ_over_ds);
	return dQ_over_ds;
}

/* [3.1.1.5] */
/*  Momentum flux variation along the plume axis　*/
double func13(double rho_a_tmp, double rho_c_tmp, double M_tmp, double theta_tmp, double Q_tmp){
	double dM_over_ds;

	dM_over_ds = GRAVITY * (rho_a_tmp - rho_c) * Q_tmp * Q_tmp / (rho_c_tmp * M_tmp) * sin(theta_tmp);
	dM_over_ds = dM_over_ds + 2 * rho_a_tmp * Q_tmp / sqrt(rho_c_tmp * M_tmp) * Ue * V * cos(theta_tmp);

	return dM_over_ds;
}

/* [3.1.1.6] */
/*  Plume axis angle variation along the plume axis　*/
double func14(double M_tmp, double Q_tmp, double rho_a_tmp, double rho_c_tmp, double theta_tmp){
	double dtheta_over_ds;

	dtheta_over_ds = GRAVITY * (rho_a_tmp - rho_c_tmp) * Q_tmp * Q_tmp * cos(theta_tmp) / (rho_c_tmp * M_tmp * M_tmp);
	dtheta_over_ds = dtheta_over_ds - 2 * rho_a_tmp * Q_tmp * Ue * V * sin(theta_tmp) / (M * sqrt(rho_c_tmp * M_tmp));

	//printf("dtheta = %1.4f\n", dtheta_over_ds);
	return dtheta_over_ds;
}

/* [3.1.1.7] */
/*  Enthalpy flux variation along the plume axis　*/
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

/* [3.1.1.8] */
/*  Entrainment velocity along the plume boundary　*/
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

/* [3.1.1.9] */
/*  Bulk plume density　*/
double func17(double n_tmp, double p_tmp, double Rg_tmp, double T_tmp){			// plume density
	double rho_tmp;

	rho_tmp = (1 - n_tmp) / rho_s + n_tmp * Rg_tmp * T_tmp / p_tmp;
	rho_c = 1 / rho_tmp;
	//printf("rho     = %1.4f\n", rho);
	return rho_c;
}

/* [3.1.1.10] */
/*  Gas fraction of plume　*/
double func18(double Q_tmp){			// solid content in the plume
	double n_tmp;
	n_tmp = 1 - (1 - n0) * Q0 / Q_tmp;

	return n_tmp;
}

/* [3.1.1.11] */
/*  Bulc gas constant　*/
double func19(double n_tmp){
	double Rg_tmp;
	Rg_tmp = Ra + (Rg0 - Ra) * n0 * (1 - n_tmp) / (n_tmp * (1 - n0));

	return Rg_tmp;
}

/* End of Woodhouse functions*/


/* [3.1.2.]
 * Define particle source points along the plume trajectory.
 *
 * Source points are placed at regular intervals along plume-axis distance
 * using S_DELTA_FOR_FALL_CALC. Their position, radius, and travel time are
 * obtained by linear interpolation from the plume trajectory calculated at
 * intervals of S_DELTA_FOR_PLUME_CALC.
 *
 * Outputs:
 * - sourceX/Y/Z : coordinates of particle source points [m]
 * - sourceR     : plume radius at each source point [m]
 * - sourceT     : elapsed time from vent to each source point [s]
 */
void set_source_points_on_plume(double *sourceX, double *sourceY, double *sourceZ, double *sourceR, double *sourceT, double *plume_trajX, double *plume_trajY, double *plume_trajZ, double *plume_trajR, double *plume_trajT){
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

/* [3.1.3.]
 * Build WIND structure array from input atmospheric data arrays.
 *
 * This function converts separate arrays (altitude, wind speed,
 * direction, temperature, pressure) into an array of WIND structures.
 */
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

/////////////////////[END OF THE PART 03]/////////////////////

/////////////////////////////////[PART 04]///////////////////////////////// 
//
//
// ==============================
// [4] Atmosphere & Wind
// ==============================
// [4.1] build_atmosphere_tables
//   [4.1.1] interpolate_atmosphere_and_wind
//       [4.1.1.1] calc_Tatm
//       [4.1.1.2] calc_Patm
//       [4.1.1.3] compute_air_density
//       [4.1.1.4] compute_pressure_gradient
//       [4.1.1.5] interpolate_wind_speed
//       [4.1.1.6] interpolate_wind_direction_across_360


/* [4.1.]
 *
 * Build atmospheric and wind profiles on the simulation vertical grid.
 *
 * This function interpolates raw atmospheric data (read by read_wind_file)
 * onto a regular vertical grid used in the simulation.
 *
 * The grid has uniform spacing Z_DELTA except for the top level,
 * where h[zmax] is explicitly set to the plume height Ht.
 *
 * Inputs:
 * - wind_v, wind_dir, wind_tmp, wind_pres : raw atmospheric data
 *
 * Outputs:
 * - h[z]     : height level [m]
 * - atmT[z]  : atmospheric temperature [K]
 * - atmP[z]  : atmospheric pressure [Pa]
 * - windX[z] : wind velocity in X-direction [m/s]
 * - windY[z] : wind velocity in Y-direction [m/s]
 */
void build_atmosphere_tables(
    int windlinenum,
    double Ht,
    double *wind_v,
    double *wind_dir,
    double *wind_tmp,
    double *wind_pres,
    double **h,
    double **atmT,
    double **atmP,
    double **windX,
    double **windY,
    int *zmax
){
	// define vertical grid
	// z = 0      : sea level
	// z = zmax-1 : highest regular grid level below Ht
	// z = zmax   : plume height Ht
	//
	// h[0 ... zmax-1] are regular grid points with spacing Z_DELTA.
	// h[zmax] is explicitly set to Ht.
    *zmax = ceil(Ht / Z_DELTA);
    ZDIM = *zmax + 1;

	// allocate arrays for interpolated atmospheric and wind profiles
	// defined on a uniform vertical grid with spacing Z_DELTA
	//
	// h[z]     : height grid [m]
	// atmT[z]  : temperature profile [K]
	// atmP[z]  : pressure profile [Pa]
	// windX[z] : wind velocity in X-direction [m/s]
	// windY[z] : wind velocity in Y-direction [m/s]
	//
	// dimension: ZDIM (number of vertical grid points)
    *h     = (double *)malloc(ZDIM * sizeof(double));
    *atmT  = (double *)malloc(ZDIM * sizeof(double));
    *atmP  = (double *)malloc(ZDIM * sizeof(double));
    *windX = (double *)malloc(ZDIM * sizeof(double));
    *windY = (double *)malloc(ZDIM * sizeof(double));

    for(int z = 0; z < *zmax; z++){	// uniform vertical grid with spacing Z_DELTA
        (*h)[z] = z * Z_DELTA;
    }
    (*h)[*zmax] = Ht;	// set the top grid level to the exact plume height Ht

	// interpolate input atmospheric data onto the vertical grid
    interpolate_atmosphere_and_wind(
        windlinenum,
        *h, *atmT, *atmP, *windX, *windY,		// interpolated data
        wind_v, wind_dir, wind_tmp, wind_pres	// input data (argv[2])
    );
}


/* [4.1.1.]
 * Interpolate atmospheric and wind data onto the simulation height grid.
 *
 * For each height level h[z], this function computes:
 * - atmospheric temperature atmT[z]
 * - atmospheric pressure atmP[z]
 * - wind velocity components windX[z], windY[z]
 *
 * Input wind data are read from the atmospheric file and stored in
 * wind_v, wind_dir, wind_tmp, and wind_pres.
 *
 * Notes:
 * - Input pressure is converted from hPa to Pa.
 * - Wind direction and speed are converted to X/Y components.
 * - If WRITE_COLUMN_FILES is enabled, the interpolated table is written
 *   to atmosphere_used.txt.
 */
void interpolate_atmosphere_and_wind(int windlinenum, double *h, double *atmT, double *atmP, double *windX, double *windY, double *wind_v, double *wind_dir, double *wind_tmp, double *wind_pres){
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
		v = interpolate_wind_speed(h[z], windlinenum);
		dir = interpolate_wind_direction_across_360(h[z], windlinenum);
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

/* [4.1.1.1.]
 * Interpolate atmospheric temperature at height h
 * from discrete atmospheric data.
 *
 * Linear interpolation is used between adjacent height levels.
 * If h is above the highest level, the temperature at the highest
 * level is returned.
 */
double calc_Tatm(double h, int total){
	int i=1;
	double t_atm = W1[total-1].t_atm;

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

/* [4.1.1.2]
 * Estimate atmospheric pressure at height h from discrete atmospheric data.
 *
 * Temperature is linearly interpolated between input height levels.
 * Pressure is then estimated from the pressure at the lower level
 * using a hydrostatic approximation with the mean temperature
 * between the lower level and height h.
 *
 * Input pressure is assumed to be in hPa and converted to Pa.
 */
double calc_Patm(double h, int total){
	int i=1;
	double p_atm;
	double a;

	double t_atm = W1[total-1].t_atm;
	double t_atm0 = W1[total-1].t_atm;
	double p_atm0 = W1[total-1].p_atm;

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

/* [4.1.1.3.]
 * Compute air density using the ideal gas law.
 * (Equation 23 in Woodhouse et al.)
 *
 *   rho = p / (R * T)
 */
double compute_air_density(double p_tmp, double t_tmp){	// atmospheric density
	double rho_tmp;

	rho_tmp = p_tmp / (Ra * t_tmp);
	return rho_tmp;
}

/* [4.1.1.4.]
 * Compute vertical pressure gradient under hydrostatic balance.
 * (Equation 22 in Woodhouse et al.)
 *
 *   dp/ds = - (g * p) / (R * T)
 */
double compute_pressure_gradient(double p_tmp, double t_tmp){	// atmospheric pressure
	double dp_over_ds;

	dp_over_ds = -1 * (GRAVITY * p_tmp) / (Ra * t_tmp);
	return dp_over_ds;
}

/* [4.1.1.5.]
 * Interpolate wind speed at height h from discrete wind data.
 *
 * If h is above the highest wind-data level, the wind speed at the
 * highest level is returned.
 */
double interpolate_wind_speed(double h, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double v;

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

/* [4.1.1.6.]
 * Interpolate wind direction at height h from discrete wind data.
 *
 * Wind direction is circular data, so this function handles wrap-around
 * across 0/360 degrees before linear interpolation.
 *
 * If h is above the highest wind-data level, the wind direction at the
 * highest level is returned.
 */
double interpolate_wind_direction_across_360(double h, int total){	//return wind direction based on discrete wind data
	int i=1;
	double dir;
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
} // end of 4.1.1.6.

/////////////////////[END OF THE PART 04]/////////////////////

/////////////////////////////////[PART 05]///////////////////////////////// 
//
//
/// ==============================
// [5] Mass Loading Core
// ==============================
// [5.1] calculate_massloading

// calculate_massloading helper
// [5.1.1.] write_phi_s_table
// [5.1.2.] store_profile_for_phi

// --- Release from plume ---
// [5.2.1.]  compute_mass_release_along_plume
// [5.2.1.1.] calc_pdf_fraction
// [5.2.2.]  compute_total_released_mass
// [5.2.3.] write_cloud_trajectory_and_mass
// [5.2.4.]  get_sdimcutoff

// --- Diffusion and drift in the air ---
// [5.3.1.]  compute_falltime_and_drift_profile
// [5.3.2.]  drift_from_a_certain_source
// [5.3.3.]  calc_cloud_sigma2
// [5.3.4.]  calc_particle_terminal_velocity

// --- CPU mass loading ---
// [5.4.]  calc_mass_loading

// --- Flatten ---
// [5.4.1] idx_ps
// [5.4.2] idx_psz

// --- CUDA backend ---
// [5.5.0.]  funcD01_direct	merged funcD01a and funcD01b
// [5.5.1.]  funcD01a		compute partial mass loading for each (location, source, phidec)
// [5.5.2.]  funcD01b		reduce over (source, phidec) to obtain total per location

// --- GPU buffer / wrapper ---
// [5.6.1.]  prepare_mass_loading
// [5.6.2.]  allocate_device_buffers_struct
// [5.6.3.]  pack_fixed_data_buffers
// [5.6.4.]  copy_fixed_data_to_device_buffers
// [5.6.5.]  compute_mass_loading
// [5.6.6.]  pack_location_data_buffers
// [5.6.7.]  copy_location_data_to_device_buffers
// [5.6.8.]  launch_mass_loading_kernels_buffers
// [5.6.9.]  copy_result_to_host_buffers
// [5.6.10.] cleanup_buffers
// [5.6.11.] accumulate_massloading_for_phi

// CPU functions
// [5.5.1b] calc_mass_loading_element
// [5.5.2b] calc_mass_loading_location


/*[5.1.]
 * Perform main simulation steps for each grain-size class:
 * - fall          : particle settling through atmosphere
 * - drift         : horizontal transport by wind
 * - release       : mass release along plume axis
 * - mass loading  : deposition at ground locations
 * - accumulation  : summation over sources and grain sizes
 */
void calculate_massloading(
    int zmax,
    double *h,
    double *atmP,
    double *atmT,
    double *windX,
    double *windY,
    double *locX,
    double *locY,
    double *locZ,
    double *sourceX,
    double *sourceY,
    double *sourceZ,
    double *sourceRadius,
    double *ttlfalltime,
    double *driftX,
    double *driftY,
    double *ttlfalltime_phiint,
    double *ttldriftX_phiint,
    double *ttldriftY_phiint,
    double *ttlfalltime_phidec,
    double *ttldriftX_phidec,
    double *ttldriftY_phidec,
    double *massreleased_per_ds_and_phidec,
    SEG *massreleased_per_ds,
    double *cloud_center_x,
    double *cloud_center_y,
    double *cloud_sigma2,
    double *tmpmassloading,
    double *ttlmassloading,
    double *cummassphi,
    DEP *location_properties,
    RELEASE *r
){
    /*
     * Main computation loop over integer phi classes.
     *
     * For each grain-size class, this loop:
     * - computes fall time and drift for decimal phi intervals
     * - computes particle release along the plume axis
     * - computes drift centers and diffusion width from each source
     * - computes mass loading at ground locations (contribution of each source and decimal phi size classes)
     * - accumulates deposit information for output (mass loading of each phi size classes and total mass)
     *
     */

    char string[30];
    int phisize;
    double grainsize, phi;

    /*
	* Main loop over grain-size classes (phi) to obtain mass loading at each location.
	*
	* For each integer phi (outer loop) and its decimal subdivisions (inner loop):
	*
	* A result file (massloading.txt) shows total and 1phi-interval mass loadings but
	* mass loading calculation is done for each 0.1phi classes 
	*
	* 1. Compute fall time and horizontal drift for each grain size (F10)
	*    → ttlfalltime[z], driftX[z], driftY[z]
	*
	* 2. Store results for decimal phi classes (phidec) and integer phi classes
	*
	* 3. Compute particle release (segregation) along the plume axis
	*    → massreleased_per_ds_and_phidec
	*
	* After finishing decimal phi loop (for one integer phi):
	*
	* 4. Compute cloud center positions and dispersion from each source (F20)
	*    → cloud_center_x/y, cloud_sigma2
	*
	* 5. Output cloud-center trajectory and sum released mass over decimal phi bins
    *    to obtain per-source mass for the current integer phi interval
	*
	* 6. Compute mass loading at ground locations
	*    → tmpmassloading → accumulate into ttlmassloading
	*
	* 7. Store results in location_properties and reset temporary arrays
	*/

	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		for(int phidecimal = 0; phidecimal < PHIDECDIM; phidecimal++){

			phi = phiint + MAX_GRAINSIZE + 1 - phidecimal * INTERVAL_DECIMAL_PHI;
			//printf("INTERVAL_DECIMAL_PHI=%1.1f\tPHIDECDIM = %d\tphi = %1.1f\n", INTERVAL_DECIMAL_PHI, PHIDECDIM, phi);
			grainsize = pow(2, -phi) * 0.001;	// grain size in mm

			/*
			* 1. Compute fall time and horizontal drift for each grain size (interval_fall_class; F10)
			*
			* Compute fall time and horizontal drift for a given grain size.
			*
			* Inputs:
			*   zmax, phidecimal, grainsize
			*   h, atmP, atmT, windX, windY
			*
			* Outputs (written in-place):
			*   driftX[z]      : horizontal drift in X-direction [m]
			*   driftY[z]      : horizontal drift in Y-direction [m]
			*   ttlfalltime[z] : elapsed fall time to height z [s]
			*/
			compute_falltime_and_drift_profile(zmax, phidecimal, grainsize, h, atmP, atmT, windX, windY, driftX, driftY, ttlfalltime);
			
			// 2. Store results for decimal phi classes (phidec) and integer phi classes
			//    Map ttlfalltime[z] and driftXY[z] to ttlfalltime_phidec[z, phidec] and ttldriftXYphidec[z, phidec]
			store_profile_for_phi(phidecimal, ttlfalltime, ttlfalltime_phidec);
			store_profile_for_phi(phidecimal, driftX, ttldriftX_phidec);
			store_profile_for_phi(phidecimal, driftY, ttldriftY_phidec);
			
			/*printf("grainsize = %1.4e m\n", grainsize);
			for(int t = 0; t < ZDIM; t++){
				printf("%1.4f,%1.4f,%1.4f\n", ttlfalltime[t], driftX[t], driftY[t]);
			}*/
			// end of 20240728
			
			// output 1 phi interval fallout
			if(phidecimal == 0){
				store_profile_for_phi(phiint, ttlfalltime, ttlfalltime_phiint);
				store_profile_for_phi(phiint, driftX, ttldriftX_phiint);
				store_profile_for_phi(phiint, driftY, ttldriftY_phiint);
			}
			/* 3. Compute particle release (segregation) along the plume axis */
			compute_mass_release_along_plume(zmax, phidecimal, phi, h, atmP, atmT, windX, windY, massreleased_per_ds_and_phidec);
		}// END OF DECIMAL PHI LOOP
		
		// 4. Compute cloud center positions and dispersion from each source (F20)
		drift_from_a_certain_source(sourceX, sourceY, sourceZ, sourceRadius, ttlfalltime_phidec, ttldriftX_phidec, ttldriftY_phidec, cloud_center_x, cloud_center_y, cloud_sigma2);
		// 5. Convert particle release into per-source (s) representation
		write_cloud_trajectory_and_mass(phiint, cloud_center_x, cloud_center_y, cloud_sigma2, massreleased_per_ds_and_phidec, massreleased_per_ds);  // mass released for 1phi interval is also calculated from 0.1 phi interval data
		get_sdimcutoff(cloud_sigma2, massreleased_per_ds, phiint);	// obtain SDIMCUTOFF
		
		// 20240728 output 0.1 phi interval fallout to each 1 phi interval file
		if(WRITE_DECIMAL_FALL_TRAJ){
			sprintf(string, "decimal_falltime_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			write_vertical_profiles_for_phi(string, h, ttlfalltime_phidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
			sprintf(string, "decimal_falldriftX_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			write_vertical_profiles_for_phi(string, h, ttldriftX_phidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
			sprintf(string, "decimal_falldriftY_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			write_vertical_profiles_for_phi(string, h, ttldriftY_phidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);	
		}  // 20240728 END
		
		if(WRITE_FALL_INFO_FILES){
			sprintf(string, "particle_segregation_%1.0f.txt", phiint + MAX_GRAINSIZE + 1);
			write_phi_s_table(string, massreleased_per_ds_and_phidec, phiint + MAX_GRAINSIZE + 1, PHIDECDIM, INTERVAL_DECIMAL_PHI);
		}

		r[phiint].phi = phiint + MAX_GRAINSIZE + 1;
		r[phiint].actual = compute_total_released_mass(massreleased_per_ds_and_phidec);
		
		// print calculation status to standard output (per integer phi)
		phisize = phiint + MAX_GRAINSIZE + 1;
		if(PRINT_PROGRESS) printf("PHI = %d\tSDIMCUTOFF = %d\tMASS_RELEASED = %1.4e\n",  phisize, SDIMCUTOFF, r[phiint].actual);
		
		// 6. Compute mass loading at ground locations
		if(SDIMCUTOFF > 0){
			double *massloading_loc_source_phi; // =  massloading data for each combination of location, s (position in plume) and phi (grain size, decimal phi)
			massloading_loc_source_phi = (double*)calloc(PHIDECDIM * SDIMCUTOFF * LOCDIM, sizeof(double));
		
#ifdef CUDA
			calc_mass_loading(sourceZ, cloud_center_x, cloud_center_y, cloud_sigma2, locX, locY, locZ, massloading_loc_source_phi, tmpmassloading, massreleased_per_ds_and_phidec);
			accumulate_massloading_for_phi(phiint, tmpmassloading, ttlmassloading, cummassphi);
#else
			calc_mass_loading_element(phisize, sourceZ, cloud_center_x, cloud_center_y, cloud_sigma2, locX, locY, locZ, massloading_loc_source_phi, massreleased_per_ds_and_phidec);
			calc_mass_loading_location(phiint, massloading_loc_source_phi, tmpmassloading, ttlmassloading, cummassphi);
#endif
			
#ifdef TEST
			write_massloading_loc_source_phi(phiint, massloading_loc_source_phi);
#endif
			if(WRITE_DECIMAL_MASSLOADING){write_massloading_per_phidec_at_locations(location_properties, phiint, massloading_loc_source_phi);}
			free(massloading_loc_source_phi);
		}
		
		store_massloading_for_phi(phiint, location_properties, tmpmassloading);
		clear_array(LOCDIM, tmpmassloading);
		
	}// END OF INTEGER PHI LOOP
} /// END OF 5.1.


/* [5.1.1.] write_phi_s_table
 *
 * Store a vertical profile for a phi class.
 *
 * Copies a height-dependent quantity src[z]
 * into the storage array associated with phi.
 *
 * Used for caching intermediate profiles during
 * mass-loading calculations, including:
 *   - fall time
 *   - horizontal drift
 *   - diffusion parameters
 */
void write_phi_s_table(
    const char *filename,
    double *massrelease_by_phi_and_source,
    double base,
    int stepnum,
    double stepdelta
){
    FILE *outfile = fopen(filename, "w");

    // header
    fprintf(outfile, "s");
    for(int phiint = 0; phiint < stepnum; phiint++){
        fprintf(outfile, "\tphi%1.1f", base - phiint * stepdelta);
    }
    fprintf(outfile, "\n");

    // data
    for(int s = 0; s < SDIM_FOR_FALL_CALC; s++){
        fprintf(outfile, "%d", s);
        for(int phiint = 0; phiint < stepnum; phiint++){
            fprintf(outfile, "\t%1.4f",
                massrelease_by_phi_and_source[s + phiint * SDIM_FOR_FALL_CALC]);
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
} // end of 5.1.1.

/*
 * [5.1.2.] Mass-loading helper
 * Store a vertical 1D profile into a phi-indexed 2D table.
 *
 * Copies a height-dependent profile src[z]
 * into a flattened array dst[phi, z].
 *
 * Layout:
 *   dst[z + phi * ZDIM]
 *
 * Used to cache per-phi vertical quantities such as:
 *   - fall time
 *   - horizontal drift
 *   - diffusion parameters
 *
 * This enables later lookup by (phi, z)
 * during mass-loading calculations.
 *
 * 
 */
void store_profile_for_phi(int phi, double *src, double *dst){
    for(int z = 0; z < ZDIM; z++){
        dst[z + phi * ZDIM] = src[z];
    }
} // end of 5.1.2.

/* [5.2.1.]
 * Compute released particle mass along the plume axis for one decimal phi class.
 *
 * For a given grain size (phi), this function estimates how much particle mass
 * is released from each source interval s along the plume axis.
 *
 * The release distribution is controlled by:
 * - terminal fall velocity of the particle
 * - wind speed at plume height
 * - plume thickness
 * - grain-size PDF fraction
 *
 * Output:
 * - massreleased[phidecimal][s] : released mass from source interval s
 */
void compute_mass_release_along_plume(int zmax, int phidecimal, double phi, double *h, double *atmP, double *atmT, double *windX, double *windY, double *massreleased){
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
}	// END OF 5.2.1.

/* [5.2.1.1.]
 * Compute mass fraction for a decimal phi bin.
 *
 * The grain-size distribution is approximated by a normal distribution
 * in phi scale with:
 *   - MEDIAN_GRAINSIZE as the mean
 *   - STD_GRAINSIZE as the standard deviation
 *
 * The probability density is evaluated at the center of the phi bin:
 *
 *   phi_center = phi - INTERVAL_DECIMAL_PHI / 2
 *
 * and multiplied by the bin width INTERVAL_DECIMAL_PHI.
 */
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
} // end of 5.2.1.1.

/* [5.2.2.]
 * Compute total released mass from the plume for the current phi class.
 *
 * Sum massreleased_per_ds_and_phidec over all (s, phidec),
 * but only for source indices s < SDIMCUTOFF.
 */
double compute_total_released_mass(double *massreleased_per_ds_and_phidec){
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
}	// END OF 5.2.2.

/* [5.2.3.]
 * Write cloud-center trajectory and summarize released mass for one integer phi class.
 *
 * For each source interval s along the plume axis, this function:
 * - sums released mass over decimal phi bins within the current integer phi class
 * - stores the summed mass in massreleased_per_ds[s].mass_from_ds[phiint]
 * - optionally writes cloud-center coordinates and dispersion to depcenttraj*.txt
 *
 * Inputs:
 * - x, y   : cloud-center coordinates indexed by (phidec, s, z)
 * - sig    : cloud dispersion variance indexed by (phidec, s, z)
 * - r      : released mass indexed by (phidec, s)
 *
 * Note:
 * The output trajectory uses z = 0, i.e. the cloud center at ground level.
 */
void write_cloud_trajectory_and_mass(int phiint, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *massreleased_per_ds_and_phidec, SEG *massreleased_per_ds){
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
			fprintf(outfile, "i\tx_from_vent\ty_from_vent\tx_coord\ty_coord\tcloud_sigma2\tmass_released\n");
		}

		for(int i = 0; i < SDIM_FOR_FALL_CALC; i++){
			released = 0.0;
			for(int gsize = 0; gsize < PHIDECDIM; gsize++){
				released += massreleased_per_ds_and_phidec[i + gsize * SDIM_FOR_FALL_CALC];
			}
			
			massreleased_per_ds[i].mass_from_ds[phiint] = released; 
			
			idz = i * ZDIM;
			if(WRITE_DEPCENT_TRAJECTORY) fprintf(outfile, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4e\n", i, cloud_center_x[idz], cloud_center_y[idz], cloud_center_x[idz] + VENT_EASTING, cloud_center_y[idz] + VENT_NORTHING, cloud_sigma2[idz], released);
		}
		if(WRITE_DEPCENT_TRAJECTORY) fclose(outfile);
} // End of 5.2.3.

/* [5.2.4.]
 * Determine SDIMCUTOFF for the current integer phi class.
 *
 * SDIMCUTOFF is the effective upper limit of source points used in
 * mass-loading calculation. Source points farther than this cutoff are
 * ignored when their estimated contribution becomes smaller than the
 * configured minimum threshold.
 *
 * The contribution is estimated from:
 *   released mass from source s / cloud_sigma2 at ground level
 *
 * Inputs:
 * - cloud_sigma2[phidec][s][z]          : cloud dispersion variance
 * - massreleased_per_ds[s].mass_from_ds : released mass per source and phi
 */
void get_sdimcutoff(
    double *cloud_sigma2,
    SEG *massreleased_per_ds,
    int phiint
){
    double estimated_contribution;

    for(int s = 0; s < SDIM_FOR_FALL_CALC; s++){
        /*
         * Use phidec = 0 and z = 0 as representative values
         * for the current source point s.
         */
        estimated_contribution =
            massreleased_per_ds[s].mass_from_ds[phiint]
            / cloud_sigma2[s * ZDIM];

        if(estimated_contribution < MINIMUM_CONTRIBUTION * S_DELTA_FOR_FALL_CALC){
            SDIMCUTOFF = s;
            break;
        }
    }
}	// END OF 5.2.4.

/* [5.3.1.]
 * Compute fall time and wind drift for a particle of a given grain size.
 *
 * Starting from plume height Ht (= h[zmax]), this function integrates
 * particle settling downward through each vertical interval.
 *
 * Outputs:
 * - ttlfalltime[z] : elapsed fall time from Ht to height level z [s]
 * - driftX[z]      : accumulated X-direction drift from Ht to height level z [m]
 * - driftY[z]      : accumulated Y-direction drift from Ht to height level z [m]
 *
 * Notes:
 * - The top interval is h[zmax] → h[zmax-1].
 * - Lower intervals have thickness Z_DELTA.
 * - In the default formulation, fall velocity and wind velocity are averaged
 *   between the top and bottom of each interval.
 * - Under TEPHRA2, the upper-level wind and terminal velocity are used.
 */
void compute_falltime_and_drift_profile(int zmax, int phidecimal, double grainsize, double *h, double *atmP, double *atmT, double *windX, double *windY, double *driftX, double *driftY, double *ttlfalltime){
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
} // END OF 5.3.1.

/* [5.3.2.]
 * drift_from_a_certain_source (F20)
 * Compute the center position and dispersion of a particle cloud
 * released from each point along the plume axis.
 *
 * For each source point s and height interval z, this function calculates:
 * - cloud_center_x[phidec][s][z]     : X-coordinate of cloud center
 *                                (wind drift + transport along the plume)
 * - cloud_center_y[phidec][s][z]     : Y-coordinate of cloud center
 *                                (wind drift + transport along the plume)
 * - cloud_sigma2[phidec][s][z] : variance of horizontal dispersion of the cloud
 *
 * A "cloud" is a group of particles that:
 * - share the same grain size (same phidec)
 * - are released from the same source point s
 */
void drift_from_a_certain_source(double *source_x, double *source_y, double *source_height, double *sourceRadius, double *TotalFallTime, double *driftX, double *driftY, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2){
	int s, z, phidec, idz, idz_s;
	int z_source;		// z_source means interval count of z axis of just above the source height
	double residue_up, fall_time_residue_up;
	double falltime, ttldriftX, ttldriftY;

	for(int idx = 0; idx < ZDIM * SDIM_FOR_FALL_CALC * PHIDECDIM; idx++){		// idx is count for driftXY_s and cloud_sigma2
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

			cloud_center_x[idx] = source_x[s] + ttldriftX;
			cloud_center_y[idx] = source_y[s] + ttldriftY;
			cloud_sigma2[idx] = calc_cloud_sigma2(sourceRadius[s] * PLUME_RADIUS_CORRECTION, falltime); // F21
		}
	}
} // End of the function (F20) [5.3.2.]

/* [5.3.3.]
 * Compute variance (sigma^2) of particle cloud dispersion.
 *
 * This function estimates the horizontal spread of a particle cloud
 * originating from a source of finite radius, based on fall time.
 *
 * Two regimes are used following Bonadonna et al. (2005):
 *
 * - Coarse particles (short fall time):
 *     Diffusion-dominated spreading (Eq. 6)
 *
 * - Fine particles (long fall time):
 *     Turbulent eddy diffusion dominates (Eq. 8)
 *
 * Inputs:
 *   source_radius :
 *     radius of the source (plume cross-section) [m]
 *
 *   falltime :
 *     particle fall time from release height [s]
 *
 * Output:
 *   cloud_sigma2 :
 *     variance (sigma^2) of horizontal dispersion [m^2]
 *
 * Notes:
 * - virtual_falltime represents the time required for a point source
 *   to spread to the initial plume radius.
 * - FALL_TIME_THRESHOLD separates coarse and fine particle regimes.
 */
double calc_cloud_sigma2(double source_radius, double falltime) {
	// time needed for point source to diffuse until sigma equals to the plume radius
	double virtual_falltime = 0.0;
	double cloud_sigma2 = 0.0;

	if (falltime < FALL_TIME_THRESHOLD){
		// coarse particle	Bonadonna+(2005) Eq.6
		virtual_falltime = source_radius * source_radius / (4 * DIFFUSION_COEFFICIENT);
		cloud_sigma2 = 4 * DIFFUSION_COEFFICIENT * (falltime + virtual_falltime);
	}else{
		// fine particle	Bonadonna+(2005) Eq.8
		virtual_falltime = pow((5 * source_radius * source_radius) / (8 * EDDY_CONST), 0.4);
		cloud_sigma2 = 8 * EDDY_CONST / 5 * pow(falltime + virtual_falltime, 2.5);
	}
	return(cloud_sigma2);
} // END OF 5.3.3.


// [5.3.4.]
/* Compute particle terminal fall velocity based on Reynolds-number-dependent drag regimes */
double calc_particle_terminal_velocity(double h, double ashdiam, double part_density, double p, double t) {
  // Modified from function “particle_fall_time” in tephra_calc.c of tephra2
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
} // end of 5.3.4.

#ifdef CUDA
/* [5.4.]
 * Compute mass loading at all ground locations using CUDA.
 * This function manages the complete GPU workflow for mass-loading calculation.
 *
 * Inputs
 * -------
 * sourceZ
 *     Height of each particle source along the plume axis.
 *
 * cloud_center_x, cloud_center_y
 *     Horizontal coordinates of the particle cloud center
 *     for every (phidec, source, height).
 *
 * cloud_sigma2
 *     Horizontal variance (σ²) of the particle cloud
 *     for every (phidec, source, height).
 *
 * locX, locY, locZ
 *     Coordinates of all ground locations.
 *
 * massreleased
 *     Released particle mass for every (phidec, source).
 *
 * Outputs
 * --------
 * ttlml
 *     Total mass loading at each ground location.
 *
 * massloading_loc_source_phi
 *     (Optional / diagnostic)
 *     Mass loading for each combination of
 *     (location, source, decimal-phi).
 *
 * Notes
 * -----
 * Because the number of ground locations (LOCDIM) may exceed available
 * GPU memory, locations are processed in chunks. *******Fixed plume-related
 * data are copied to the GPU once for each phi interval, whereas
 * location data are transferred for each chunk.
 */
void calc_mass_loading(double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *locX, double *locY, double *locZ, double *massloading_loc_source_phi, double *ttlml, double *massreleased){
	/* 
	* D in the name of parameter (e.g. ttlmlD) comes from "Device", which means such parameters
	* are used in GPU calculation
	* F in the name of parameter (e.g. ttlmlF) means such parameters are temporaly ones in CPU
	*/

	int chunk_locdim = 8192;  // Empirically tuned on NVIDIA GeForce RTX 3060; output verified by diff

	/* 1. Define size of arrays used in GPU */
	/*
	* PSZ : size of arrays indexed by (phidec, source, height_interval)
	* LSP : size of arrays indexed by (location, phidec, source)
	*/
	size_t LSP;
	size_t PSZ; 
	
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

	PSZ = PHIDECDIM * SDIMCUTOFF * ZDIM;	//PSZ = PHIDECDIM * SDIM_FOR_FALL_CALC* ZDIM;
	LSP = (size_t)chunk_locdim * SDIMCUTOFF * PHIDECDIM;	//LSP = LOCDIM * SDIM_FOR_FALL_CALC* PHIDECDIM;

	/* end of 1.*/

	/* 2. Host Buffer Allocation */

	float *ttlmlF;
	ttlmlF = (float *)malloc(chunk_locdim * sizeof(float));
	if (!ttlmlF) {
    fprintf(stderr, "Error: malloc failed for ttlmlF\n");
    exit(EXIT_FAILURE);
	}

	// Allocate host buffers used for double-to-float conversion
	float *sourceZF, *centXF, *centYF, *sigsqF, *locXF, *locYF, *locZF, *massreleasedF;
	sourceZF = (float *)malloc(SDIMCUTOFF * sizeof(float));
	centXF = (float *)malloc(PSZ * sizeof(float));
	centYF = (float *)malloc(PSZ * sizeof(float));
	sigsqF = (float *)malloc(PSZ * sizeof(float));

	locXF = (float *)malloc(chunk_locdim * sizeof(float));
	locYF = (float *)malloc(chunk_locdim * sizeof(float));
	locZF = (float *)malloc(chunk_locdim * sizeof(float));

	massreleasedF = (float *)malloc(PHIDECDIM * SDIMCUTOFF * sizeof(float));
	
	if (!ttlmlF || !sourceZF || !centXF || !centYF || !sigsqF ||
    !locXF || !locYF || !locZF || !massreleasedF) {
    fprintf(stderr, "Error: malloc failed for host buffers\n");
    exit(EXIT_FAILURE);
	}	

	/* end of 2.*/

	/* 3-4. Device Buffer Allocation */
	Buffers b;

	b.PSZ = PSZ;
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

	prepare_mass_loading(	//Steps 3-5
		&b,
		sourceZ,
		cloud_center_x,
		cloud_center_y,
		cloud_sigma2,
		massreleased,
    	chunk_locdim
	);
	
	/* Loop for compute mass loading */
	for (int loc0 = 0; loc0 < LOCDIM; loc0 += chunk_locdim) {

		int locN = chunk_locdim;
		if (loc0 + locN > LOCDIM) {
			locN = LOCDIM - loc0;
		}

		/*
		* Step 6 of calc_mass_loading:
		* Per-chunk GPU processing pipeline.
		*/
		compute_mass_loading(&b, locX, locY, locZ, ttlml, loc0, locN);
	}
		/* Step 7 Clean up*/
		cleanup_buffers(&b);
	/* end of host bridge */

} // END OF 5.4.
#endif

#ifdef CUDA
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
/* [5.4.1.]
 * Flatten 3D index (phidec, s, z) into 1D array index.
 *
 * Data layout:
 *   data[phidec][s][z] is stored as a contiguous 1D array:
 *     index = phidec * (sdim * zdim) + s * zdim + z
 *
 * Dimensions:
 *   phidec : grain size class (decimal phi)
 *   s      : source index along plume
 *   z      : vertical level
 */
 
static HOST_DEVICE inline int
idx_psz(
	int phidec, int s, int z, int sdim, int zdim)
{
    return phidec * sdim * zdim + s * zdim + z;
} // END OF 5.4.1.

/* [5.4.2.]
 * Flatten 2D index (phidec, s) into 1D array index.
 *
 * Data layout:
 *   data[phidec][s] is stored as:
 *     index = phidec * sdim + s
 *
 * Dimensions:
 *   phidec : grain size class (decimal phi)
 *   s      : source index along plume
 */
static HOST_DEVICE inline int
idx_ps(
    int phidec,
    int s,
    int sdim
){
    return phidec * sdim + s;
}// END OF 5.4.2.

#ifdef CUDA
/* [5.5.0.]
 * merged funcD01a and funcD01b to calculate massloading for a certain phi class of 1 phi interval for location
 * direct: directly computes ttlmlD without storing massloading_loc_source_phiD
*/

__global__ void funcD01_direct(
    int locdim,
    int zdim,
    int sdim,
    int phidecdim,
    float zdelta,
    float *ttlmlD,
    float *sourceZD,
    float *centX,
    float *centY,
    float *cloud_sigma2,
    float *locX,
    float *locY,
    float *locZ,
    float *massreleased
){
    int j = blockIdx.x;      // 1 block = 1 location
    int t = threadIdx.x;

    extern __shared__ float shsum[];

    float local_sum = 0.0f;
    int n = sdim * phidecdim;

    for (int k = t; k < n; k += blockDim.x) {
        int s = k / phidecdim;
        int phidec = k % phidecdim;

        if (locZ[j] >= sourceZD[s]) continue;

        int z = (int)(locZ[j] / zdelta);

        int ipsz = idx_psz(phidec, s, z, sdim, zdim);
        int ips  = idx_ps(phidec, s, sdim);

        float depcentX =
            centX[ipsz + 1]
            + (centX[ipsz] - centX[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        float depcentY =
            centY[ipsz + 1]
            + (centY[ipsz] - centY[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        float sigma2 =
            cloud_sigma2[ipsz + 1]
            + (cloud_sigma2[ipsz] - cloud_sigma2[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        float square_distance =
            (depcentX - locX[j]) * (depcentX - locX[j])
            + (depcentY - locY[j]) * (depcentY - locY[j]);

#ifndef TEPHRA2
        local_sum +=
            1.0f / (M_2PI * sigma2)
            * exp(-square_distance / (2.0f * sigma2))
            * massreleased[ips];
#else
        local_sum +=
            1.0f / (M_PI * sigma2)
            * exp(-square_distance / sigma2)
            * massreleased[ips];
#endif
    }

    shsum[t] = local_sum;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            shsum[t] += shsum[t + stride];
        }
        __syncthreads();
    }

    if (t == 0) {
        ttlmlD[j] = shsum[0];
    }
}

/**
/* [5.5.1.]
 * @brief Compute partial mass loading for each (location, source, phi) element (= lsmpl).
 *
 * Each CUDA thread computes one flattened element of:
 *
 *     massloading_loc_source_phiD[location, source, phidec]
 *
 * Flattened layout:
 *
 *     tid = ((location * sdim) + source) * phidecdim + phidec
 *
 * @param N Total number of flattened elements:
 *          locdim * sdim * phidecdim
 * @param zdim Number of vertical grid intervals
 * @param sdim Number of source points
 * @param phidecdim Number of grain-size subdivisions
 * @param zdelta Vertical grid spacing
 * @param massloading_loc_source_phiD Output partial mass loading per location/source/phi
 * @param ttlmlD Output total mass loading buffer, used by funcD01b
 * @param sourceZD Source height array
 * @param centX Deposit center X array indexed by (phidec, source, z)
 * @param centY Deposit center Y array indexed by (phidec, source, z)
 * @param cloud_sigma2 Variance array indexed by (phidec, source, z)
 * @param locX Location X array for the current chunk
 * @param locY Location Y array for the current chunk
 * @param locZ Location Z array for the current chunk
 * @param massreleased Released mass indexed by (phidec, source)
 */
__global__ void funcD01a(
    int N,
    int zdim,
    int sdim,
    int phidecdim,
    float zdelta,
    float *massloading_loc_source_phiD,
    float *ttlmlD,
    float *sourceZD,
    float *centX,
    float *centY,
    float *cloud_sigma2,
    float *locX,
    float *locY,
    float *locZ,
    float *massreleased
){
    int j, s, z, phidec, ips, ipsz;
    float depcentX, depcentY, sigma2, square_distance;

    unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if(tid < N){
        massloading_loc_source_phiD[tid] = 0.0f;

        /*
         * Decode flattened thread index.
         *
         * Layout:
         *   tid = ((j * sdim) + s) * phidecdim + phidec
         *
         * j      : location index within the current chunk
         * s      : source index
         * phidec : grain-size subdivision index
         */
        j      = tid / (phidecdim * sdim);
        s      = (tid / phidecdim) % sdim;
        phidec = tid % phidecdim;

		/*
		* Determine vertical layer index z such that
		*   z * zdelta <= locZ[j] < (z + 1) * zdelta
		
         * Arrays centX, centY, and cloud_sigma2 are interpolated
         * between z and z+1.
         */
        z = locZ[j] / zdelta;

        /*
         * ipsz indexes arrays flattened from:
         *   [phidec][source][z]
         *
         * ips indexes arrays flattened from:
         *   [phidec][source]
         */
		 
		ipsz = idx_psz(phidec, s, z, sdim, zdim);
		ips = idx_ps(phidec, s, sdim);

        depcentX =
            centX[ipsz + 1]
            + (centX[ipsz] - centX[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        depcentY =
            centY[ipsz + 1]
            + (centY[ipsz] - centY[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        sigma2 =
            cloud_sigma2[ipsz + 1]
            + (cloud_sigma2[ipsz] - cloud_sigma2[ipsz + 1])
            * (zdelta * (z + 1) - locZ[j]) / zdelta;

        square_distance =
            (depcentX - locX[j]) * (depcentX - locX[j]) 
            + (depcentY - locY[j]) * (depcentY - locY[j]);
		/*square_distance =
            (depcentX - locX[j]) * (depcentX - locX[j])
            + (depcentY - locY[j]) * (depcentY - locY[j]);*/

        /*
         * Only sources above the current location (locZ[j]) contribute
         * to mass loading at that location.
         */
        if(locZ[j] < sourceZD[s]){
            /* Original formulation: Bonadonna et al. (2005) */
            massloading_loc_source_phiD[tid] =
                1 / (M_2PI * sigma2)
                * exp(-square_distance / (2 * sigma2))
                * massreleased[ips];

#ifdef TEPHRA2
            /* Formulation used in Tephra2 and WT */
            massloading_loc_source_phiD[tid] =
                1 / (M_PI * sigma2)
                * exp(-square_distance / sigma2)
                * massreleased[ips];
#endif
        }

        /*
         * Note:
         * This increment has no effect unless this if-block is changed
         * to a while-loop. It is kept here to avoid changing behavior.
         */
        tid += blockDim.x * gridDim.x;
    }
} // END OF 5.5.1.


/**
 * [5.5.2]
 * @brief Reduce partial mass loading over source and phi for each location.
 *
 * funcD01a produces:
 *
 *     massloading_loc_source_phiD[location, source, phidec]
 *
 * This kernel sums all source/phi contributions for each location:
 *
 *     ttlmlD[location] = sum over source and phidec
 *
 * @param N Total number of flattened massloading_loc_source_phiD elements:
 *          locdim * sdim * phidecdim
 * @param locdim Number of locations in the current chunk
 * @param massloading_loc_source_phiD Partial mass loading array
 * @param ttlmlD Output total mass loading per location
 */
__global__ void funcD01b(
    int N,
    int locdim,
    float *massloading_loc_source_phiD,
    float *ttlmlD
){
    unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if(tid < locdim){
        int n_per_location = N / locdim;

        ttlmlD[tid] = 0.0f;

        for(int i = 0; i < n_per_location; i++){
            ttlmlD[tid] += massloading_loc_source_phiD[tid * n_per_location + i];
        }
    }
} // END OF 5.5.2.

/* [5.6.1.] */
/* Step 3-5 of calc_mass_loading:
 * Preparation stage of calc_mass_loading.
 * Allocate device buffers and initialize fixed data.
 *
 * This stage performs:
 *   Step 3: allocate device (GPU) buffers (structures)
 *   Step 4: pack and copy fixed data to device
 *   Step 5: allocate device (GPU) buffers (fixed data)
 *
 * This function is called once before chunk-based processing.
 */
static void prepare_mass_loading(
    Buffers *b,
    double *sourceZ,
    double *cloud_center_x,
    double *cloud_center_y,
    double *cloud_sigma2,
    double *massreleased,
    int chunk_locdim){
	/* Step 3 of calc_mass_loading */
    allocate_device_buffers_struct(b, chunk_locdim);

	/* Step 4 of calc_mass_loading */
    pack_fixed_data_buffers(
        b,
        sourceZ,
        cloud_center_x,
        cloud_center_y,
        cloud_sigma2,
        massreleased
    );
	/* Step 5 of calc_mass_loading */
    copy_fixed_data_to_device_buffers(b);
}

/***************************
/* [5.6.2.]
 * Step 3 of calc_mass_loading:
 * Allocate device (GPU) buffers.
 *
 * This step allocates all GPU memory required for mass-loading
 * calculation, based on problem size and chunking strategy.
 *
 * Memory layout:
 *
 *   massloading_loc_source_phiD :
 *     (location, source, phidec)
 *       = LSP = LOCDIM * SDIM * PHIDECDIM
 *
 *   ttlmlD :
 *     total mass loading per location (current chunk)
 *       = chunk_locdim
 *
 *   sourceZD :
 *     source height for each plume segment (s)
 *       = SDIMCUTOFF
 *
 *   centXD, centYD, sigsqD :
 *     cloud properties indexed by (phidec, s, z)
 *       = PSZ = PHIDECDIM * SDIM * ZDIM
 *
 *   locXD, locYD, locZD :
 *     location coordinates for current chunk
 *       = chunk_locdim
 *
 *   massreleasedD :
 *     released mass per (phidec, s)
 *       = PHIDECDIM * SDIMCUTOFF
 *
 * Notes:
 * - chunk_locdim controls the number of locations processed per kernel launch.
 * - PSZ and LSP must be initialized before calling this function.
 */
static void allocate_device_buffers_struct(Buffers *b, int chunk_locdim)
{
    CUDA_CHECK(cudaMalloc((void**)&b->device.massloading_loc_source_phiD,
        b->LSP * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.ttlmlD,
        chunk_locdim * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.sourceZD,
        SDIMCUTOFF * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.centXD,
        b->PSZ * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.centYD,
        b->PSZ * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.sigsqD,
        b->PSZ * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locXD,
        chunk_locdim * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locYD,
        chunk_locdim * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.locZD,
        chunk_locdim * sizeof(float)));

    CUDA_CHECK(cudaMalloc((void**)&b->device.massreleasedD,
        SDIMCUTOFF * PHIDECDIM * sizeof(float)));
}
/* end of step 3*/

/***************************
/* [5.6.3.]
 * Step 4 of calc_mass_loading:
 * Pack fixed data into host buffers (CPU → GPU preparation).
 *
 * This step converts input data from double precision (CPU-side)
 * to single precision (float) and packs them into contiguous
 * host buffers for efficient transfer to GPU.
 *
 * Data packed in this step are invariant across location chunks:
 *
 *   sourceZF :
 *     source height along plume (s dimension)
 *       size = SDIMCUTOFF
 *
 *   centXF, centYF, sigsqF :
 *     cloud properties indexed by (phidec, s, z)
 *       layout = [phidec][s][z] → flattened to 1D
 *       size = PSZ = PHIDECDIM * SDIM * ZDIM
 *
 *   massreleasedF :
 *     released mass per (phidec, s)
 *       layout = [phidec][s]
 *       size = PHIDECDIM * SDIMCUTOFF
 *
 * Notes:
 * - CUDA device memory is linear; multi-dimensional arrays are flattened.
 * - Conversion to float reduces memory footprint and improves GPU throughput.
 * - These buffers are reused across chunks (only location data changes per chunk).
 */
static void pack_fixed_data_buffers(
    Buffers *b,
    double *sourceZ,
    double *cloud_center_x,
    double *cloud_center_y,
    double *cloud_sigma2,
    double *massreleased
){                                                 
	/*
	* CUDA device memory is handled as linear memory.
	* Therefore, multi-dimensional indices are flattened into
	* one-dimensional array indices before data transfer to GPU.
	*/

    for(int i = 0; i < SDIMCUTOFF; i++){
        b->host.sourceZF[i] = (float)sourceZ[i];
    }

    for(size_t ipszc = 0; ipszc < b->PSZ; ipszc++){
        b->host.centXF[ipszc] = (float)cloud_center_x[ipszc];
        b->host.centYF[ipszc] = (float)cloud_center_y[ipszc];
        b->host.sigsqF[ipszc] = (float)cloud_sigma2[ipszc];
    }

    for(int ipsc = 0; ipsc < SDIMCUTOFF * PHIDECDIM; ipsc++){
        b->host.massreleasedF[ipsc] = (float)massreleased[ipsc];
    }
}
/* End of Step 4*/


/**************************
/* [5.6.4.]
 * Step 5 of calc_mass_loading:
 * Copy fixed host buffers to device buffers.
 *
 * This step transfers data that are reused for all location chunks:
 * - source height
 * - cloud center positions
 * - cloud dispersion variance
 * - released mass per (phidec, s)
 *
 * These arrays are copied only once before the chunk loop.
 * Location data are copied separately for each chunk.
 */
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
		b->PSZ * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.centYD,
		b->host.centYF,
		b->PSZ * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.sigsqD,
		b->host.sigsqF,
		b->PSZ * sizeof(float),
		cudaMemcpyHostToDevice
	));

	CUDA_CHECK(cudaMemcpy(
		b->device.massreleasedD,
		b->host.massreleasedF,
		SDIMCUTOFF * PHIDECDIM * sizeof(float),
		cudaMemcpyHostToDevice
	));
}
/* end of Step 5 */


/* 5.6.5. */
/* Step 6 of calc_mass_loading:
 * Compute mass loading for a chunk of locations using GPU.
 *
 * This function executes the per-chunk processing pipeline:
 *
 *   Step 6a: pack location data (CPU, double → float)
 *   Step 6b: copy location data to GPU
 *   Step 6c: compute mass loading on GPU (CUDA kernels)
 *   Step 6d: copy results back to CPU and store them
 *
 * Inputs:
 *   locX, locY, locZ :
 *     full arrays of location coordinates
 *
 *   ttlml :
 *     output array of total mass loading per location
 *
 *   loc0 :
 *     starting index of the current chunk
 *
 *   locN :
 *     number of locations in this chunk
 *
 * Notes:
 * - Fixed data (plume properties, released mass, etc.) must already
 *   be copied to the device before calling this function.
 * - This function operates only on a subset (chunk) of locations.
 */
static void compute_mass_loading(
    Buffers *b,
    double *locX,
    double *locY,
    double *locZ,
    double *ttlml,
    int loc0,
    int locN
){
	/* Per-chunk GPU processing pipeline */
	pack_location_data_buffers(b, locX, locY, locZ, loc0, locN);  // Step 6a: pack location data
	copy_location_data_to_device_buffers(b, locN);                // Step 6b: copy location data to GPU
	launch_mass_loading_kernels_buffers(b, locN);                 // Step 6c: compute mass loading on GPU
	copy_result_to_host_buffers(b, ttlml, loc0, locN);            // Step 6d: copy results back to host
}


/***************************
/* [5.6.7.]
 * Step 6a (1/4) of calc_mass_loading:
 * Pack location data into host buffers for the current chunk.
 *
 * This step extracts a subset of location coordinates from the full
 * location arrays and converts them from double to float for GPU use.
 *
 * Inputs:
 *   locX, locY, locZ :
 *     full arrays of location coordinates (size = LOCDIM)
 *
 *   loc0 :
 *     starting index of the current chunk
 *
 *   locN :
 *     number of locations in this chunk (<= chunk_locdim)
 *
 * Output:
 *   locXF, locYF, locZF :
 *     packed location coordinates (size = locN)
 *
 * Notes:
 * - Only location data change per chunk; other data are reused.
 * - Conversion to float reduces memory transfer and improves performance.
 */
static void pack_location_data_buffers(
    Buffers *b,
    double *locX,
    double *locY,
    double *locZ,
    int loc0,
    int locN
){
    for(int j = 0; j < locN; j++){
        b->host.locXF[j] = (float)locX[loc0 + j];
        b->host.locYF[j] = (float)locY[loc0 + j];
        b->host.locZF[j] = (float)locZ[loc0 + j];
    }
}
/* End of Step 6a */

/**************************
/* [5.6.7.]
 * Step 6b (2/4) of calc_mass_loading:
 * Copy location data for the current chunk from host to device.
 *
 * This step transfers location coordinates (X, Y, Z) for a subset
 * of locations processed in the current chunk.
 *
 * Inputs:
 *   locN :
 *     number of locations in this chunk (<= chunk_locdim)
 *
 * Notes:
 * - This transfer is performed for each chunk.
 * - Only location data change per chunk; other data are already on device.
 * - Keeping this separate from fixed-data transfer (Step 6) avoids
 *   redundant memory copies and improves performance.
 */
static void copy_location_data_to_device_buffers(Buffers *b, int locN)
{
    CUDA_CHECK(cudaMemcpy(
        b->device.locXD,
        b->host.locXF,
        locN * sizeof(float),
        cudaMemcpyHostToDevice
    ));

    CUDA_CHECK(cudaMemcpy(
        b->device.locYD,
        b->host.locYF,
        locN * sizeof(float),
        cudaMemcpyHostToDevice
    ));

    CUDA_CHECK(cudaMemcpy(
        b->device.locZD,
        b->host.locZF,
        locN * sizeof(float),
        cudaMemcpyHostToDevice
    ));
}
/* end of Step 6b*/

/**************************
 * [5.6.8.] 
 * Step 6c of calc_mass_loading:
 * Launch CUDA kernels for the current location chunk.
 *
 * In normal mode, funcD01_direct computes total mass loading
 * directly for each location using block-level reduction.
 *
 * In diagnostic mode, funcD01a and funcD01b are used to retain
 * massloading_loc_source_phiD for decimal-phi/source debugging output.
 *
 * 1. funcD01a:
 *    Compute partial mass loading for each flattened element:
 *      (location, source, phidec)
 *
 * 2. funcD01b:
 *    Reduce partial values over (source, phidec)
 *    to obtain total mass loading for each location.
 *
 * Inputs:
 *   locN :
 *     number of locations in the current chunk
 *
 * Notes:
 * - LSP_chunk is the number of (location, source, phidec) elements
 *   processed in this chunk.
 * - N is passed as int because CUDA kernels currently use int indexing.
 * - blocksize controls CUDA thread grouping; chunk_locdim/locN controls
 *   the amount of location data processed per launch.
 */
static void launch_mass_loading_kernels_buffers(
    Buffers *b,
    int locN
)
{
    size_t LSP_chunk = (size_t)locN * SDIMCUTOFF * PHIDECDIM;

    if (LSP_chunk > INT_MAX) {
        fprintf(stderr, "Error: LSP_chunk too large: %zu\n", LSP_chunk);
        exit(EXIT_FAILURE);
    }

    int N = (int)LSP_chunk;

    int blocksize = 128; /* conventional CUDA block size */
    dim3 block(blocksize, 1, 1);
    dim3 grid((N + block.x - 1) / block.x, 1, 1);

	if (!(WRITE_DECIMAL_MASSLOADING)) {
		int blocksize = 128;
		dim3 block(blocksize);
		dim3 grid_reduce(locN);
		size_t shmem = blocksize * sizeof(float);

		funcD01_direct<<<grid_reduce, block, shmem>>>(
			locN, ZDIM, SDIMCUTOFF, PHIDECDIM, (float)Z_DELTA,
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

	} else {
		funcD01a<<<grid, block>>>(
			N, ZDIM, SDIMCUTOFF, PHIDECDIM, (float)Z_DELTA,
			b->device.massloading_loc_source_phiD,
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

		dim3 grid2((locN + block.x - 1) / block.x, 1, 1);

		funcD01b<<<grid2, block>>>(
			N,
			locN,
			b->device.massloading_loc_source_phiD,
			b->device.ttlmlD
		);
		CUDA_KERNEL_CHECK();
	}

CUDA_CHECK(cudaDeviceSynchronize());
}
/* end of Step 6c*/

/**************************
 * [5.6.9]
 * Step 6d (4/4) of calc_mass_loading (per chunk):
 * Copy computed results from device to host and store them.
 *
 * This step retrieves total mass loading per location from the GPU,
 * converts data from float (device) to double (CPU), and writes them
 * back to the global output array.
 *
 * Inputs:
 *   loc0 :
 *     starting index of the current chunk in the full location array
 *
 *   locN :
 *     number of locations in this chunk
 *
 * Notes:
 * - This is the final step of per-chunk processing.
 * - Device results (ttlmlD) are stored in float for performance,
 *   and converted back to double for CPU-side consistency.
 */
static void copy_result_to_host_buffers(Buffers *b, double *ttlml, int loc0, int locN)
{
    CUDA_CHECK(cudaMemcpy(
        b->host.ttlmlF,
        b->device.ttlmlD,
        locN * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    for(int j = 0; j < locN; j++){
        ttlml[loc0 + j] = (double)b->host.ttlmlF[j];
    }
}
/* end of Step 6d */

/* [5.6.10.]
 * Step 7 of calc_mass_loading:
 * Cleanup buffers.
 *
 * This step releases all host and device memory allocated
 * during the mass-loading calculation.
 */
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

	CUDA_CHECK(cudaFree(b->device.massloading_loc_source_phiD));
	CUDA_CHECK(cudaFree(b->device.ttlmlD));
}
/* end of Step 7 */


/* 
 * [5.6.11]
 * Accumulate mass loading for the current integer phi class.
 *
 * For each ground location j:
 * - add current-phi mass loading to total mass loading
 * - accumulate mass-weighted phi (tmp * phi) for later mean grain-size calculation
 */
void accumulate_massloading_for_phi(int phiint, double *tmp, double *ttl, double *cummassphi){
	double phi;

	for(int j = 0; j < LOCDIM; j++){
			phi = phiint + MAX_GRAINSIZE + 1;
			ttl[j] += tmp[j];
			cummassphi[j] += tmp[j] * phi;
	}
}
#endif

/* NON CUDA FUNCTIONS (START)*/
//
// [5.5.1b.]
//
// D01a
// Calculate mass loading of a certain grain size on a certain point on the ground (Sloc) from a certain source: Sloc(phi, s)
void calc_mass_loading_element(int phisize, double *sourceZ, double *cloud_center_x, double *cloud_center_y, double *cloud_sigma2, double *locX, double *locY, double *locZ, double *massloading_loc_source_phi, double *massreleased){
	int j, s, phidec, z;
	int ips;    // counter for arrays having grainsize(phidec) - plumelength(s; non cut off) order such as massreleased
	int ipsz;	// counter for arrays having grainsize(phidec) - plumelength(s; non cut off) - height(z) order such as cloud_center_x 
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
		ipsz = idx_psz(phidec, s, z, SDIM_FOR_FALL_CALC, ZDIM);
		//ips = s + phidec * SDIM_FOR_FALL_CALC;

		ips = idx_ps(phidec, s, SDIM_FOR_FALL_CALC);

		if(locZ[j] < sourceZ[s]){
			depcentX = cloud_center_x[ipsz+1] + (cloud_center_x[ipsz] - cloud_center_x[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			depcentY = cloud_center_y[ipsz+1] + (cloud_center_y[ipsz] - cloud_center_y[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			sigma2 = cloud_sigma2[ipsz+1] + (cloud_sigma2[ipsz] - cloud_sigma2[ipsz+1]) * (Z_DELTA * (z + 1) - locZ[j]) / Z_DELTA;
			square_distance = pow((depcentX - locX[j]), 2) + pow((depcentY - locY[j]), 2);
			
			massloading_loc_source_phi[idx] = 1 / (M_2PI * sigma2) * exp(-square_distance / (2 * sigma2)) * massreleased[ips];
#ifdef TEPHRA2
			massloading_loc_source_phi[idx] = 1 / (M_PI * sigma2) * exp(-square_distance / (sigma2)) * massreleased[ips];	// Formulation used in Tephra2 and WT
#endif
			//if(j == 0 && phidec == 0){fprintf(outfile, "%d\t%d\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.6e\t%1.6e\t%1.6e\n", idx, s, j, phisize - phidec * 0.1, depcentX, depcentY, sourceZ[s], depcentX - locX[j], depcentY - locY[j], sigma2, square_distance, massloading_loc_source_phi[idx], massreleased[ips]);}
		}
	}
} // End of the function

// D01b
//
// [5.5.2b]
//
// Calculate mass loading on a certain grain size on a certain point on the ground.
// All grainsizes and sources are summed up from the "elements", which is calculated by D01a.
void calc_mass_loading_location(int phiint, double *massloading_loc_source_phi, double *massloading, double *ttl, double *cummassphi){
	int idx;
	double phi;

	for(int j = 0; j < LOCDIM; j++){
		for(idx = PHIDECDIM * SDIMCUTOFF * j; idx < PHIDECDIM * SDIMCUTOFF * (j + 1); idx++){
			phi = phiint + MAX_GRAINSIZE + 1; // - phidec * INTERVAL_DECIMAL_PHI;
			massloading[j] += massloading_loc_source_phi[idx];
			ttl[j] += massloading_loc_source_phi[idx];
			cummassphi[j] += massloading_loc_source_phi[idx] * phi;
		}
	}
	//printf("LINE302 j = %d\n", j);
}

/* NON CUDA FUNCTIONS (END)*/



void set_coordinates_to_location_properties(DEP *location_properties, double *locX, double *locY, double *locZ){
  for(int j = 0; j < LOCDIM ; j++){
    location_properties[j].x = locX[j];
		location_properties[j].y = locY[j];
		location_properties[j].z = locZ[j];
		location_properties[j].dist = sqrt(pow(locX[j], 2) + pow(locY[j], 2));
  }
}

/*
 * Store mass loading for the current integer phi class.
 *
 * This function is intended to be called inside the phi loop,
 * so that mass loading is accumulated for each phi class.
 */
void store_massloading_for_phi(int size, DEP *location_properties, double *loading){
	int phi;
	phi = size + MAX_GRAINSIZE + 1;
	
	for(int j = 0; j < LOCDIM; j++){
		location_properties[j].dep[size] = loading[j];
		if(phi > 0){location_properties[j].smallerthan1mm += loading[j];}	
	}
}

/////////////////////////////////[PART 06]///////////////////////////////// 
// ==============================
// [6] Post-processing
// ==============================
// [6.1.] store_total_massloading_and_mean_phi
// [6.2.] createisopachdata
// [6.2.1.] generate_isopach_analysis


/*
 * Store total mass loading and mean grain size at each location.
 *
 * [6.1.]
 *
 * For each ground location j:
 * - assign total mass loading (ttlmassloading)
 * - compute mean grain size (in phi) as mass-weighted average:
 *     mean phi = cummassphi / ttlmassloading
 *
 * If total mass loading is zero, mean diameter is set to -9999.
 */
void store_total_massloading_and_mean_phi(DEP *location_properties, double *ttlmassloading, double *cummassphi){
	for(int j = 0; j < LOCDIM; j++){
		location_properties[j].ttlmassloading = ttlmassloading[j];
		if(ttlmassloading[j] > 0){
			location_properties[j].meandiameter = cummassphi[j] / ttlmassloading[j];
		}else{
			location_properties[j].meandiameter = -9999;
		}
	}
}

/*
 * Perform isopach analysis on deposit data.
 *
 * [6.2.]
 *
 * This function:
 * - sorts locations by total mass loading (descending)
 * - generates isopach (thickness–area) relationship
 *
 * Input:
 *   l :
 *     array of deposit data (locations)
 *
 * Notes:
 * - Sorting modifies the input array in place.
 * - The sorted order is required for cumulative area analysis.
 */
void createisopachdata(DEP *l){
	qsort(l, LOCDIM, sizeof(DEP), compare_ttlmassloading);
	generate_isopach_analysis(l);  // print out S-A relation

	//qsort(l, LOCDIM, sizeof(DEP), compare_Md);
	//countmeandiameter(l);
}

/* [6.2.1.]
 * Generate isopach analysis from deposition results.
 *
 * Prerequisite:
 * The input array l (location_properties) must be sorted in descending
 * order of total mass loading (ttlmassloading).
 *
 * This function analyzes the final mass-loading distribution and writes
 * several summary files for isopach-based interpretation.
 *
 * Analysis workflow:
 *
 * 1. Extract the distribution axis.
 *    - Select locations inside the map domain.
 *    - Traverse the sorted array and keep points that are progressively
 *      farther from the vent.
 *    - Store these axis points in l2.
 *
 * 2. Write raw S-A data to S_vs_Area.txt.
 *    - S: mass loading
 *    - A: estimated isopach area (based on grid count)
 *    - Also outputs sqrt(A), distance, coordinates, direction,
 *      mean grain size, and fine fraction.
 *
 * 3. Interpolate S-A data at powers-of-two mass-loading levels.
 *    - Results are written to S_vs_Area_summary.txt.
 *
 * 4. Compute D-F indicators.
 *    - D: area at 0.01 * MaxS
 *    - F: fine fraction at 0.1 * MaxS
 *    - Results are written to DF.txt.
 *
 * 5. Interpolate using sqrt(A) spacing.
 *    - Resample the distribution at uniform sqrt(A) intervals.
 *    - Results are written to Area_vs_S_summary.txt.
 */
void generate_isopach_analysis(DEP *l){
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

	// Count grid points on the dispersal axis.
	// For each mass-loading level, the axis is defined by the farthest
	// grid point from the vent within the map domain.
	for(j=0; j<LOCDIM; j++){
		if(l[j].x > MAPENDW && l[j].x < MAPENDE && l[j].y > MAPENDS && l[j].y < MAPENDN){
					if(l[j].dist > distance){
						distance = l[j].dist;
						jmax++;
					} 
		}
	}
	
	maxS = l[0].ttlmassloading;
	
	// Store grid points on the dispersal axis into l2.
	// For each mass-loading level in the sorted location array,
	// select the farthest grid point from the vent within the map domain.
	// At the same time, count isopach area and stores massloading (S) vs isopach area (A) relation in the S_vs_Area.txt
	//
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
				dir = compute_direction_from_vent(l2[j2].x, l2[j2].y);
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

	/*
	* Resample S-A relationship at logarithmically spaced mass-loading levels.
	*
	* Using the axis points stored in l2, interpolate the distribution
	* at mass-loading values spaced by powers of two (S = 2^n).
	*
	* For each target S:
	* - find two neighboring axis points that bracket S
	* - linearly interpolate area (A), sqrt(A), coordinates, distance,
	*   mean grain size, and fine fraction
	*
	* Results are written to S_vs_Area_summary.txt.
	*/	
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
  			dir = compute_direction_from_vent(l2[j].x, l2[j].y);
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
  

  /*
	* Generate datasets for D-F plot (Walker, 1984).
	*
	* DF.txt:
	* - MaxS      : maximum mass loading
	* - 0.1MaxS   : mass-loading level used to estimate F
	* - 0.01MaxS  : mass-loading level used to estimate D
	* - D         : interpolated isopach area at 0.01 * MaxS
	* - F         : interpolated fine-particle fraction at 0.1 * MaxS
	*
	* Area_vs_S_summary.txt:
	* - Resample S-A data at integer sqrt(A) intervals.
	* - This dataset is used to construct Walker-style plots
	*   relating isopach area and mass loading.
	*/
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
  
  for (sqrtA = sqrtAinit; sqrtA < (int)maxsqrtA; sqrtA++){
	for(j=0; j<j2; j++){
		if(l2[j].dep[1] <= (double)sqrtA && l2[j+1].dep[1] > (double)sqrtA){
			ratio = ((double)sqrtA - l2[j].dep[1]) / (l2[j+1].dep[1] - l2[j].dep[1]);
			intA = (l2[j+1].dep[0] - l2[j].dep[0]) * ratio + l2[j].dep[0];
			intS = (l2[j+1].ttlmassloading - l2[j].ttlmassloading) * ratio + l2[j].ttlmassloading;
			intx = (l2[j + 1].x - l2[j].x) * ratio + l2[j].x;
			inty = (l2[j + 1].y - l2[j].y) * ratio + l2[j].y;
			dir = compute_direction_from_vent(l2[j].x, l2[j].y);
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
			fprintf(outfile, "%1.4e\t%1.1f\t%1.1f\t%1.4f\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.4f\n", intS, intA, (double)sqrtA, intdist, intx, inty, dir, intmean, intF);
			break;
		}
	}	// end of j loop
  } 		// end of sqrtA loop
  fclose(outfile);
}

/////////////////////[END OF THE PART 06]/////////////////////


/////////////////////////////////[PART 07]///////////////////////////////// 
//
//
// ==============================
// [7] Output
// ==============================q
// Particle deposit on the ground
// [7.1.1.] write_deposit_summary					// massloading.txt
// [7.1.2.] write_total_massloading	// debug use
// [7.1.3.] write_massloading_loc_source_phi
// [7.1.4.] write_massloading_per_phidec_at_locations

// Particle release and diffusion
// [7.2.1.] write_vertical_profiles_for_phi
// [7.2.2.] write_particle_release_theoretical_vs_actual
// [7.2.3.] write_massrelease_per_source_phi
// [7.2.4.] write_plume_files						// plumetraj.txt and plumesourceposition.txt

/*
 * Write final deposit results at all locations to file.
 *
 * [7.1.1.]
 *
 * Output file:
 *   massloading.txt
 *
 * Each row corresponds to one location j and contains:
 * - coordinates (x, y, z) [m]
 * - distance from vent [m]
 * - total mass loading [kg/m^2]
 * - fraction of particles smaller than 1 mm [%]
 * - mean grain size (phi)
 * - mass loading for each integer phi class
 *
 * Notes:
 * - Coordinates are shifted by vent position (VENT_EASTING/NORTHING).
 * - Phi classes are written from coarse to fine.
 */
void write_deposit_summary(DEP *location_properties){
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
		fprintf(outfile, "%1.0f\t%1.0f\t%1.0f\t%1.1f", location_properties[j].x + VENT_EASTING, location_properties[j].y + VENT_NORTHING, location_properties[j].z, location_properties[j].dist);
		double fraction = 0.0;
		if(location_properties[j].ttlmassloading > 0){
			fraction = location_properties[j].smallerthan1mm /
					location_properties[j].ttlmassloading * 100.0;
}
		fprintf(outfile, "\t%1.4e\t%1.4f\t%1.4f", location_properties[j].ttlmassloading, fraction, location_properties[j].meandiameter);
		for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
			fprintf(outfile, "\t%1.4e", location_properties[j].dep[phiint]);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}


/* [7.1.2.]
 * Debug utility:
 * Write total mass loading per location for the current phi class.
 * (Not used in normal work flow)
 */
void write_total_massloading(int phiint, double *ttlml){
	int phi;
	char string[20];
	FILE *outfile;
	
	phi = phiint + MAX_GRAINSIZE + 1; 
	sprintf(string, "ttlml%d.txt", phi);
	outfile = fopen(string, "w");

	fprintf(outfile, "j\tttlml\n");
	for(int j = 0; j < LOCDIM; j++){
		fprintf(outfile, "%d\t%1.4e\n", j, ttlml[j]);
	}
	fclose(outfile);
}

/* [7.1.3.]
 * Write mass loading contribution for each (location, source, phidec).
 *
 * This debug/output file shows how much each source point and decimal phi
 * class contributes to mass loading at each ground location.
 */
void write_massloading_loc_source_phi(int phiint, double *massloading_loc_source_phi){
	int j, s, phidec, phi;
	char string[64];
	FILE *outfile;
	
	phi = phiint + MAX_GRAINSIZE + 1; 
	sprintf(string, "massloading_loc_source_phi%d.txt", phi);
	outfile = fopen(string, "w");

	fprintf(outfile, "j\ts\tphi\tphidec\tmassloading_loc_source_phi\n");
	for(int i = 0; i < LOCDIM * SDIM_FOR_FALL_CALC * PHIDECDIM; i++){
		phidec = i % PHIDECDIM;
		s = (i / PHIDECDIM) % SDIM_FOR_FALL_CALC;
		j = i / (PHIDECDIM * SDIM_FOR_FALL_CALC);
		fprintf(outfile, "%d\t%d\t%1.1f\t%d\t%1.4e\n", j, s, phi - phidec * 0.1, phidec, massloading_loc_source_phi[i]);
	}
	fclose(outfile);
}

/* [7.1.4.]
 * Write mass loading at each location for decimal phi classes.
 *
 * For the current integer phi interval, this function sums
 * massloading_loc_source_phi over all source points s and outputs
 * mass loading for each decimal phi bin (phidec) at each ground location.
 *
 * Input:
 * - location_properties[j]      : properties of ground location j
 * - massloading_loc_source_phi  : mass loading indexed by (location, source, phidec)
 *
 * Output:
 * - S_decimalphi_#.txt
 */
void write_massloading_per_phidec_at_locations(
    DEP *location_properties,
    int phiint,
    double *massloading_loc_source_phi
){
    int idx;
    int phi;
    char string[64];
    FILE *outfile;

    phi = phiint + MAX_GRAINSIZE + 1;

    sprintf(string, "S_decimalphi_%d.txt", phi);
    outfile = fopen(string, "w");

    fprintf(outfile, "X\tY\tZ");
    for(int phidec = 0; phidec < PHIDECDIM; phidec++){
        fprintf(outfile, "\t%1.1f", phi - phidec * INTERVAL_DECIMAL_PHI);
    }
    fprintf(outfile, "\n");

    for(int j = 0; j < LOCDIM; j++){
        double massloading_each_phidec[PHIDECDIM];

        for(int phidec = 0; phidec < PHIDECDIM; phidec++){
            massloading_each_phidec[phidec] = 0.0;
        }

        for(int s = 0; s < SDIMCUTOFF; s++){
            for(int phidec = 0; phidec < PHIDECDIM; phidec++){
                idx = ((j * SDIMCUTOFF) + s) * PHIDECDIM + phidec;
                massloading_each_phidec[phidec] += massloading_loc_source_phi[idx];
            }
        }

        fprintf(outfile, "%1.1f\t%1.1f\t%1.1f",
                location_properties[j].x + VENT_EASTING,
                location_properties[j].y + VENT_NORTHING,
                location_properties[j].z);

        for(int phidec = 0; phidec < PHIDECDIM; phidec++){
            fprintf(outfile, "\t%1.4e", massloading_each_phidec[phidec]);
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
}

/*
 * Write vertical profiles indexed by (phi, z).
 *
 * [7.2.1.]
 *
 * Outputs a table where:
 * - rows correspond to height index z
 * - columns correspond to phi classes
 *
 * Data layout:
 *   profile_by_phi[z + phi * ZDIM]
 *
 * Used for fall time, drift, and other height-dependent quantities.
 */
void write_vertical_profiles_for_phi(
    const char *filename,
    double *h,
    double *profile_by_phi,
    double base,
    int stepnum,
    double stepdelta
){
    FILE *outfile = fopen(filename, "w");

    // write header
    fprintf(outfile, "z\th(m)");
    for(int phiint = 0; phiint < stepnum; phiint++){
        fprintf(outfile, "\tphi%1.1f", base - phiint * stepdelta);
    }
    fprintf(outfile, "\n");

    // write data (top → bottom)
    for(int z = ZDIM - 1; z >= 0; z--){
        fprintf(outfile, "%d\t%1.1f", z, h[z]);
        for(int phiint = 0; phiint < stepnum; phiint++){
            fprintf(outfile, "\t%1.4f", profile_by_phi[z + phiint * ZDIM]);
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
}

/* [7.2.2]
 * Due to the upper limit of SDIM, particles that would fall beyond
 * the maximum source distance are not included in the calculation.
 *
 * This effect is more significant for fine particles, which travel farther.
 * As a result, the actual released mass can be smaller than the
 * prescribed (theoretical) release amount.
 */
void write_particle_release_theoretical_vs_actual(RELEASE *r){
	FILE *outfile;

	outfile = fopen("particle_release_theoretical_vs_actual.txt", "w");
	fprintf(outfile, "#i\tFraction(phi)\tTheoretical(kg)\tWtPercent\tActual(kg)\tActual/Theoretical\n");

	for(int phiint = MIN_GRAINSIZE - MAX_GRAINSIZE - 1; phiint >= 0; phiint--){
		fprintf(outfile, "%d\t%1.1f\t%1.4e\t%1.4f\t%1.4e\t%1.4f\n", phiint, r[phiint].phi, r[phiint].theoretical, r[phiint].theoretical/ERUPTION_MASS*100, r[phiint].actual, r[phiint].actual/r[phiint].theoretical);
	}

	fclose(outfile);
}

// [7.2.3.]
/* Write mass release per source point for each integer phi class */
void write_massrelease_per_source_phi(SEG *r){
	FILE *outfile;

	outfile = fopen("particle_released_per_ds.txt", "w");
	
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

/*
 * Write plume trajectory and particle source positions to files.
 *
 * [7.2.4.]
 *
 * plumetraj.txt:
 *   Centerline of the plume obtained by solving the plume differential equations.
 *   Each entry corresponds to a point along the plume axis.
 *
 * plumesourceposition.txt:
 *   Discrete particle release points distributed along the plume axis.
 *   These serve as sources of particle emission for the fall calculation.
 */
void write_plume_files(
    int WRITE_COLUMN_FILES,
    int SDIM_FOR_PLUME_CALC,
    int SDIM_FOR_FALL_CALC,
    double *plume_trajX,
    double *plume_trajY,
    double *plume_trajZ,
    double *plume_trajR,
    double *plume_trajT,
    double *sourceX,
    double *sourceY,
    double *sourceZ,
    double *sourceRadius,
    double *sourceT
){
    if(WRITE_COLUMN_FILES){
        FILE *outfile = fopen("plumetraj.txt", "w");
        const char *header = "calc_step\tx\ty\tz\tR\ttime\n";
        printxyzq(outfile, header, SDIM_FOR_PLUME_CALC,
                  plume_trajX, plume_trajY, plume_trajZ,
                  plume_trajR, plume_trajT);
        fclose(outfile);
    }

    if(WRITE_COLUMN_FILES){
        FILE *outfile = fopen("plumesourceposition.txt", "w");
        const char *header = "source\tx\ty\tz\tR\ttime\n";
        printxyzq(outfile, header, SDIM_FOR_FALL_CALC,
                  sourceX, sourceY, sourceZ,
                  sourceRadius, sourceT);
        fclose(outfile);
    }
}
/////////////////////[END OF THE PART 07]/////////////////////


/////////////////////////////////[PART 08]///////////////////////////////// 
//
//
// ==============================
// [8] Utilities
// ==============================
// [8.1] clear_array
// [8.2] idx_psz
// [8.3] idx_ps
// [8.4] get_line_number
// [8.5] get_wind_line_number

// --- write x-y-z-parameter table ---
// [8.6] printxyz
// [8.7] printxyzq
// [8.8] printxyze

// --- analysis ---
// [8.9.]  countmeandiameter
// [8.10.] compare_ttlmassloading
// [8.11.] compare_Md
// [8.12.] compute_direction_from_vent


/*
 * [8.1] Initialize a double array with zeros.
 */
void clear_array(int dim, double *ary){
    for(int i = 0; i < dim; i++){
        ary[i] = 0.0;
    }
}

/*
 * [8.4] Count non-comment lines in a text file.
 *
 * Lines beginning with '#' are ignored.
 */
int get_line_number(FILE *f){

    char line[1000];

    int i = 0;

    while(NULL != fgets(line, sizeof(line), f)){

        if(line[0] == '#'){
            continue;
        }

        i++;
    }

    return(i);
}

/*
 * [8.5.] Count wind-profile levels in a wind file.
 *
 * Reads vertical wind-profile data of the form:
 *
 *   height  speed  direction  temperature  pressure
 *
 * Comment lines beginning with '#' are ignored.
 *
 * If the first valid level is not located at z = 0,
 * one additional near-ground level is reserved.
 *
 * Returns:
 *   Total number of wind levels including the
 *   additional ground level when required.
 */
int get_wind_line_number(FILE *f){

    char line[1000];

    int i = 0;
    int additional = 1;

    double wind_height;
    double wind_speed;
    double wind_dir;
    double wind_temp;
    double wind_pres;

    while(NULL != fgets(line, sizeof(line), f)){

        /* Skip comments */
        if(line[0] == '#'){
            continue;
        }

        /* Parse one wind-profile line */
        int ret = sscanf(line,
                         "%lf %lf %lf %lf %lf",
                         &wind_height,
                         &wind_speed,
                         &wind_dir,
                         &wind_temp,
                         &wind_pres);

        /* Ignore malformed lines */
        if(ret != 5){
            continue;
        }

        /*
         * If the first valid level is not z = 0,
         * allocate one extra ground-level slot.
         */
        if(wind_height == 0.0 && i == 0){
            additional = 0;
        }

        i++;
    }

    return(i + additional);
}

// [8.6.]
/* 1. x, y, z */
void printxyz(FILE *in, const char *header, int imax, double *sourceX, double *sourceY, double *sourceZ){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\n", i, sourceX[i], sourceY[i], sourceZ[i]);
  }
}

// [8.7.]
void printxyzq(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q, double *t){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", i, x[i], y[i], z[i], q[i], t[i]);
  }
}

// [8.8.]
void printxyze(FILE *in, const char *header, int imax, double *x, double *y, double *z, double *q){
  fprintf(in, "%s", header);
  for(int i = 0; i < imax; i++){
    fprintf(in, "%d\t%1.4f\t%1.4f\t%1.4f\t%1.4e\n", i, x[i], y[i], z[i], q[i]);
  }
}

/* [8.9.]
 * Compute area–mean grain size relationship of deposits.
 *
 * This function analyzes deposited mass at all locations and computes:
 * - cumulative area exceeding a given mean grain size (phi)
 * - maximum distance from vent for that area
 *
 * Two output files are generated:
 *
 *   Mean_vs_Area_summary.txt :
 *     Summary for discrete phi thresholds (0.1 phi bins)
 *     Columns:
 *       mean(phi), area [km^2], max distance [m],
 *       x, y coordinates of max distance, direction from vent [deg]
 *
 *   Mean_vs_Area_all.txt :
 *     Continuous record for all locations
 *     Columns:
 *       mean(phi), cumulative area [km^2], max distance [m]
 *
 * Notes:
 * - Only locations with ttlmassloading >= MINIMUM_DEPOSIT_FOR_MD_CALC are used.
 * - Locations are assumed to be processed in order of increasing distance.
 * - Area is computed as count × mesh area (MESH_SIZE_IN_KM^2).
 * - The loop terminates when mean diameter exceeds phi > 5.
 */
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
				dir = compute_direction_from_vent(x, y);
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

/* [8.10.]
 * Comparison function for sorting DEP by total mass loading.
 *
 * Sort order:
 *   descending (largest ttlmassloading first)
 *
 * Used with qsort() in isopach analysis.
 */
int compare_ttlmassloading(const void *a, const void *b){
    double z1 = ((DEP *)a)->ttlmassloading;
    double z2 = ((DEP *)b)->ttlmassloading;

    if (z1 < z2) return 1;
    if (z1 > z2) return -1;
    return 0;
}

/* [8.11.]
 * Comparison function for sorting DEP by mean grain size (phi).
 *
 * Sort order:
 *   ascending (smaller phi first)
 */
int compare_Md(const void *a, const void *b){
    double z1 = ((DEP *)a)->meandiameter;
    double z2 = ((DEP *)b)->meandiameter;

    if (z1 > z2) return 1;
    if (z1 < z2) return -1;
    return 0;
}

/* [8.12.]
 * Compute direction from vent to (x, y) in degrees.
 *
 * Angle is measured clockwise from north (y-axis),
 * consistent with typical geographic convention.
 */
double compute_direction_from_vent(double x, double y){
    double dir;

    dir = atan2(x, y) * 180.0 / M_PI;

    if(dir < 0){
        dir += 360.0;
    }

    return dir;
}

/////////////////////[END OF THE PART 08]/////////////////////



/////////////////////////////////[PART 08]///////////////////////////////// 
//
//
// ==============================
// [9] Cleanup
// ==============================
/* Free all dynamically allocated memory used in the simulation.
 *
 * [9.1.]
 *
 * This includes:
 * - input data arrays
 * - plume and source arrays
 * - atmospheric profiles
 * - working arrays for fall, drift, and mass loading
 * - location properties and release data structures
 *
 * Centralizing deallocation helps prevent memory leaks and
 * keeps resource management consistent.
 */
void free_all(
    double *wind_alt,
    double *wind_v,
    double *wind_dir,
    double *wind_tmp,
    double *wind_pres,
    double *locX,
    double *locY,
    double *locZ,
    double *plume_trajX,
    double *plume_trajY,
    double *plume_trajZ,
    double *plume_trajR,
    double *plume_trajT,
    double *sourceX,
    double *sourceY,
    double *sourceZ,
    double *sourceRadius,
    double *sourceT,
    double *h,
    double *atmT,
    double *atmP,
    double *windX,
    double *windY,
    double *ttlfalltime,
    double *driftX,
    double *driftY,
    double *ttlfalltime_phiint,
    double *ttldriftX_phiint,
    double *ttldriftY_phiint,
    double *ttlfalltime_phidec,
    double *ttldriftX_phidec,
    double *ttldriftY_phidec,
    double *massreleased_per_ds_and_phidec,
    SEG *massreleased_per_ds,
    double *cloud_center_x,
    double *cloud_center_y,
    double *cloud_sigma2,
    double *tmpmassloading,
    double *ttlmassloading,
    double *cummassphi,
    DEP *location_properties,
    RELEASE *r
){
    /*
     * Free all dynamically allocated arrays used in the simulation.
     * Centralizing deallocation helps prevent memory leaks and
     * keeps resource management consistent.
     */

    free(wind_alt);
    free(wind_v);
    free(wind_dir);
    free(wind_tmp);
    free(wind_pres);

    free(locX);
    free(locY);
    free(locZ);

    free(plume_trajX);
    free(plume_trajY);
    free(plume_trajZ);
    free(plume_trajR);
    free(plume_trajT);

    free(sourceX);
    free(sourceY);
    free(sourceZ);
    free(sourceRadius);
    free(sourceT);

    free(h);
    free(atmT);
    free(atmP);
    free(windX);
    free(windY);

    free(ttlfalltime);
    free(driftX);
    free(driftY);

    free(ttlfalltime_phiint);
    free(ttldriftX_phiint);
    free(ttldriftY_phiint);

    free(ttlfalltime_phidec);
    free(ttldriftX_phidec);
    free(ttldriftY_phidec);

    free(massreleased_per_ds_and_phidec);
    free(massreleased_per_ds);

    free(cloud_center_x);
    free(cloud_center_y);
    free(cloud_sigma2);

    free(tmpmassloading);
    free(ttlmassloading);
    free(cummassphi);

    free(location_properties);
    free(r);
}

/////////////////////[END OF THE PART 09]/////////////////////
