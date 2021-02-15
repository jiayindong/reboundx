/**
 * Planetary embryos embedded in a protoplanetary disk with an outer Jupiter migrating inward
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

#define EMBRYO_MAX 1000

void heartbeat(struct reb_simulation* sim);
double tmax;

int main(int argc, char* argv[]){

    /* import planetary embryos from Dawson+16 */
    double mass[EMBRYO_MAX]; // Mass in solar
	double a[EMBRYO_MAX]; // Semi-major axis
	double e[EMBRYO_MAX]; // Eccentricity
	double inc[EMBRYO_MAX]; // Inc (deg)
	double omega[EMBRYO_MAX]; // Argument of peri (deg)
	double Omega[EMBRYO_MAX]; // Longitude of ascending node (Deg)
	double M[EMBRYO_MAX]; // Mean anomaly (deg)

	FILE* file = fopen("/Users/dj/Desktop/migration/embryos.txt", "r");
    if (file == NULL){
        printf("Error: Could not open embryos.txt!");
        exit(1);
    }

    memset(mass, 0, sizeof(mass));
    memset(a, 0, sizeof(a));
    memset(e, 0, sizeof(e));
    memset(inc, 0, sizeof(inc));
    memset(omega, 0, sizeof(omega));
    memset(Omega, 0, sizeof(Omega));
    memset(M, 0, sizeof(M));

    int n_embryos = 0;
    while (!feof(file) && (n_embryos < EMBRYO_MAX)){
    	fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", &(mass[n_embryos]), &(a[n_embryos]), &(e[n_embryos]),
    		&(inc[n_embryos]), &(omega[n_embryos]), &(Omega[n_embryos]), &(M[n_embryos]));
    	n_embryos++;
    }
	printf("Number of planetary embryos: %d\n", n_embryos);
    fclose(file);

    /* start the rebound simulation here */
    struct reb_simulation* sim = reb_create_simulation();
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->ri_ias15.epsilon = 0.;
    //sim->dt = 0.001*2.*M_PI;

    sim->collision            = REB_COLLISION_DIRECT;
    sim->collision_resolve    = reb_collision_resolve_merge;       // Choose merger collision routine.

    sim->heartbeat = heartbeat;

    // add the host star
    struct reb_particle star = {0};
    star.m      = 1.;
    star.r      = 0.005;
    reb_add(sim, star);

    // add the Jupiter
    struct reb_particle planet = {0};
    planet.m    = 1.e-3;
    planet.r    = 0.0005;
    planet.x    = 2.;
    planet.vy   = sqrt(sim->G*(star.m+planet.m)/planet.x);
    reb_add(sim, planet);

    // add planetary embryos
    int last = 38;
    for (int i = last; i < n_embryos; i++){
    	double f = reb_tools_M_to_f(e[i], M[i]/(2.*M_PI));
    	struct reb_particle p = reb_tools_orbit_to_particle(sim->G, star, mass[i], a[i], e[i], 
    		inc[i]/(2.*M_PI), Omega[i]/(2.*M_PI), omega[i]/(2.*M_PI), f);
    	p.r = 8.5e-5; // 2 earth radii in au
    	p.lastcollision = 0;
    	reb_add(sim, p);
    }

    reb_move_to_com(sim);
    printf("Total number of particles: %d\n", sim->N);

    tmax = 2.e4*2*M_PI;

	struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_operator* mo = rebx_load_operator(rebx, "modify_orbits_direct");
    rebx_add_operator(rebx, mo);
    rebx_set_param_double(rebx, &mo->ap, "p", 1); // e-damping at constant angular momentum

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -tmax);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_e", -tmax/100.);

    for (int i = 0; i < (n_embryos-last); i++){
    	rebx_set_param_double(rebx, &sim->particles[1+i].ap, "tau_e", -tmax/100.);
    }

    system("rm -v orbits.txt");
    system("rm -v collisions.txt");

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);

    return 0;
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 20.*M_PI)){
        reb_output_timing(sim, tmax);
    }
    if(reb_output_check(sim, 40.)){
        reb_output_orbits(sim,"orbits.txt");
    }
}