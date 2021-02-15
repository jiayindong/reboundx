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

    /* start the rebound simulation here */
    struct reb_simulation* sim = reb_create_simulation();
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->ri_ias15.epsilon = 0.;
    //sim->dt = 0.001*2.*M_PI;

    sim->collision            = REB_COLLISION_DIRECT;
    sim->collision_resolve    = reb_collision_resolve_merge;       // Choose merger collision routine.

    sim->heartbeat = heartbeat;

    sim->usleep = 100;

    // add the host star
    struct reb_particle star = {0};
    star.m      = 1.;
    star.r      = 0.005;
    reb_add(sim, star);

    // add the Jupiter
    struct reb_particle planet = {0};
    planet.m    = 1.e-3;
    planet.r    = 0.0005;
    planet.x    = 3.;
    planet.vy   = sqrt(sim->G*(star.m+planet.m)/planet.x);
    reb_add(sim, planet);

    // Add planets
    int N_planets = 7;
    for (int i=0;i<N_planets;i++){
        double a = 1.+(double)i/(double)(N_planets-1);        // semi major axis in AU
        printf("%f\n", a);
        double v = sqrt(1./a);                     // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4; 
        planet.r = 4e-2;                     // radius in AU (it is unphysically large in this example)
        planet.lastcollision = 0;                // The first time particles can collide with each other
        planet.x = a; 
        planet.vy = v;
        reb_add(sim, planet); 
    }
    reb_move_to_com(sim);

    printf("Total number of particles: %d\n", sim->N);

    tmax = 2.e4*2*M_PI;

	struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_operator* mo = rebx_load_operator(rebx, "modify_orbits_direct");
    rebx_add_operator(rebx, mo);
    rebx_set_param_double(rebx, &mo->ap, "p", 1); // e-damping at constant angular momentum

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -tmax);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_e", -tmax/1000.);

    system("rm -v orbits.txt");

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