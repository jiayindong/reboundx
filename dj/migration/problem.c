/**
Reproduce the study of GJ876 by Lee & Peale 2002
**/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    //sim->integrator            = REB_INTEGRATOR_IAS15;
    //sim->ri_ias15.epsilon      = 0.;            // accuracy parameter
    sim->integrator            = REB_INTEGRATOR_WHFAST;
    sim->dt = 1.e-2*2*M_PI;

    sim->heartbeat      = heartbeat;

    struct reb_particle star;
    star.m      = 0.32;
    reb_add(sim, star);

    // struct reb_particle primary = sim->particles[0];
    // double m1 = 1.e-3;
    // double m2 = 1.e-3;
    // double a1 = 0.5;
    // double a2 = 1.;
    // double e = 0.;
    // double inc = 0.;
    // double Omega = 0.;
    // double omega = 0.;
    // double f = 0.;

    //struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, primary, m1, a1, e, inc, Omega, omega, f);
    //struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, primary, m2, a2, e, inc, Omega, omega, f);

    struct reb_particle p1 = {0};    // Planet 1
    p1.x      = 0.5;
    p1.m      = 0.56e-3;
    p1.vy     = sqrt(sim->G*(star.m+p1.m)/p1.x);
    reb_add(sim, p1); 
    
    struct reb_particle p2 = {0};    // Planet 2
    p2.x      = 1.;
    p2.m      = 1.89e-3;
    p2.vy     = sqrt(sim->G*(star.m+p2.m)/p2.x);
    reb_add(sim, p2); 

    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_attach(sim);

    // There are two options for how to modify orbits.  You would only choose one (comment the other out).
    // You can't set precession separately with modify_orbits_forces (eccentricity and inclination damping induce pericenter and nodal precession).

    struct rebx_operator* mo = rebx_load_operator(rebx, "modify_orbits_direct");                                     // directly update particles' orbital elements each timestep
    rebx_add_operator(rebx, mo);
    //struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");                                         // add forces that orbit-average to give exponential a and e damping
    //rebx_add_force(rebx, mo);

    // Set the timescales for each particle.
    tmax = 2.e4*2.*M_PI;

    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", -tmax);         // add semimajor axis damping on outer planet (e-folding timescale)
    //rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_omega", -tmax/100.); // add linear precession (set precession period). Won't do anything for modify_orbits_forces
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_e", -tmax/100.);     // add eccentricity damping on outer planet(e-folding timescale)

    /* One can also adjust a coupling parameter between eccentricity and semimajor axis damping.  We use the parameter p
     * as defined by Deck & Batygin (2015).  The default p=0 corresponds to no coupling, while p=1 corresponds to e-damping
     * at constant angular momentum.  This is only implemented for modify_orbits_direct.
     * modify_orbits_forces damps eccentricity at constant angular momentum.
     *
     * Additionally, the damping by default is done in Jacobi coordinates.  If you'd prefer to use barycentric
     * coordinates, or coordinates referenced to a particular particle, set a coordinates parameter in the effect
     * parameters returned by rebx_add to REBX_COORDINATES_BARYCENTRIC or REBX_COORDINATES_PARTICLE.
     * If the latter, add a 'primary' flag to the reference particle (not neccesary for barycentric):
     */

    rebx_set_param_double(rebx, &mo->ap, "p", 1); // e-damping at constant angular momentum
    //rebx_set_param_int(rebx, &mo->ap, "coordinates", REBX_COORDINATES_PARTICLE);
    //rebx_set_param_int(rebx, &sim->particles[0].ap, "primary", 1);

    system("rm -v orbits_direct.txt");

    reb_integrate(sim, tmax);
    rebx_free(rebx);    // Free all the memory allocated by rebx
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 20.*M_PI)){
        reb_output_timing(sim, tmax);
    }
    if(reb_output_check(sim, 40.)){
        reb_output_orbits(sim,"orbits_direct.txt");
    }
}