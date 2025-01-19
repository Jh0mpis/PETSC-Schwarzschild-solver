/* 
  Schwarzschild geodesics solver:

  The following program solves the differential equations given by the Schwarzschild metric in General Relativity (GR).

  The differential equations are the following ones:

  t' = E r / (r - 2 M),
  r'' = -k M / r^2 + L_z^2 / r^3 - 3 L_z M / r^4,
  phi' = L_z / r^2.

  In order to solve the differential equation as a first order differential equation, we use the following substitution:

  r' = R,

  then, the system to solve is the following one:

  t' = E r / (r - 2 M),
  r' = R,
  R' = - k M / r^2 + L_z^2 / r^3 - 3 L_z M / r^4,
  phi' = L_z / r^2.

  Using the following solution vector:

  t
  r
  R
  phi
(2.0 * system->M * system->k * r*r - 3.0 * system->L_z * system->L_z * r + 12.0 * system->L_z * system->M) / (r*r*r*r*r)

  The Jacobian writes as follow:

  0            -2 M / (r - 2 M)^2             0 0
  0                    0                      1 0
  0 (2 M k r^2 - 3 L_z^2 r + 12 L_z M) / r^5  0 0
  0               -2 L_z / r^3                0 0

  We are going to solve the system

  X' = F(X,t)

*/
#include <petsc.h>
#include <stdlib.h>

// Structure for save the system parameters
typedef struct {
  PetscReal E;      // Particle energy
  PetscReal L_z;    // Angular momentum
  PetscReal M;      // Black hole mass
  PetscReal k;      // Particle mass parameter
} Parameters;

// Function for define the right hand side term F(X,t).
// ts: TS handler.
// t: current time (not the proper time)
// curr_x: Current solution.
// diff_out: Right hand side output.
// ctx: For additional data, like the system struct.
PetscErrorCode diff_func(TS ts, PetscReal t, Vec curr_x, Vec diff_out, void *ctx){
  const PetscScalar *curr_x_ptr;
  Parameters *system = (Parameters *) ctx;
  PetscReal r, R;

  PetscFunctionBeginUser;

  // Get the vec data pointer
  PetscCall(VecGetArrayRead(curr_x, &curr_x_ptr));

  // Defining the useful coordinates
  r = curr_x_ptr[1];
  R = curr_x_ptr[2];
  PetscCall(VecRestoreArrayRead(curr_x, &curr_x_ptr));

  PetscCall(VecSetValue(diff_out, 0,
                        system->E * r / (r - 2.0 * system->M), INSERT_VALUES));   // Differential equation for time
  PetscCall(VecSetValue(diff_out, 1, R, INSERT_VALUES));                          // Differential equation for the radius
  PetscCall(VecSetValue(diff_out, 2, - system->M * system->k / (r*r) + system->L_z*system->L_z / (r*r*r) - 3.0 * system->L_z * system->M / (r*r*r*r) , INSERT_VALUES));   // Differential equation for the radial speed
  PetscCall(VecSetValue(diff_out, 3, system->L_z / (r*r), INSERT_VALUES));        // Differential equation for the angle

  PetscCall(VecAssemblyBegin(diff_out));
  PetscCall(VecAssemblyEnd(diff_out));


  PetscFunctionReturn(PETSC_SUCCESS);
}

// Function for defining the right hand side Jacobian.
// ts : TS handler
// t : current time (not the proper time)
// curr_x : Current solution.
// diff_J : Jacobian
// diff_P : Pre-conditioner 
// ctx: For additional data, like the system struct.
PetscErrorCode diff_Jacobian(TS ts, PetscReal t, Vec curr_x, Mat diff_J, Mat diff_P, void *ctx){

  const PetscScalar *curr_x_ptr;
  Parameters *system = (Parameters *) ctx;
  PetscReal r;

  PetscFunctionBeginUser;

  // Get the vec data pointer
  PetscCall(VecGetArrayRead(curr_x, &curr_x_ptr));

  // Defining the useful coordinates
  r = curr_x_ptr[1];

  PetscCall(VecRestoreArrayRead(curr_x, &curr_x_ptr));

  // Computing the Jacobian Matrix
  PetscCall(MatSetValue(diff_J, 0, 1, -2.0 * system->E * system->M/((r - 2.0 * system->M)*(r - 2.0 * system->M)), INSERT_VALUES));
  PetscCall(MatSetValue(diff_J, 1, 2, 1.0, INSERT_VALUES));
  PetscCall(MatSetValue(diff_J, 2, 1, (2.0 * system->M * system->k * r*r - 3.0 * system->L_z * system->L_z * r + 12.0 * system->L_z * system->M) / (r*r*r*r*r), INSERT_VALUES));
  PetscCall(MatSetValue(diff_J, 3, 1, -2.0 * system->L_z / (r*r*r), INSERT_VALUES));

  PetscCall(MatAssemblyBegin(diff_J, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(diff_J, MAT_FINAL_ASSEMBLY));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// Function for return the solution in polar coordinates
PetscErrorCode polar_monitor(TS ts, PetscInt stepid, PetscReal tau, Vec curr_x, void *ctx){
  const PetscReal *curr_x_ptr;
  Parameters* system = (Parameters *) ctx;
  PetscReal t, r, R, phi;

  PetscFunctionBeginUser;

  // Get the vec data pointer
  PetscCall(VecGetArrayRead(curr_x, &curr_x_ptr));

  // Defining the useful coordinates
  t = curr_x_ptr[0];
  r = curr_x_ptr[1];
  R = curr_x_ptr[2];
  phi = curr_x_ptr[3];

  PetscCall(VecRestoreArrayRead(curr_x, &curr_x_ptr));

  // Print the solution
  if(r > 2 * system->M){
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "%g\t%g\t%g\t%g\t%g\n", t, r, R, phi, 0.5 * R * R + system->k * system->M / r + system->L_z * system->L_z/(2*r*r) - system->L_z * system->M/(r*r*r)));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// Function for return the solution in cartesian coordinates
PetscErrorCode cartesian_monitor(TS ts, PetscInt stepid, PetscReal tau, Vec curr_x, void *ctx){
  const PetscReal *curr_x_ptr;
  Parameters* system = (Parameters *) ctx;
  PetscReal t, r, R, phi;

  PetscFunctionBeginUser;

  // Get the vec data pointer
  PetscCall(VecGetArrayRead(curr_x, &curr_x_ptr));

  // Defining the useful coordinates
  t = curr_x_ptr[0];
  r = curr_x_ptr[1];
  R = curr_x_ptr[2];
  phi = curr_x_ptr[3];

  PetscCall(VecRestoreArrayRead(curr_x, &curr_x_ptr));

  // Print the solution
  if(r > 2 * system->M){
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "%g\t%g\t%g\t%g\t%g\n", t, r * PetscCosReal(phi), r * PetscSinReal(phi), R, 0.5 * R * R + system->k * system->M / r + system->L_z * system->L_z/(2*r*r) - system->L_z * system->M/(r*r*r)));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// main function
int main(int argc, char *argv[])
{
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  Parameters system;
  PetscReal M = 1.0, E = 0.1, L_z = 3.0;
  PetscReal r_0, R_0 = 0.0, phi_0 = 0.0;
  PetscInt total_time = 100, solver = 0;
  PetscBool k = 1, monitor = 0;
  Mat J;
  Vec x;
  TS ts;

  // ------- Receiving the console parameters ------ //
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-monitor", &monitor, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-total_time", &total_time, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-solver", &solver, NULL));


  // ------- Receiving the system parameters ------ //
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-M", &M, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-E", &E, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-L_z", &L_z, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-k", &k, NULL));


  r_0 = 20.0 * M;
  // ------- Receiving the initial conditions ------ //
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-r_0", &r_0, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-vel_0", &R_0, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-phi_0", &phi_0, NULL));

  // ------- Setting up the system ------ //
  system.E = E;
  system.M = M;
  system.L_z = L_z;

  if (k){
    system.k = 1.0;
  }else{
    system.k = 0.0;
  }

  PetscCall(PetscPrintf(PETSC_COMM_SELF, "Solving the system with: \n\tM \t=\t %g\n\tE \t=\t %g\n\tL_z \t=\t %g\n\tk \t=\t %g\n\n",
                        system.M, system.E, system.L_z, system.k));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "Initial conditions: \n\tr(0) \t=\t %g\n\tr'(0) \t=\t %g\n\tphi(0) \t=\t %g\n\n",
                        r_0, R_0, phi_0));


  // ------- Creating the data for solve the diff equation ------ //
  PetscCall(VecCreate(PETSC_COMM_SELF, &x));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecSetSizes(x, PETSC_DECIDE, 4));
  PetscCall(VecSetUp(x));

  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, 4, 4, NULL, &J));

  // ------- Setting up the initial conditions ------ //
  PetscCall(VecSetValue(x, 1, r_0, INSERT_VALUES));
  PetscCall(VecSetValue(x, 2, R_0, INSERT_VALUES));
  PetscCall(VecSetValue(x, 3, phi_0, INSERT_VALUES));
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));


  // ------- Create the TS handler ------ //
  PetscCall(TSCreate(PETSC_COMM_SELF, &ts));

  // ------- Setting the solver ------ //
  // Set the RHS function
  PetscCall(TSSetRHSFunction(ts, NULL, diff_func, &system));
  // Set the RHSJacobian
  PetscCall(TSSetRHSJacobian(ts, J, J, diff_Jacobian, &system));
  // Set the monitor
  if( monitor == 1 ){
    PetscCall(TSMonitorSet(ts, polar_monitor, &system, NULL));
  }else{
    PetscCall(TSMonitorSet(ts, cartesian_monitor, &system, NULL));
  }

  // ------- Choosing the solver ------ //
  if( solver == 1){
    PetscCall(TSSetType(ts, TSEULER));
  }else if (solver == 2){
    PetscCall(TSSetType(ts, TSSSP));
  }else{
    PetscCall(TSSetType(ts, TSRK));
    PetscCall(TSRKSetType(ts, TSRK4));
  }

  // ------- Setting the total time for the solver ------ //
  PetscCall(TSSetMaxTime(ts, total_time));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER));

  PetscCall(TSSetFromOptions(ts));

  // ------- Solve ------ //
  PetscCall(TSSolve(ts, x));

  
  // ------- Free ------ //
  PetscCall(MatDestroy(&J));
  PetscCall(VecDestroy(&x));
  PetscCall(TSDestroy(&ts));
  
  PetscCall(PetscFinalize());
  return 0;
}
