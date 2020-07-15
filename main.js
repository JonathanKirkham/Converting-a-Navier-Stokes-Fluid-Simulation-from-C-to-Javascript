const main = () => {

    let verbose = 1;          /* Verbosity level */
    let xlength = 22.0;     /* Width of simulated domain */
    let ylength = 4.1;      /* Height of simulated domain */
    let imax = 330;           /* Number of cells horizontally */
    let jmax = 60;           /* Number of cells vertically */

    let outname;
    let output = 0;
    let output_frequency = 0;
    let t_end = 40; //2.1       /* Simulation runtime */

    let del_t = 0.003;      /* Duration of each timestep */
    let tau = 0.5;          /* Safety factor for timestep control */

    let itermax = 10;        /* Maximum number of iterations in SOR */
    let eps = 0.001;        /* Stopping error threshold for SOR */
    let omega = 1.7;        /* Relaxation parameter for SOR */
    let gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    let Re = 150.0;         /* Reynolds number */
    let ui = 1.0;           /* Initial X velocity */
    let vi = 0.13;           /* Initial Y velocity */


    let t, delx, dely;
    let i, j, itersor = 0, ifluid = 0, ibound = 0;
    let res;

    let flag;
    let init_case, iters = 0;
    let show_help = 0, show_usage = 0, show_version = 0;


    delx = xlength / imax;
    dely = ylength / jmax;

    /* Allocate arrays */
    u = alloc_matrix(imax + 2, jmax + 2);
    v = alloc_matrix(imax + 2, jmax + 2);
    f = alloc_matrix(imax + 2, jmax + 2);
    g = alloc_matrix(imax + 2, jmax + 2);
    p = alloc_matrix(imax + 2, jmax + 2);
    rhs = alloc_matrix(imax + 2, jmax + 2);
    flag = alloc_matrix(imax + 2, jmax + 2);




    let checker = 0;
    let checker1 = 0.0;

    // Set up initial values
    for (i = 0; i <= imax + 1; i++) {
        for (j = 0; j <= jmax + 1; j++) {
            checker += (i * jmax) + j + 1;
            checker1 += (i * jmax) + j + 1.0;
            u[i][j] = ui;
            v[i][j] = vi;
            p[i][j] = 0.0;
        }
    }

    // init_flag sets the domain's boundary shape and apply_boundary_conditions
    // ensures u and v adhere to this shape 
    ibound = init_flag(flag, imax, jmax, delx, dely, ibound);

    // Enforces the velocity conditions at the boundaries
    apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);




    // Main loop is called asynchronously and computes the fluid velocity and draws to screen
    // The function establishes a loop by calling setTimeout.
    var mainLoop = function (t, iters) {

        // This function computes the maximum timestep for del_t such that it
        //satisfies the CFL condition, which avoids instability in the simulation	
        del_t = set_timestep_interval(del_t, imax, jmax, delx, dely, u, v, Re, tau);

        // Calculate the number of empty cells (non-boundary cells)
        ifluid = (imax * jmax) - ibound;

        // Compute initial guess for the velocity (f and g)
        compute_tentative_velocity(u, v, f, g, flag, imax, jmax,
            del_t, delx, dely, gamma, Re);

        compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);

        // Calculates the pressure in each cell from the rhs
        if (ifluid > 0) {
            Arr = poisson(p, rhs, flag, imax, jmax, delx, dely,
                eps, itermax, omega, res, ifluid);
            itersor = Arr[0];
            res = Arr[1];
        } else {
            itersor = 0;


        }

        console.log(iters, t + del_t, del_t, itersor, res, ibound);

        // Calculates the velocities (u and v) from the tentative velocities and the pressures
        update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);

        //Makes sure that the fluid is interacting with the boundary cells
        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);


        write_ppm(u, v, p, flag, imax, jmax, xlength, ylength, outname,
            iters, output_frequency, document.querySelector('canvas'));


        //  Termination Condition, set loop is time is not exhausted
        if (t < t_end)
            setTimeout(function () { mainLoop(t + del_t, iters + 1); }, 0);

    }

    // Asynchronous Initialisation
    setTimeout(function () {
        set_output(flag, imax, jmax, document.querySelector('canvas'));
        mainLoop(0.0, iters);
    }, 100);
};
main();