function compute_tentative_velocity(u, v, f, g, flag, imax,
				    jmax, del_t, delx, dely, gamma, Re) {

    var  i, j;
    var du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {

                du2dx = ((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])+
			 gamma*Math.abs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])-
                    (u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j])-
			 gamma*Math.abs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j]))
                    /(4.0*delx);
                duvdy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])+
                    gamma*Math.abs(v[i][j]+v[i+1][j])*(u[i][j]-u[i][j+1])-
                    (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j])-
                    gamma*Math.abs(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]-u[i][j]))
                    /(4.0*dely);
                laplu = (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/delx/delx+
                    (u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dely/dely;
   
                f[i][j] = u[i][j]+del_t*(laplu/Re-du2dx-duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }

    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {

	      duvdx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])+
                    gamma*Math.abs(u[i][j]+u[i][j+1])*(v[i][j]-v[i+1][j])-
                    (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j])-
		    gamma*Math.abs(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]-v[i][j]))
                    /(4.0*delx);
                dv2dy = ((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])+
                    gamma*Math.abs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])-
                    (v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j])-
			 gamma*Math.abs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j]))
                    /(4.0*dely);

                laplv = (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/delx/delx+
                    (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dely/dely;

                g[i][j] = v[i][j]+del_t*(laplv/Re-duvdx-dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (j=1; j<=jmax; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    for (i=1; i<=imax; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }
}


/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
function set_timestep_interval(del_t, imax, jmax,  delx,
    dely, u, v, Re, tau)

{
	var del_t;
    var i, j;
    var umax, vmax, deltu, deltv, deltRe; 



    /* del_t satisfying CFL conditions */
	if (tau >= 1.0e-10) { /* else no time stepsize control */
        umax = 1.0e-10;
        vmax = 1.0e-10; 
        for (i=0; i<=imax+1; i++) {
            for (j=1; j<=jmax+1; j++) {
                umax = Math.max(Math.abs(u[i][j]), umax);

            }
        }
        for (i=1; i<=imax+1; i++) {
            for (j=0; j<=jmax+1; j++) {
                vmax = Math.max(Math.abs(v[i][j]), vmax);

            }

        }


        deltu = delx/umax;
        deltv = dely/vmax; 
        deltRe = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.0;


        if (deltu<deltv) {
            del_t = Math.min(deltu, deltRe);
        } else {
            del_t = Math.min(deltv, deltRe);
        }
        del_t = tau * (del_t); /* multiply by safety factor */

    }
	return del_t;

}
/* Calculate the right hand side of the pressure equation * */
function compute_rhs(f, g, rhs, flag,  imax,
    jmax, del_t,  delx,  dely)
{
    var i, j;

    for (i=1;i<=imax;i++) {
        for (j=1;j<=jmax;j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = (
                             (f[i][j]-f[i-1][j])/delx +
                             (g[i][j]-g[i][j-1])/dely
                            ) / del_t;
            }
        }
    }
}

/*Red/Black SOR to solve the poisson equation */
function poisson(p, rhs, flag, imax, jmax,
    delx, dely, eps, itermax, omega,
    res, ifull)
{
    var i, j, iter;
    var add, beta_2, beta_mod;
    var p0 = 0.0;
    
    var rb; /* Red-black value. */

    var rdx2 = 1.0/(delx*delx);
    var rdy2 = 1.0/(dely*dely);
    beta_2 = -omega/(2.0*(rdx2+rdy2));

    /* Calculate sum of squares */
    for (i = 1; i <= imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j]*p[i][j]; }
        }
    }
   
    p0 = Math.sqrt(p0/ifull);
    if (p0 < 0.0001) { p0 = 1.0; }


    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (rb = 0; rb <= 1; rb++) {
            for (i = 1; i <= imax; i++) {
                for (j = 1; j <= jmax; j++) {
                    if ((i+j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.-omega)*p[i][j] - 
                              beta_2*(
                                    (p[i+1][j]+p[i-1][j])*rdx2
                                  + (p[i][j+1]+p[i][j-1])*rdy2
                                  -  rhs[i][j]
                              );
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */                 
							var eps_E = (flag[i+1][j] & C_F)?1:0;
							var eps_W = (flag[i-1][j] & C_F)?1:0;
							var eps_N = (flag[i][j+1] & C_F)?1:0;
							var eps_S = (flag[i][j-1] & C_F)?1:0;
							  beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2); 
							  p[i][j] = (1.-omega)*p[i][j] -
                            beta_mod*(
                                  (eps_E*p[i+1][j]+eps_W*p[i-1][j])*rdx2
                                + (eps_N*p[i][j+1]+eps_S*p[i][j-1])*rdy2
                                - rhs[i][j]
                            );
                    }
                } /* end of j */
            } /* end of i */
        } /* end of rb */
        
        /* Partial computation of residual */
        res = 0.0;
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {
                    /* only fluid cells */
							var eps_E = (flag[i+1][j] & C_F)?1:0;
							var eps_W = (flag[i-1][j] & C_F)?1:0;
							var eps_N = (flag[i][j+1] & C_F)?1:0;
							var eps_S = (flag[i][j-1] & C_F)?1:0;
		  add = (eps_E*(p[i+1][j]-p[i][j]) - 
			 eps_W*(p[i][j]-p[i-1][j])) * rdx2  +
		    (eps_N*(p[i][j+1]-p[i][j]) -
		     eps_S*(p[i][j]-p[i][j-1])) * rdy2  -  rhs[i][j];
                    res += add*add;
                }
            }
     
        res = Math.sqrt((res)/ifull)/p0;

        /* convergence? */
        if (res<eps) break;
    } 
	/* end of iter */
    return[iter,res];
}

	/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
	
	function update_velocity(u, v, f, g, p,
    flag, imax, jmax, del_t, delx, dely)
{
    var i, j;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j]-(p[i+1][j]-p[i][j])*del_t/delx;
            }
        }
    }
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
	      v[i][j] = g[i][j]-(p[i][j+1]-p[i][j])*del_t/dely;
            }
        }
    }
}