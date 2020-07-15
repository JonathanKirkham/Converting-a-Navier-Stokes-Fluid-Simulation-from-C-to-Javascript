/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
function init_flag(flag, imax, jmax, delx, dely, ibound) {
    let i, j;
    let mx, my, x, y, rad1;

    /* Mask of a circular obstacle */
    mx = 20.0 / 41.0 * jmax * dely;
    my = mx;
    rad1 = 5.0 / 41.0 * jmax * dely;
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            x = (i - 0.5) * delx - mx;
            y = (j - 0.5) * dely - my;
            flag[i][j] = (x * x + y * y <= rad1 * rad1) ? C_B : C_F;
        }
    }
    // To return the project back to it's original state, comment out lines 24-48	
    //  This code creates a canvas and draws shapes to it.
    //  The canvas is then used to define the obstacle cells
    let canvas = document.createElement('canvas');
    // canvas.width = imax;
    // canvas.height = jmax;
    // let context = canvas.getContext('2d');

    // context.fillStyle = "rgb(255,255,255)";
    // BEGIN drawing obstacles to context
    //context.fillRect(0.2*imax, 0.3*jmax, 0.1*imax, 0.1*jmax);
    // Drawing a circle is more complex (drawing an ellipse is even more so)
    // context.beginPath();
    // context.arc(20, 8, 1, 0, Math.PI * 2); // The first two parameters define the centre, the third defines the radius
    // context.arc(20, 28, 2, 0, Math.PI * 2);
    // context.arc(20, 48, 3, 0, Math.PI * 2);
    // context.fill();

    // // END drawing obstacles to context

    // var imgData = context.getImageData(0, 0, imax, jmax);
    // for (i = 0; i < imax; i++) {
    //     for (j = 0; j < jmax; j++) {
    //         flag[i][j] = (imgData.data[4 * (j * imax + i)] == 0) ? C_F : C_B
    //     }
    // }

    // End of canvas code */

    /* Mark the north & south boundary cells */
    for (i = 0; i <= imax + 1; i++) {
        flag[i][0] = C_B;
        flag[i][jmax + 1] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (j = 1; j <= jmax; j++) {
        flag[0][j] = C_B;
        flag[imax + 1][j] = C_B;
    }

    /* flags for boundary cells */
    var ibound = 0;
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            if (!(flag[i][j] & C_F)) {
                (ibound)++;
                if (flag[i - 1][j] & C_F) flag[i][j] |= B_W;
                if (flag[i + 1][j] & C_F) flag[i][j] |= B_E;
                if (flag[i][j - 1] & C_F) flag[i][j] |= B_S;
                if (flag[i][j + 1] & C_F) flag[i][j] |= B_N;
            }

        }
    }
    return ibound;
}