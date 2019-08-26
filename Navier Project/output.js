// Sets a square rectangle of image data to the colours r,g,b
var setPixleRect = function(data,x,y,rowlength, r,g,b) {
 var i,j;
 for (i = 0; i < pixelSize; i++) {
     for (j = 0; j < pixelSize; j++) {
  var pixelIdx = (y*pixelSize + j)*rowlength + (x*pixelSize + i)
  data[4*pixelIdx + 0] = r;
  data[4*pixelIdx + 1] = g;
  data[4*pixelIdx + 2] = b;
  data[4*pixelIdx + 3] = 255;
     }
 }
    };
var m_canvas, m_context;
var outmode = 0;
var psi;
var zeta;
var pixelSize = 3;
	
//  Initialises the buffer canvas, allocates the zeta array and draws the boundary obstacle at startup
function set_output(flag, imax, jmax, cvs){
	zeta = alloc_matrix(imax+2, jmax+2);
	m_canvas = document.createElement('canvas');
	m_canvas.width = cvs.width;
	m_canvas.height = cvs.height;
	m_context = m_canvas.getContext('2d');
	m_context.fillStyle = "rgb(0,255,0)";
  for (j = 1; j < jmax+1 ; j++) {
       for (i = 1; i < imax+1 ; i++) {
         
           if (!(flag[i][j] & C_F)) {
       m_context.fillRect( i*pixelSize, j*pixelSize, pixelSize, pixelSize );
       }
	   
   }
  }
   cvs.getContext('2d').drawImage(m_canvas,0,0);
   
	
}

//  draws the vorticity to the buffer canvas, and draws the buffer canvas to the screen
function write_ppm(u, v, p, flag,
	     imax, jmax, xlength, ylength, 
	       outname, iters, freq, cvs) {
  var i, j;
  var zmax = -1e10, zmin = 1e10;
  var pmax = -1e10, pmin = 1e10;
  var umax = -1e10, umin = 1e10;
  var vmax = -1e10, vmin = 1e10;
  

 
  var delx = xlength/imax;
  var dely = ylength/jmax;
  

//  Computes the vorticity from the velocities   
calc_psi_zeta(u, v, psi, zeta, flag, imax, jmax, delx, dely);



// Gets Image Data From The Canvas and set its colour bytes to vorticity value
  var imageData = m_context.getImageData(0,0,imax*pixelSize,jmax*pixelSize);
    for (j = 1; j < jmax+1 ; j++) {
        for (i = 1; i < imax+1 ; i++) {
  var r, g, b;
     if ((flag[i][j] & C_F)) {
  var z = (i < imax && j < jmax)?zeta[i][j]:0.0;
  r = g = b = Math.round(Math.pow(Math.abs(z/12.6),.4) * 255);
  setPixleRect(imageData.data, i - 1,j - 1, imax*pixelSize, r,g,b);
		}
	}
}
    m_context.putImageData(imageData,0,0,0,0,imax*pixelSize,jmax*pixelSize);
	
	
// Uses context functions to draw rectangles on the canvas (less efficient than editing colour bytes directly)	
  /* for (j = 1; j < jmax+1 ; j++) {
       for (i = 1; i < imax+1 ; i++) {
           var r, g, b;
           if ((flag[i][j] & C_F)) {
		
                  var z = zeta[i][j];
                 r = g = b = Math.round(Math.pow(Math.abs(z/12.6),.4) * 255);
				 m_context.fillStyle = "rgb("+r+","+g+","+b+")";
					m_context.fillRect( i*pixelSize, j*pixelSize, pixelSize, pixelSize );
      }
	   
       }
   } */
   // Copy Buffer to the main canvas
   cvs.getContext('2d').drawImage(m_canvas,0,0);
   
 }
 

 


/* Computation of stream function and vorticity */
function calc_psi_zeta(u, v, psi, zeta, flag, imax,  jmax,  delx,  dely)
{
    var i, j;

    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
    for (i=1;i<=imax-1;i++) {
        for (j=1;j<=jmax-1;j++) {
            if ( (flag[i][j] & C_F) && (flag[i+1][j] & C_F) &&
                (flag[i][j+1] & C_F) && (flag[i+1][j+1] & C_F)) {
                 zeta[i][j] = (u[i][j+1]-u[i][j])/dely
                            -(v[i+1][j]-v[i][j])/delx;
            } else {
                zeta[i][j] = 0.0;
            	}
        	}
    	}
		
	}	
