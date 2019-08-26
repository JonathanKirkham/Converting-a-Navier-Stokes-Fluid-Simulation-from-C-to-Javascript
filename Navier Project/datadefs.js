var C_B =    0x0000;   /* This cell is an obstacle/boundary cell */
var B_N =    0x0001;  /* This obstacle cell has a fluid cell to the north */
var B_S =    0x0002;   /* This obstacle cell has a fluid cell to the south */
var B_W =    0x0004;   /* This obstacle cell has a fluid cell to the west */
var B_E =    0x0008;   /* This obstacle cell has a fluid cell to the east */
var B_NW =    (B_N | B_W);
var B_SW =   (B_S | B_W);
var B_NE =    (B_N | B_E);
var B_SE =    (B_S | B_E);
var B_NSEW =  (B_N | B_S | B_E | B_W);

var C_F =    0x0010;   /* This cell is a fluid cell */

/* Macros for poisson(), denoting whether there is an obstacle cell
 * adjacent to some direction*/

eps_E = function(flag,i,j) { return ((flag[i+1][j] & C_F)?1:0); }
eps_W = function(flag,i,j) { return ((flag[i-1][j] & C_F)?1:0); }
eps_N = function(flag,i,j) { return ((flag[i][j+1] & C_F)?1:0); }
eps_S = function(flag,i,j) { return ((flag[i][j-1] & C_F)?1:0); }
