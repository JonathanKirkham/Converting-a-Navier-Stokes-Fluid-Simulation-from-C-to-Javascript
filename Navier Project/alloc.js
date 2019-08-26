// Initialise Float32Array.
// Allocate Array of Floating Point Arrays
function alloc_matrix(rows, cols) {
    var matrix = [];
    var i, j;
    for (i = 0; i < rows; i++) {
	//var row = new Float32Array(cols);
	var row = new Float64Array(cols); 
	matrix.push(row);
    }
    return matrix;
}
// Allocate Array of Javascript Arrays
/* function alloc_matrix(rows, cols) {
    var matrix = [];
    var i, j;
    for (i = 0; i < rows; i++) {
	var row = [];
	for (j = 0; j < cols; j++) {
	    row.push(0);
	}
	matrix.push(row);
    }
return matrix; } */