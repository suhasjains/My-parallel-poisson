
//Inputs


//box dimension
boxdim = 1;

//cell size
gridsize = boxdim/11;


//defining points
Point(1) = {0,0,0,gridsize};
Point(2) = {boxdim,0,0,gridsize};
Point(3) = {boxdim,boxdim,0,gridsize};
Point(4) = {0,boxdim,0,gridsize};

//making lines using points
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};


//making a loop out of lines
Line Loop(9) = {5,6,7,8};

//making a surface out of line loop
Plane Surface(10) = 9;

//Make it structured
Transfinite Line{5,6,7,8} = boxdim/gridsize;
Transfinite Surface{10};
Recombine Surface{10};

//Creating physical lines
Physical Line("Bottom") = {5};
Physical Line("Right") = {6};
Physical Line("Top") = {7};
Physical Line("Left") = {8};
Show "*";
