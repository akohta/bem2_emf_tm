// -------------------------------------
lc=0.02; // element size 

rc=0.5;  // radius of circle

xc=0.0;   
yc=0.0;  // center of circle (xc,yc)
// -------------------------------------


Point(1)={xc   ,yc   ,0.0,lc}; // center of circle
Point(2)={xc+rc,yc   ,0.0,lc};
Point(3)={xc   ,yc+rc,0.0,lc};
Point(4)={xc-rc,yc   ,0.0,lc};
Point(5)={xc   ,yc-rc,0.0,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};

//Physical Line(99)={-1,-2,-3,-4}; // Domain 0 (open region)
Physical Line( 1)={ 1, 2, 3, 4}; // Domain 1 
Physical Line(99)={-1,-2,-3,-4}; // Domain 0 (open region)

