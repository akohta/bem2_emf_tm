//--------------------
lc=0.02; // element size 

lx=0.5;  // width of x-direction
ly=0.5;  // width of y-direction

xc=0.0;  //
yc=0.0;  // center of rectangle
// -------------------

Point(1)={xc+0.5*lx,yc+0.5*ly,0.0,lc}; 
Point(2)={xc-0.5*lx,yc+0.5*ly,0.0,lc};
Point(3)={xc-0.5*lx,yc-0.5*ly,0.0,lc};
Point(4)={xc+0.5*lx,yc-0.5*ly,0.0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Physical Line(99)={-1,-2,-3,-4}; // open region
Physical Line( 1)={ 1, 2, 3, 4}; // domain 1

