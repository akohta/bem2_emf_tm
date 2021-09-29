// ------------------------------
lc=0.02; // element size 

ri=0.2;  // radius of inner circle
xi=0.25;   
yi=0.0;  // center of inner circle

ro=0.5;  // radius of outer circle
xo=0.0;
yo=0.0;  // center of outer circle
// -------------------------------

Point( 1)={xi   ,yi   ,0.0,lc}; // center
Point( 2)={xi+ri,yi   ,0.0,lc};
Point( 3)={xi   ,yi+ri,0.0,lc};
Point( 4)={xi-ri,yi   ,0.0,lc};
Point( 5)={xi   ,yi-ri,0.0,lc};
Point( 6)={xo   ,yo   ,0.0,lc}; // center
Point( 7)={xo+ro,yo   ,0.0,lc};
Point( 8)={xo   ,yo+ro,0.0,lc};
Point( 9)={xo-ro,yo   ,0.0,lc};
Point(10)={xo   ,yo-ro,0.0,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};
Circle(5)={7,6,8};
Circle(6)={8,6,9};
Circle(7)={9,6,10};
Circle(8)={10,6,7};

Physical Line(99)={-5,-6,-7,-8}; // open region
Physical Line( 1)={ 1, 2, 3, 4}; // Domain 1
Physical Line( 2)={-1,-2,-3,-4, 5, 6, 7, 8}; // Domain 2

