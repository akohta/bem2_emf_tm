// ---------------------------
lc0=0.01;  // element size 0
lc1=0.005; // element size 1

rc=1.0; // radius of half circle

xc=0.001;
yc=0.0; // center of half circle

t=0.05; // thickness of layer 
l0=0.1; // center hole size
l1=0.2; // outside hole size 
// --------------------------

Point(1)={xc   ,yc   ,0.0,lc1}; // center of circle
Point(2)={xc   ,yc-rc,0.0,lc1};
Point(3)={xc-rc,yc   ,0.0,lc0};
Point(4)={xc   ,yc+rc,0.0,lc1};
Point(5)={xc+t ,yc+rc,0.0,lc1};
Point(6)={xc+t ,yc+0.5*l1,0.0,lc1};
Point(7)={xc   ,yc+0.5*l0,0.0,lc1};
Point(8)={xc   ,yc-0.5*l0,0.0,lc1};
Point(9)={xc+t ,yc-0.5*l1,0.0,lc1};
Point(10)={xc+t, yc-rc,0.0,lc1};

Circle(1)={3,1,2};
Circle(2)={4,1,3};
Line(3)={5,4};
Line(4)={6,5};
Line(5)={7,6};
Line(6)={8,7};
Line(7)={9,8};
Line(8)={10,9};
Line(9)={2,10};
Line(10)={7,4};
Line(11)={2,8};

Physical Line(99)={-1,-2,-3,-4,-5,-6,-7,-8,-9}; // open region
Physical Line( 1)={ 1, 2,10, 6,11}; // domain 1
Physical Line( 2)={ 3, 4, 5,-10, 7, 8, 9,-11}; // Domain 2
