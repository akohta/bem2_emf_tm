// -------------------------------------
lc=0.02; // element size 

rc0= 0.50;  // radius of circle 0
xc0=-0.55;   
yc0= 0.00;  // center of circle 0(xc,yc)

rc1= 0.25;  // radius of circle 1
xc1= 0.25;
yc1= 0.00;  // center of circle 1(xc,yc)
// -------------------------------------

Point(1)={xc0    ,yc0    ,0.0,lc}; // center of circle
Point(2)={xc0+rc0,yc0    ,0.0,lc};
Point(3)={xc0    ,yc0+rc0,0.0,lc};
Point(4)={xc0-rc0,yc0    ,0.0,lc};
Point(5)={xc0    ,yc0-rc0,0.0,lc};
Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};

Point( 6)={xc1    ,yc1    ,0.0,lc}; // center of circle
Point( 7)={xc1+rc1,yc1    ,0.0,lc};
Point( 8)={xc1    ,yc1+rc1,0.0,lc};
Point( 9)={xc1-rc1,yc1    ,0.0,lc};
Point(10)={xc1    ,yc1-rc1,0.0,lc};
Circle(5)={ 7,6, 8};
Circle(6)={ 8,6, 9};
Circle(7)={ 9,6,10};
Circle(8)={10,6,7};

Physical Line(99)={-1,-2,-3,-4,-5,-6,-7,-8}; // open region 
Physical Line( 1)={ 1, 2, 3, 4}; // Domain 1
Physical Line( 2)={ 5, 6, 7, 8}; // Domain 2

