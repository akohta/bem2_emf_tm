// ------------------------
lc=0.02; // element size;

lx0= 0.5; // width of x-direction
ly0= 0.5; // width of y-direction
cx0=-0.3; //  
cy0= 0.0; // center of rectangle 0

lx1= 0.3; // width of x-direction
ly1= 0.3; // width of y-direction
cx1= 0.3; //          
cy1= 0.0; // center of rectangle 1
// ------------------------

Point(1)={cx0+0.5*lx0,cy0+0.5*ly0,0.0,lc};
Point(2)={cx0-0.5*lx0,cy0+0.5*ly0,0.0,lc};
Point(3)={cx0-0.5*lx0,cy0-0.5*ly0,0.0,lc};
Point(4)={cx0+0.5*lx0,cy0-0.5*ly0,0.0,lc};

Point(5)={cx1+0.5*lx1,cy1+0.5*ly1,0.0,lc};
Point(6)={cx1-0.5*lx1,cy1+0.5*ly1,0.0,lc};
Point(7)={cx1-0.5*lx1,cy1-0.5*ly1,0.0,lc};
Point(8)={cx1+0.5*lx1,cy1-0.5*ly1,0.0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};

Physical Line(99)={-1,-2,-3,-4,-5,-6,-7,-8}; // open region
Physical Line( 1)={ 1, 2, 3, 4}; // domain 1
Physical Line( 2)={ 5, 6, 7, 8}; // domain 2

