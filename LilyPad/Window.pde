/********************************************
  Class to hold display window parameters
    and functions to convert from
    pixels to other units

Example code:
void setup(){
  size(300,450);
  Window mine = new Window(1,1,10,15,40,80,200,300);
  rect(40,80,200,300);
  stroke(255,0,0);fill(0,255,255);
  rect(mine.px(0.5),mine.py(0.5),mine.x.r,mine.y.r);
  rect(mine.px(9.5),mine.py(14.5),mine.x.r,mine.y.r);
  line(mine.px(1),mine.py(1),mine.px(10),mine.py(15));
  rect(mine.px(mine.ix(0)),mine.py(mine.iy(449)),50,-50);
  rect(mine.px(mine.ix(299)),mine.py(mine.iy(0)),-50,50);
}
*********************************************/

class Window
{
	Scale x,y;
	int x0,y0,dx,dy;

	Window()
	{
		this( 0., 0., 1., 1., 0, 0, width, height );
	}

	// called by the basic example - only accept the no. cells in both directions
	// and start the inner grid from index 1. read the size of the window using
	// global variables width and height, assign both axes to start at 0 for the
	// display
	Window( int n, int m)
	{
		this( 1, 1, n-2, m-2, 0, 0, width, height );
	}

	Window( int n0, int m0, int dn, int dm)
	{
		this( n0, m0, dn, dm, 0, 0, width, height );
	}

	Window( int n0, int m0, int dn, int dm, int x0, int y0, int dx, int dy)
	{
		this(n0-0.5,m0-0.5,dn,dm,x0,y0,dx,dy);
	}

	Window( int n0, float m0, int dn, float dm, int x0, int y0, int dx, int dy)
	{
		this(n0-0.5,m0,dn,dm,x0,y0,dx,dy);
	}

	// construct from components
	Window( float n0, float m0, float dn, float dm, int x0, int y0, int dx, int dy)
	{
		// create the scale objects for both axes allowing interpolation onto the
		// display to be performed; inner scale is the CFD grid
		x = new Scale(n0,n0+dn,x0,x0+dx);
		y = new Scale(m0,m0+dm,y0,y0+dy);
		// assign the size and beginning of the display area
		this.x0 = x0;
		this.y0 = y0;
		this.dx = dx;
		this.dy = dy;
	}

	// convert from the cell index to pixel coordinate
	float ix(int i){ return x.in((float)i);}
	float iy(int i){ return y.in((float)i);}
	
	// ??? inverse conversion, used for something to do with mouse coordinates
	float idx(int i){ return i/x.r;}
	float idy(int i){ return i/y.r;}
	
	// convert a pixel location to a cell index
	int px(float i){ return (int)(x.out(i));}
	int py(float i){ return (int)(y.out(i));}
	
	// ??? inverse transfrom
	int pdx(float i){ return (int)(x.r*i);}
	int pdy(float i){ return (int)(y.r*i);}
	
	// Checks if the passed point fits within the limits of the display window.
	boolean inside( int x, int y )
	{
		return( x>=x0 && x<=x0+dx && y>=y0 && y<=y0+dy );
	}
}

// holds the scale for each axis, allowing for convertion between the actual grid
// and the display window
class Scale
{
	float inS,inE, // start and end of the inner scale
		outS,outE, // extents of the outer scale
		r; // ratio of units between the outer and inner scale

	Scale( float outS, float outE )
	{
		this(0,1,outS,outE);
	}
	
	Scale( float inS, float inE, float outS, float outE )
	{
		this.inS  = inS; // start of the inner scale
		this.inE  = inE; // end of the inner scale
		this.outS = outS; // start of the outer scale
		this.outE = outE; // end of the outer scale
		r = (outE-outS)/(inE-inS); // ratio between the units of the outer and inner scales
	}
	
	// converts from the inner to outer scale boundign the result between start and end of the outer scale
	float outB( float in ){ return out(min(max(in,inS),inE));}
	
	// converts from the inner to outer scale
	float out( float in ){ return (in-inS)*r+outS;} 
	
	// converts from the outer to inner scale
	float in( float out ){ return (out-outS)/r+inS;}
}
