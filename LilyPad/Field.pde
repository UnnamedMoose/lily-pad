/**********************************
Field class

Holds the values and operators for a
fluid dynamics field

example code:
Field p,u,v;

void setup(){
  size(400,400);
  background(1,0,0);
  p = new Field(100,150);
  p.eq(1.0,40,50,40,75);
  u = new Field(100,150,1,0.25);
  v = new Field(100,150,2,0.2);
}
void draw(){
  p.advect( 1., u, v );
  p.display(0,1);
}

***********************************/

class Field
{
	// internal values of the field
	float[][] a;  
	
	int n,m, // number of grid points in each direction (x and y, respectively)
		btype = 0; // defines where the variable is defined in the segregated storage scheme
			// 0 - cell centre (e.g. pressure)
			// 1 - negative x-face (e.g. x-velocity)
			// 2 - negative y-face (e.g. y-velocity)

	float bval = 0; // BC value for Dirichlet BC
	boolean gradientExit = false; // whether to apply special outlet treatment at the RHS boundary (staggered variable only)

	// construct from components
	Field( int n, int m, int btype, float bval )
	{
		this.n = n;
		this.m = m;
		this.a = new float[n][m];
		this.btype = btype;
		this.bval = bval;
		// sets the internal field values to bval
		this.eq(bval);
	}
	
	// create the field with zero btype and bval
	Field( int n, int m)
	{
		this(n,m,0,0);
	}

	// create this field as a copy of a given field b
	// BUT field b is of different size than this field
	// and so skip the first is/js of its indices when
	// transferring the values
	Field(Field b, int is, int n, int js, int m )
	{
		this.n = n;
		this.m = m;
		a = new float[n][m];
		btype = b.btype;
		bval = b.bval;
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
		  		a[i][j] = b.a[i+is][j+js];
			}
		}
	}
	
	// create this field as a copy of field b
	Field(Field b)
	{
		n = b.n;
		m = b.m;
		a = new float[n][m];
		btype = b.btype;
		bval = b.bval;
		this.eq(b);
	}

	// guess what? this takes the Laplacian
	Field laplacian ()
	{
		// create an empty field d filled with zeros
		Field d = new Field( n, m );
		
		// loop over all INTERNAL cells (skip end stencils)
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				// use a 2nd order linear scheme in x- and y-directions
				// recall that the code is dimensionless with
				// the grid size being the length-scale
				// this means we skip the delta_x and delta_y
				// in the denominator
				d.a[i][j] = -4*a[i][j] + a[i+1][j]
					+ a[i-1][j] + a[i][j+1] + a[i][j-1];
			}
		}
		// leave the boundary values equal to zero
		return d;
	}

	// compute the gradient of the field
	VectorField gradient()
	{
		// this will hold the grad(a) values, initialise with zeros
		VectorField g = new VectorField(n,m,0,0);
		
		// loop over all internal faces
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				// compute the x and y derivatives using a
				// 1st order backward scheme
				g.x.a[i][j] = a[i][j]-a[i-1][j];
				g.y.a[i][j] = a[i][j]-a[i][j-1];
			}
		}
		
		// call the vectorField method to set the boundary values to zero in both directions (gradientExit=false)
		// this calls the Field::setBC() function on both x- and
		// y-components of the gradient function; x-dir is given
		// btype 1 and y-dir btype 2 - this causes the field values to be stored at the faces
		// as per staggered mesh formulation for a vector field
		g.setBC(); // issues?
		return g;
	}

	// advect the field - unwrap the vector fields into scalar components
	// u is the velocity field at the new time step, u0 is the old velocity
	void advect(float step, VectorField u, VectorField u0)
	{
		advect(step, u.x, u.y, u0.x, u0.y );
	
	}
	
	void advect(float step, Field u, Field v, Field u0, Field v0)
	{
		/* advect the field with the u, v velocity fields.
		Use a first order lagrangian method:
		Da/Dt = 0
		which translates into:
		a(t=dt,x,y) = a(t=0, x-dt*u(x,y), y-dt*v(x,y))
		
		The example code shows how diffusive this is 
		EDIT: by using an RK2 step to find x0 and using
		a quadratic interpolation, this method is second
		order and nondiffuse.
		
		EDIT2: by treating the old and new velocities 
		seperately, the method is second order in time.*/
		
		// create a copy of this field, keeping the current values
		Field a0 = new Field(this);
		
		// loop over all internal cells
		for( int i=1; i<n-1; i++)
		{
			for( int j=1; j<m-1; j++)
			{
				// select the interpolation point to be the +ve face of cell 0,
				// iterating up to -ve face of the last cell - i.e. all internal faces
				float x = i;
				float y = j;
				// correct for boundary values if necessary
				// Note: linear and quadratic interpolation methods will modify x and y by +0.5
				// if the field is stored at the -ve face centres because we want to convert the x-
				// and y-location into indices in the matrix. The interpolation method is generic and
				// so also makes sure the point at which interpolation is performed is valid, bounded inside
				// the domain and so on. If btype>0 the variable is stored at face centre, and hence 0.5\delta
				// before/below the (x,y) of the (i,j) cell
				if(btype==1) x -= 0.5;
				if(btype==2) y -= 0.5;
				// get the velocities at this time step
				float ax = -step*u.linear( x, y );
				float ay = -step*v.linear( x, y );
				// get the velocities at the previous time step
				// account for the fact that the interpolation point would have
				// "moved" back then as well
				float bx = -step*u0.linear( x+ax, y+ay );
				float by = -step*v0.linear( x+ax, y+ay );
				// assume that at the new time step field value is equal to
				// that at the previous time step at a location -dT*(u,v)
				// perform Runge-Kutta interpolation using three data points
				// to get the location in the field which would have been the current
				// value one time step ago
				// then use quadratic interpolation to get the field value there and
				// assign it to be the field value after update
				a[i][j] = a0.quadratic( x+0.5*(ax+bx), y+0.5*(ay+by) );
			}
		}
		// correct the BCs
		setBC();
	}
	
	// advect the field at the first time step given the current
	// velocity field; unwrap vectors into scalars
	void advect(float step, VectorField u)
	{
		advect(step, u.x, u.y);
	}
	
	void advect(float step, Field u, Field v)
	{
		/* advect the field with the u,v velocity field
		here we use a first order lagrangian method
		Da/Dt = 0
		a(t=dt,x,y) = a(t=0,x-dt*u(x,y),y-dt*v(x,y))
		the example code shows how diffusive this is 
		EDIT: by using an RK2 step to find x0 and using
		a quadratic interpolation, this method is second
		order and nondiffuse.
		EDIT2: RK2 step is not needed for first step
		in time accurate soln.*/
		
		// create a copy of the field
		Field a0 = new Field(this);
		
		// loop over all internal values
		for( int i=1; i<n-1; i++)
		{
			for( int j=1; j<m-1; j++)
			{
				// select interpolation point to be this cell
				float x = i;
				float y = j;
				// correct if the field is stored at face centres
				if(btype==1) x -= 0.5;
				if(btype==2) y -= 0.5;
				// get the velocities at this point, multiply by -dt
				float ax = -step*u.linear( x, y );
				float ay = -step*v.linear( x, y );
				// assume that at the new time step field value is equal to
				// that at the previous time step at a location -dT*(u,v)
				// perform quadratic interpolation for accuracy
				a[i][j] = a0.quadratic( x+ax, y+ay );
			}
		}
		// correct the BCs
		setBC();
	}

	// perform quadratic interpolation at (x0,y0)
	float quadratic( float x0, float y0)
	{
		// correct (x0,y0) if the field is stored at face centres
		float x = x0, y = y0;
		if(btype==1) x += 0.5;
		if(btype==2) y += 0.5;
	
		// get the indices of cells we need
		int i = round(x), j = round(y);
	
		// make sure we won't exceed the grid dimensions
		// perform linear interpolation if we do
		if( i>n-2 || i<1 || j>m-2 || j<1 )
			return linear( x0, y0 );
	
		// to interpolate we need x and y expressed about
		// the i-th and j-th cell
		x -= i;
		y -= j;
	
		// interpolate in the x-direction at 3 y-positions
		float e = quadratic1D(x,a[i-1][j-1],a[i][j-1],a[i+1][j-1]);
		float f = quadratic1D(x,a[i-1][j  ],a[i][j  ],a[i+1][j  ]);
		float g = quadratic1D(x,a[i-1][j+1],a[i][j+1],a[i+1][j+1]);
	
		// use the interpolated values to interpolate in the y-direction
		return quadratic1D(y,e,f,g);
	}

	// perform quadratic interpolation in 1D given a distance
	// x from the central point, as well as field values at
	// indices (i-1), i, and (i+1);
	float quadratic1D(float x, float e, float f, float g)
	{
		// compute the value
		float x2 = x*x;
		float fx = f*(1.-x2);
		fx += (g*(x2+x)+e*(x2-x))*0.5;
		// bound the scheme between the specified field values
		fx = min(fx,max(e,f,g));
		fx = max(fx,min(e,f,g));
		return fx;
	} 

	// Interpolate the field at point (x0,y0)
	float linear(float x0, float y0)
	{
		// take the point of interest and make it bounded between
		// 0.5 and x_max-0.5 - i.e. cell centres of 0th and N-1st cell (boundary cells)
		float x  = min(max(0.5,x0), n-1.5);
		// correct if the field is stored at face centres
		if(btype==1) x += 0.5;
	
		// index of the cell around which we interpolate
		// make sure we do not pick N-1st cell (the last one)
		// casting <float>x to <int> will round down to the nearest integer
		// so that x should be between i and i+1
		int i = min( (int)x, n-2 ); 
		// weighting factor in the x-direction
		float s = x-i;
	
		float y  = min(max(0.5,y0), m-1.5);
		// correct if the field is stored at face centres
		if(btype==2) y += 0.5;
	
		int j = min( (int)y, m-2 );
		float t = y-j;
	
		// if weighting factors are zero then just pick the field value (i,j)
		if(s==0 && t==0)
		{
			return a[i][j];
		}
		// perform linear interpolation
		// first do in the y-direction at both x-locations we need
		// then take the results and interpolate in the x-direction
		else
		{
			return  s*(t*a[i+1][j+1]+(1-t)*a[i+1][j])+
				(1-s)*(t*a[i  ][j+1]+(1-t)*a[i  ][j]);
		}
	}
  
	// simple wrapper around the linear interpolation method
	float interp( float x0, float y0 )
	{
		return linear(x0,y0);
	}

	// show this field as an image given the maximum and minimum values between
	// which to assign the colour scale
	void display( float low, float high )
	{
		PImage img = createImage(n,m,RGB);
		img.loadPixels();
		for ( int i=0 ; i<n ; i++ )
		{
			for ( int j=0 ; j<m ; j++ )
			{
				float f = a[i][j];
				int k = i+j*n;
				f = map(f,low,high,0,255);
				img.pixels[k] = color(f);
			}
		}
		img.updatePixels();
		image(img,0,0,width,height);
	}

	// correct the boundary conditions
	void setBC ()
	{
		float s=0;
	
		// go over all y-indices
		for (int j=0 ; j<m ; j++ )
		{
			// set the 0th and n-1st (first and last) field values
			// to what they are one index into the domain
			// -> apply zero-gradient (von Neumann) BC at left and right sides
			a[0][j]   = a[1][j];
			a[n-1][j] = a[n-2][j];
		
			// if the field values are stored at the faces
			if(btype==1)
			{
				// if we have a zero-gradient type right-hand-side boundary
				// but fixed-value inlet then prescribe the value only there
				if(gradientExit)
				{
					// set inlet to fixed-value
					a[1][j] = bval;
					// s now holds an integral of the field values at the RHS boundary
					if(j>0 & j<m-1)
						s += a[n-1][j];          
				}
				// set the field values to the constant held in the class
				// if we have a fixed-value type BC set
				// do this to cells 1st cell into the domain at the LHS and last one at the RHS
				// this ensures zero-gradient at the LHS is maintained later on at the LHS
				else
				{
					a[1][j]   = bval;  
					a[n-1][j] = bval;
				}
			}
		}
	
		// go over all x-indices
		for (int i=0 ; i<n ; i++ )
		{
			// apply zero-gradient in the y-direction at top and bottom
			a[i][0]   = a[i][1];
			a[i][m-1] = a[i][m-2];
		
			// if the field values are stored at face centres set boundary values
			if(btype==2)
			{
				a[i][1]   = bval;  
				a[i][m-1] = bval;
			}
		}
	
		// set the RHS BC // TODO this needs more clarity
		if(gradientExit)
		{
			// get the average values of s
			s /= float(m-2);
			// remove the mean value from the RHS cell and the boundary value
			for( int j=1; j<m-1; j++ )
				a[n-1][j] += bval-s;
		}
	}

	// compute gradient in a normal direction
	// wnx and wny are wall-normal x- and y-directions with their
	// VectorField.x holding the actual values in x- and y-dirs
	Field normalGrad(VectorField wnx, VectorField wny)
	{
		// result field
		Field g = new Field(n,m,0,0);
	
		// loop over all internal cells
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				// compute the wall-normal gradient
				// approximate d(a)/dx_i using 2nd order central scheme
				// [a(i+1,j)-a(i-1,j)]/(2 delta_x), where delta_x = 1.0
				// take a dot-product with the wall-normal direction
				g.a[i][j] = 0.5*( wnx.x.a[i][j] * (a[i+1][j]-a[i-1][j])
					+ wny.x.a[i][j] * (a[i][j+1]-a[i][j-1]) );
			}
		}
		return g;
	}
	
	//------------------------
	// Operators
	
	// multiply by a scalar
	Field times( float b )
	{
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] *= b;
			}
		}
		return c;
	}
	
	// multiply this field by a given field item-wise
	Field times( Field b )
	{
		// if the field storage types do not match throw an error
		mismatch(this.btype,b.btype); 
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] *= b.a[i][j];
			}
		}
		return c;
	}
	
	// multiply by a scalar and set this field to the result
	void timesEq( float b )
	{
		eq(times(b));
	}
	
	// multiply item-wise with a field b and set this to the result
	void timesEq( Field b )
	{
		eq(times(b));
	}
	
	// add a scalar
	Field plus( float b )
	{
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] += b;
			}
		}
		return c;
	}
	
	// add to a field item-wise
	Field plus( Field b )
	{
		mismatch(this.btype,b.btype); 
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] += b.a[i][j];
			}
		}
		return c;
	}
	
	// add a scalar to this field and set it to the result
	void plusEq( float b )
	{
		eq(plus(b));
	}
	
	// add a field item-wise and store the result
	void plusEq( Field b )
	{
		eq(plus(b));
	}
	
	// subtract field
	Field minus( Field b )
	{
		mismatch(this.btype,b.btype); 
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] -= b.a[i][j];
			}
		}
		return c;
	}
	
	// subtract a scalar
	void minusEq( Field b )
	{
		eq(minus(b));
	}
	
	// return 1/a_ij (NOT an inverse matrix)
	Field inv()
	{
		Field c = new Field(this);
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				c.a[i][j] = 1./c.a[i][j];
			}
		}
		return c;
	}
	
	// invert each element of this field
	void invEq()
	{
		eq(inv());
	}
	
	// do an inner (dot) product with a given field
	// use ONLY internal field values
	float inner(Field b)
	{
		mismatch(this.btype,b.btype); 
		float s = 0;
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				s += a[i][j]*b.a[i][j];
			}
		} 
		return s;
	}
	
	// sum all INTERNAL values of this field
	float sum()
	{
		float s = 0;
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				s += a[i][j];
			}
		}
		return s;
	}
  
	// assign the values of the field to a given value, b,
	// between the start and end indices is/js and ie/je
	void eq(float b ,int is, int ie, int js, int je)
	{
		for(int i=is; i<ie; i++)
		{
			for(int j=js; j<je; j++)
			{
				a[i][j] = b;
			}
		}
	}
	
	// set the entire field value equal to a constant b
	void eq(float b)
	{
		eq(b,0,n,0,m);
	}
	
	// transfer values of the given field b
	void eq(Field b)
	{
		for( int i=0; i<n; i++)
		{
			for( int j=0; j<m; j++)
			{
				a[i][j] = b.a[i][j];
			}
		}
	}

	// used to compare btype of different fields
	void mismatch(int i, int j)
	{
		if(i!=j)
		{
			println("You can't add or multiply fields of different types.");
			exit();
		}     
	}
}
