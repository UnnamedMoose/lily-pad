/*********************************************************
solve the BDIM equation for velocity and pressure

  u = del*F+[1-del]*u_b+del_1*ddn(F-u_b)
  
  where 
    del is the zeroth moment of the smoothing kernel integral
    del_1 is the first moment (zero if mu1=false)
    u_b is the body velocity
    F is the fluid equation of motion:
    if(QUICK):
      F(u) = u(t)+\int_t^{t+dt} grad(u*u)+\mu*laplace(u)+g-grad(p)/rho \d t
    else(SEMI-LAGRANGIAN)
      F(u) = u(t,x(t))+\int_t^{t+dt} g-grad(p)/rho \d t
      
    where x(t) is the back-casted location of the grid points
      x(t) = x-\int_t^{t+dt} u \d t

Example code:

BDIM flow;
CircleBody body;
void setup(){
  size(400,400); 
  int n=(int)pow(2,6)+2; // multigrid solver needs n = power of 2 
  Window view = new Window(n,n);
  body = new CircleBody(n/3,n/2,n/8,view);
  flow = new BDIM(n,n,1.5,body);
}
void draw(){
  flow.update();  // project
  flow.update2(); // correct
  flow.u.vorticity().display(-0.75,0.75);
  body.display();
}
*********************************************************/
class BDIM
{
	int n,m; // number of cells in uniform grid
	float dt, nu, eps=2.0; // time step, viscosity, and BDIM kernel support width
	PVector g = new PVector(0,0); // gravity force
	VectorField u, // velocity field
		del, // zeroth moment of the smoothing kernel
		del1, // first moment of the smoothing kernel
		c,
		u0, // velocity field at the previous time step
		ub, // body velocity field
		wnx, // wall-normal direction to the nearest body
		wny,
		distance, // distance from the body surface
		rhoi; // inverse of local fluid density (1 for single-phase)
	Field p; // pressure field
	boolean QUICK, // whether to use the QUICK scheme in TODO
		mu1=true, // whether to use the 1st moment of the smoothing kernel
			//(see Maertens & Weymouth, "Accurate Cartesian-grid simulations of near-body flows at
			// intermediate Reynolds numbers", 2015)
		adaptive=false;

	// construct from components
	BDIM( int n, int m, float dt, Body body, VectorField uinit, float nu, boolean QUICK )
	{
		// set the constantns
		this.n = n;
		this.m = m;
		this.dt = dt;
		this.nu=nu;
		this.QUICK=QUICK;

		u = uinit;
		if(u.x.bval!=0) u.x.gradientExit = true; // use gradient exit for the axial velocity
		u0 = new VectorField(n,m,0,0);
		p = new Field(n,m);
		if(dt==0) setDt(); // adaptive time stepping for O(2) QUICK

		ub  = new VectorField(n,m,0,0);
		distance =  new VectorField(n, m, 10, 10);    
		del = new VectorField(n,m,1,1);
		del1 = new VectorField(n,m,0,0);
		rhoi = new VectorField(del); // initialise with ones
		c = new VectorField(del);
		wnx = new VectorField(n,m,0,0);
		wny = new VectorField(n,m,0,0);
		get_coeffs(body); // computes various geometric coeffs, like distance to walls, BDIM kernels, etc.
	}

	// construct without an initial velocity field, but given a free-stream x-value
	BDIM( int n, int m, float dt, Body body, float nu, boolean QUICK, int u_inf){
		this(n,m,dt,body,new VectorField(n,m,u_inf,0),nu,QUICK);}
	
	// construct using a default unit speed in the x-direction
	BDIM( int n, int m, float dt, Body body, float nu, boolean QUICK ){
		this(n,m,dt,body,new VectorField(n,m,1,0),nu,QUICK);}
  
	// If no body is supplied, create a body outside the domain
	BDIM( int n, int m, float dt, VectorField uinit, float nu, boolean QUICK ){
		this(n,m,dt,new CircleBody(-n/2,-m/2,n/10,new Window(0,0,n,m)),uinit,nu,QUICK);}

	// use default setup - unit x-speed, unit viscosity, no QUICK scheme
	BDIM( int n, int m, float dt, Body body){
		this(n,m,dt,body,new VectorField(n,m,1,0),1,false);}
 
	// advance the flow in time - 1st order time accuracy
	void update()
	{
		// O(dt,dx^2) BDIM projection step:
		c.eq(del.times(rhoi.times(dt))); // coefficient in front of the pressure gradient term in eq 33b and 34b in M&W'15

		u0.eq(u); // keep the current velocity for future use at thet next time step
		VectorField F = new VectorField(u); // this will hold the newly assembled flow equation of motion, including u_0

		// either solve simple convection or full conv-diff problem
		if(QUICK)
			F.AdvDif( u0, dt, nu );
		else
			F.advect( dt, u0 );

		updateUP( F, c );
	}

	// advance the flow in time - 2nd order time accuracy
	void update2()
	{
		// O(dt^2,dt^2) BDIM correction step:
		// us holds the divergence-free velocity field obtained during the first Euler step, therefore 1st order
		// accurate in time
		VectorField us = new VectorField(u);
		VectorField F = new VectorField(u);
	
		// this uses the QUICK interpolation on a segregated mesh to solve the advection-diffusion problem
		if(QUICK)
		{
			// evaluate the flow using BDIM starting from the velocity computed during the first step
			F.AdvDif( u0, dt, nu );
			updateUP( F, c );
		
			// Heun's corrector - take average of U(U_0) and U(U_1)
			u.plusEq(us); 
			u.timesEq(0.5);
		
			// adjust time step to maintain stability
			if(adaptive) dt = checkCFL();
		}
		// this advects each component of this field as a simple scalar variable
		else
		{
			F.eq(u0.minus(p.gradient().times(rhoi.times(0.5*dt))));
			F.advect(dt,us,u0);
			updateUP( F, c.times(0.5) );
		}
	}

	void updateUP( VectorField R, VectorField coeff )
	{
		/*  Seperate out the pressure from the forcing
			del*F = del*R+coeff*gradient(p)
		Approximate update (dropping ddn(grad(p))) which doesn't affect the accuracy of the velocity
			u = del*R+coeff*gradient(p)+[1-del]*u_b+del_1*ddn(R-u_b)
		Take the divergence
			div(u) = div(coeff*gradient(p)+stuff) = 0
		u.project solves this equation for p and then projects onto u
		*/
		R.plusEq(PVector.mult(g,dt)); // assemble the mu_0 term with u^0 in eq. 33a in M&W 2015; add body force
		u.eq(del.times(R).minus(ub.times(del.plus(-1)))); // 1st order BDIM using only zeroth moment - first part of eq 33a
		if(mu1) u.plusEq(del1.times((R.minus(ub)).normalGrad(wnx,wny))); // add the 2nd order terms - full eq. 33a
		u.setBC(); // set inflow/outflow and other BCs
		// make sure the velocity field is divergence free by adjusting the pressure correcition - eqns. 33b and 33c
		// the overloaded implementation of project is called which neglects the divergence which could
		// be created by changes to the geometry of the immersed body
		// TODO why?
		p = u.project(coeff,p);
	}

	// overloaded update methods which update coefficients if the immersed body is moving
	void update( Body body )
	{
		if(body.unsteady)
		{
			get_coeffs(body);
		}
		else
		{
			ub.eq(0.);
		}
		update();
	}

	void update2( Body body ){ // don't need to get coeffs again
		update2();}

	// computes various coefficients
	void get_coeffs( Body body )
	{
		get_dist(body); // distance from the body
		get_del(); // 0th BDIM kernel
		get_del1(); // 1st BDIM kernel
		get_ub(body); // body velocity
		get_wn(body); // wall-normal direction to the nearest body
	}

	// calculates the velocity of the body at each location in the flow
	void get_ub( Body body )
	{
		/* Immersed Velocity Field
		ub(x) = U(x)*(1-del(x))
		where U is the velocity of the body */
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				ub.x.a[i][j] = body.velocity(1,dt,(float)(i-0.5),j);
				ub.y.a[i][j] = body.velocity(2,dt,i,(float)(j-0.5));
			}
		}
	}

	// computes wall normal direction of the closest body point
	// TODO check where this gets used and make sure the distances are from face centres
	void get_wn(Body body)
	{
		PVector wn;
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				wn = body.WallNormal((float)(i-0.5),j);
				wnx.x.a[i][j]=wn.x;
				wny.x.a[i][j]=wn.y;
				wn = body.WallNormal(i,(float)(j-0.5));
				wnx.y.a[i][j]=wn.x;
				wny.y.a[i][j]=wn.y;
			}
		}   
	}
  
	// calculate distance between the body and cell centres
	void get_dist( Body body )
	{
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				distance.x.a[i][j] = body.distance((float)(i-0.5), j);
				distance.y.a[i][j] = body.distance(i, (float)(j-0.5));
			}
		}
	}

	// compute 0th kernel moment for each cell
	void get_del()
	{
		/* BDIM zeroth order moment
		del(x) = delta0(d(x))
		where d is the distance to the interface from x */
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				del.x.a[i][j] = delta0(distance.x.a[i][j]);
				del.y.a[i][j] = delta0(distance.y.a[i][j]);
			}
		}
		del.setBC();
	}

	// compute 1st kernel moment for each cell
	void get_del1()
	{
		/* BDIM first order moment
		del(x) = delta1(d(x))
		where d is the distance to the interface from x */
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				del1.x.a[i][j] = delta1(distance.x.a[i][j]);
				del1.y.a[i][j] = delta1(distance.y.a[i][j]);
			}
		}
		del1.setBC();
	}

	// calculate the 0th BDIM kernel value given a distance from the nearest body
	float delta0( float d )
	{
		if( d <= -eps ) // point inside the body
		{
			return 0;
		}
		else if( d >= eps ) // point well inside the flow
		{
			return 1;
		}
		else // point at the body-fluid interface, need to blend the solid and fluid equations
		{
			// eq. 16a in Maertens & Weymouth 2015
			return 0.5*(1.0+d/eps+sin(PI*d/eps)/PI);
		} 
	}

	// compute the 1st BDIM kernel
	float delta1( float d )
	{
		if( abs(d) >= eps) // outside of the blending region
		{
			return 0;
		}
		else // inside the blending region use eq. 16b from Maertens & Weymouth 2015
		{
			return 0.25*(eps-sq(d)/eps)-1/TWO_PI*(d*sin(d*PI/eps)+eps/PI*(1+cos(d*PI/eps)));
		} 
	}

	// return the lesser of CFL and 1
	float checkCFL(){ 
		return min(u.CFL(nu), 1);
	}

	// adjust the time step to keep the solver within stability limits
	void setDt()
	{
		dt = checkCFL();
		if(QUICK) adaptive = true;
	}
}
