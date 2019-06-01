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

  Field laplacian (){
    Field d = new Field( n, m );
    d.btype = btype;
    for ( int i=1 ; i<n-1 ; i++ ) {
    for ( int j=1 ; j<m-1 ; j++ ) {
      d.a[i][j] = -4*a[i][j]+a[i+1][j]
        +a[i-1][j]+a[i][j+1]+a[i][j-1];
    }}
    return d;
  }

  VectorField gradient(){
    mismatch(btype,0);
    VectorField g = new VectorField(n,m,0,0);
    for ( int i=1 ; i<n-1 ; i++ ) {
    for ( int j=1 ; j<m-1 ; j++ ) {
      g.x.a[i][j] = a[i][j]-a[i-1][j];
      g.y.a[i][j] = a[i][j]-a[i][j-1];
    }}
    g.setBC(); // issues?
    return g;
  }

  VectorField curl (){
    // returns curl{this \hat z}
    mismatch(btype,3);
    VectorField g = new VectorField(n,m,0,0);
    for ( int i=1 ; i<n-1 ; i++ ) {
    for ( int j=1 ; j<m-1 ; j++ ) {
      g.x.a[i][j] = a[i][j+1]-a[i][j];
      g.y.a[i][j] = a[i][j]-a[i+1][j];
    }}
    return g;
  }

  void advect( float step, VectorField u, VectorField u0 ){
    advect( step, u.x, u.y, u0.x, u0.y );
  }
  void advect( float step, Field u, Field v, Field u0, Field v0 ){
    /* advect the field with the u,v velocity field
       here we use a first order lagrangian method
         Da/Dt = 0
         a(t=dt,x,y) = a(t=0,x-dt*u(x,y),y-dt*v(x,y))
       the example code shows how diffusive this is
       EDIT: by using an RK2 step to find x0 and using
       a quadratic interpolation, this method is second
       order and nondiffuse.
       EDIT2: by treating the old and new velocities
       seperately, the method is second order in time.*/
    Field a0 = new Field(this);
    for( int i=1; i<n-1; i++){
      for( int j=1; j<m-1; j++){
        float x = i;
        float y = j;
        if(btype==1) x -= 0.5;
        if(btype==2) y -= 0.5;
        float ax = -step*u.linear( x, y );
        float ay = -step*v.linear( x, y );
        float bx = -step*u0.linear( x+ax, y+ay );
        float by = -step*v0.linear( x+ax, y+ay );
        a[i][j] = a0.quadratic( x+0.5*(ax+bx), y+0.5*(ay+by) );
      }
    }
    setBC();
  }
  void advect( float step, VectorField u ){
    advect( step, u.x, u.y );
  }
  void advect( float step, Field u, Field v){
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
    Field a0 = new Field(this);
    for( int i=1; i<n-1; i++){
      for( int j=1; j<m-1; j++){
        float x = i;
        float y = j;
        if(btype==1|btype==3) x -= 0.5;
        if(btype==2|btype==3) y -= 0.5;
        float ax = -step*u.linear( x, y );
        float ay = -step*v.linear( x, y );
        a[i][j] = a0.quadratic( x+ax, y+ay );
      }
    }
    setBC();
  }

  float quadratic( float x0, float y0){
    float x = x0, y = y0;
    if(btype==1|btype==3) x += 0.5;
    if(btype==2|btype==3) y += 0.5;
    int i = round(x), j = round(y);
    if( i>n-2 || i<1 || j>m-2 || j<1 )
      return linear( x0, y0 );
    x -= i; y -= j;
    float e = quadratic1D(x,a[i-1][j-1],a[i][j-1],a[i+1][j-1]);
    float f = quadratic1D(x,a[i-1][j  ],a[i][j  ],a[i+1][j  ]);
    float g = quadratic1D(x,a[i-1][j+1],a[i][j+1],a[i+1][j+1]);
    return quadratic1D(y,e,f,g);
  }
  float quadratic1D(float x, float e, float f, float g){
    float x2 = x*x;
    float fx = f*(1.-x2);
    fx += (g*(x2+x)+e*(x2-x))*0.5;
    fx = min(fx,max(e,f,g));
    fx = max(fx,min(e,f,g));
    return fx;
  }
  float linear( float x0, float y0 ){
    float x  = min(max(0.5,x0), n-1.5);
    if(btype==1|btype==3) x += 0.5;
    int i = min( (int)x, n-2 );
    float s = x-i;
    float y  = min(max(0.5,y0), m-1.5);
    if(btype==2|btype==3) y += 0.5;
    int j = min( (int)y, m-2 );
    float t = y-j;
    if(s==0 && t==0){
      return a[i][j];
    }else{
      return s*(t*a[i+1][j+1]+(1-t)*a[i+1][j])+
         (1-s)*(t*a[i  ][j+1]+(1-t)*a[i  ][j]);
    }
  }
  float interp( float x0, float y0 ){return linear(x0,y0);}

  void display( float low, float high ){
    PImage img = createImage(n-2,m-2,RGB);
    img.loadPixels();
    for ( int i=0 ; i<n-2 ; i++ ) {
      for ( int j=0 ; j<m-2 ; j++ ) {
        float f = a[i+1][j+1];
        int k = i+j*(n-2);
        f = map(f,low,high,0,255);
        img.pixels[k] = color(f);
      }
    }
    img.updatePixels();
    int x0 = width/(n-2), y0 = height/(m-2);
    image(img,x0,y0,width,height);
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

  Field times( float b ){
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
        c.a[i][j] *= b;
    }}
    return c;
  }
  Field times( Field b ){
    mismatch(this.btype,b.btype);
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      c.a[i][j] *= b.a[i][j];
    }}
    return c;
  }
  void timesEq( float b ){ eq(times(b)); }
  void timesEq( Field b ){ eq(times(b)); }
  Field plus( float b ){
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      c.a[i][j] += b;
    }}
    return c;
  }
  Field plus( Field b ){
    mismatch(this.btype,b.btype);
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      c.a[i][j] += b.a[i][j];
    }}
    return c;
  }
  void plusEq( Field b ){ eq(plus(b)); }
  void plusEq( float b ){ eq(plus(b)); }
  Field minus( Field b ){
    mismatch(this.btype,b.btype);
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      c.a[i][j] -= b.a[i][j];
    }}
    return c;
  }
  void minusEq( Field b ){ eq(minus(b)); }
  Field inv(){
    Field c = new Field(this);
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      c.a[i][j] = 1./c.a[i][j];
    }}
    return c;
  }
  void invEq(){ eq(inv()); }
  float inner( Field b ){
    mismatch(this.btype,b.btype);
    double s = 0;
    for ( int i=1 ; i<n-1 ; i++ ) {
    for ( int j=1 ; j<m-1 ; j++ ) {
      s += a[i][j]*b.a[i][j];
    }}
    return (float)s;
  }
  float sum(){
    float s = 0;
    for ( int i=1 ; i<n-1 ; i++ ) {
    for ( int j=1 ; j<m-1 ; j++ ) {
      s += a[i][j];
    }}
    return s;
  }
  void eq( float b ,int is, int ie, int js, int je){
    for( int i=is; i<ie; i++){
    for( int j=js; j<je; j++){
      a[i][j] = b;
    }}
  }
  void eq( float b){eq(b,0,n,0,m);}
  void eq( Field b){
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      a[i][j] = b.a[i][j];
    }}
  }

  void mismatch(int i, int j){
    if(i!=j){
      println("You can't add or multiple fields of different types.");
      exit();
    }
  }

  float L_inf(){
    float mx = 0;
    for( int i=0; i<n; i++){
    for( int j=0; j<m; j++){
      mx = max(mx,abs(a[i][j]));
    }}
    return mx;
  }
}
