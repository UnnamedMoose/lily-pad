// Class to hold the ortho-normal project terms of a line segment

class OrthoNormal
{
	float l, // length of the line segment
		nx, // normal vector x-component
		ny, // normal vector y-component
		tx, // x-tangent - unit vector between end points of the line segment
		ty, // y-tangent
		off, // normal offset - distance from the origin to the plane defined by this line segment
		t1, // distance of the end points projected onto the tangent
		t2;
	
	PVector cen; // centre of the line segment

	// default constructor
	OrthoNormal()
	{
		this(new PVector(0,0), new PVector(0,1));
	}
  
	// set the ortho-normal values based on two points
	OrthoNormal(PVector x1, PVector x2 )
	{
		l = PVector.sub(x1,x2).mag();
		tx = (x2.x-x1.x)/l;    // x tangent
		ty = (x2.y-x1.y)/l;    // y tangent
		t1 = x1.x*tx+x1.y*ty;  // tangent location of point 1
		t2 = x2.x*tx+x2.y*ty;  // tangent location of point 2
		nx = -ty;
		ny = tx;     // normal vector
		off = x1.x*nx+x1.y*ny; // normal offset
		cen = PVector.add(x1,x2); // centriod
		cen.div(2.);
	}

	// distance to a point defined by (x,y)
	float distance( float x, float y, Boolean projected)
	{
		float d = x*nx+y*ny-off; // normal distance to line 
		if(projected) return d;  // |distance|_n (signed, fastest)
		
		float d1 = x*tx+y*ty-t1; // tangent dis to start
		float d2 = x*tx+y*ty-t2; // tangent dis to end
		
		//    return sqrt(sq(d)+sq(max(0,-d1))+sq(max(0,d2))); // |distance|_2
		return abs(d)+max(0,-d1)+max(0,d2);              // |distance|_1 (faster)
	}
	
	// default to signed dist - wrapper calling the base implementation in projected mode
	float distance( float x, float y)
	{
		return distance(x,y,true);
	}

	// tangent coordinate
	float tanCoord( float x, float y )
	{
		return min( max( (x*tx+y*ty-t1)/l, 0), 1 );
	}

	// move the line in space
	void translate(float dx, float dy)
	{
		t1 += dx*tx+dy*ty;
		t2 += dx*tx+dy*ty;
		off += dx*nx+dy*ny;
		cen.add(dx,dy,0);
	}

	void print()
	{
		println("t=["+tx+","+ty+"]");
		println("n=["+nx+","+ny+"]");
		println("offsets=["+t1+","+t2+","+off+"]");
	}
}