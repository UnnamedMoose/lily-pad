/**********************************
PoissonMatrix class

Holds values and operators for a Poisson matrix

   A*x = div{lower*grad{x}} = b      (1)

where lower is a VectorField which are the
lower diagonal components of the matrix A.
From Eq.1 the matrix is symmetric and the
diagonals are given by

  sum{a_{i,j},i} = 0 forall j       (2)

example code:

void setup(){
  colorMode(RGB,1.0);
  size(500,500);
  VectorField c = new VectorField(40,40,0,0);
  c.x.eq(1.); c.y.eq(1.); c.setBC();
  PoissonMatrix A = new PoissonMatrix(c);
  Field x = new Field(40,40,0,1);
  for(int i=0; i<40; i++ ){
  for(int j=0; j<40; j++ ){
    x.a[i][j] = i*i-j*j;
  }
  }
  Field Ax = A.times(x);
  Ax.display(-2,2);
}
***********************************/
public class PoissonMatrix
{
	int n, m; // size of the matrix
	VectorField lower; // lower diagonal coefficients
	Field diagonal,inv;

	// construct given the coefficients matrix
	PoissonMatrix( VectorField lower )
	{
		this.n = lower.n;
		this.m = lower.m;
		this.lower = new VectorField(lower);
		diagonal = new Field(n,m);
		inv = new Field(n,m,0,1);

		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				// sum values at all faces of the cell and assign them to the diagonal
				// TODO what do the numbers mean?
				float sumd = lower.x.a[i][j] + lower.x.a[i+1][j] + lower.y.a[i][j] + lower.y.a[i][j+1];
				diagonal.a[i][j] = -sumd;
				// compute inverse of the diagonal, be careful not to divide by zero
				if(sumd>1e-5) inv.a[i][j] = -1./sumd;
			}
		}
	}

	// TODO not entirely sure what this does
	Field times( Field x )
	{
		Field ab = new Field(n,m);
		for ( int i=1 ; i<n-1 ; i++ )
		{
			for ( int j=1 ; j<m-1 ; j++ )
			{
				ab.a[i][j] = x.a[i  ][j  ]*diagonal.a[i][j] // multiply diagonal terms
							+x.a[i-1][j  ]*lower.x.a[i  ][j  ] // TODO what are these
							+x.a[i+1][j  ]*lower.x.a[i+1][j  ]
							+x.a[i  ][j-1]*lower.y.a[i  ][j  ]
							+x.a[i  ][j+1]*lower.y.a[i  ][j+1];
			}
		}
		return ab;
	}

	// clalculate residual TODO check MG.pde to see how exactly this works
	Field residual( Field b, Field x ) {
		return b.plus(this.times(x).times(-1)); }

  Field residual( Field b, Field x ){
    return b.minus(this.times(x));
  }
}
