#include "Vector.h"
#include <math.h>

/******************************************************************
rotmata(x,f,a):  Produce matrix a for the rotation about a vector
				 x by an angle f (in radians). Formerly ROTMAT.
******************************************************************/

void rotmata(double v[3], double angle, double amat[3][3])
	{
		double sangle, cangle;
		double x[3];
		double axis[5] = {0.,0.,1.,0.,0.};

		int i;

		sangle = sin(angle);
		cangle = cos(angle);

		unit(v,x);

		for (i=0;i<3;++i)
			subrot(x, &axis[2-i], sangle, cangle, amat[i]);

		return;
	}

/******************************************************************
a = fdot(x,y):  a equals dot product of vectors x,y.
******************************************************************/

double fdot(double x[3], double y[3])
	{
		double temp;

		temp = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];

		return temp;
	}

int idot(int x[3], int y[3])
	{
		return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
	}

/******************************************************************
subrot:  helper subroutine for rotmata. (see above).
******************************************************************/

void subrot(double v[3], double axis[3], double sinang, double cosang, double row[3])
	{
		double x[3], y[3], f;

		cross (v, axis, x);
		if ((magnitude(x)-1.0e-8) > 0)
			{
				cross(x,v,y);
	    		conmul(x, sinang, x);
	    		conmul(y, cosang, y);
	    		fvecsum(x,y,x);

				f = fdot(v, axis);

	    		conmul(v,f,y);
	    		unit(row, row);
	    		fvecsum (x,y,row);
	    		return;
			}
    	else
			cpyvec (axis, row);

    	return;
	}

/******************************************************************
identity2(a):  a is returned as an identity matrix.
******************************************************************/

void identity2(double a[3][3])
	{
    	int i,j;

    	for (i=0;i<3;++i)
			for (j=0;j<3;++j)
	    		if (i == j) a[i][j] = 1.;
					else a[i][j] = 0.;

    	return;
	}

/******************************************************************
cross(x,y,z):	z = x (cross) y. x,y,z are vectors.
******************************************************************/

void cross(double x[3], double y[3], double z[3])
	{
    	double t[3];
    	int i;

    	t[0] = x[1]*y[2]-x[2]*y[1];
    	t[1] = x[2]*y[0]-x[0]*y[2];
    	t[2] = x[0]*y[1]-x[1]*y[0];

    	for (i=0;i<3;++i) z[i]=t[i];

    	return;
	}

/******************************************************************
normto(a,b,c):  Produces a vector c normal to vectors a-0, b-0.
******************************************************************/

void normto(double a[3], double b[3], double c[3])
	{
    	double t1, t2, t3, t4, t5;
    	int i;
    	t1 = (a[0]*b[1]-b[0]*a[1]);
    	t2 = (a[1]*b[2]-a[2]*b[1]);
    	t3 = (a[0]*b[2]-b[0]*a[2]);
    	t4 = (a[0]*b[1]-b[0]*a[1]);

    	t5 = sqrt(t2*t2+t3*t3+t4*t4);

    	if (t5 == 0.)
			{
	    		for (i=0;i<3;++i) c[i]=0.;
	    		return;
			}

    	c[2] = t1/t5;
    	c[1] = -t3/t5;
    	c[0] = t2/t5;

    	return;
	}

/******************************************************************
f=vangle(a,b,c):  f is the angle (in degrees) between the vectors
		  a-b and b-c.  Formerly, angle.
******************************************************************/

double vangle(double a[3], double b[3], double c[3])
	{
    	double x[3], y[3];
    	double q, q2, q3;

    	vecdif(a,b,x);
    	vecdif(c,b,y);

    	q = magnitude(x)*magnitude(y);

    	if (q == 0.) return 0.;

    	q2 = fdot(x,y)/q;

		q3 = 57.29578*acos(q2);

    	return q3;
	}

/******************************************************************
unit(x,y):  Returns y as a unit vector in the direction of x.
******************************************************************/

void unit(double x[3], double y[3])
	{
		double d;

		d = magnitude(x);
		if (d == 0.) return;
		
		d = 1/d;

		conmul(x, d, y);
   		return;
	}

/******************************************************************
fvecsum(x,y,z):	z = x + y.  x,y,z are vectors.
******************************************************************/

void fvecsum(double x[3], double y[3], double z[3])
	{
		int i;
	
   		for (i=0;i<3;++i)
			z[i] = x[i] + y[i];

  		return;
	}

void ivecsum(int x[3], int y[3], int z[3])
	{
		int i;
	
   		for (i=0;i<3;++i)
			z[i] = x[i] + y[i];

  		return;
	}

/******************************************************************
unitry (a,b) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
conmul(x,f,y):  Returns y = f*x, where x,y are vectors, and f is a
		scaler.
******************************************************************/

void conmul(double x[3], double d, double y[3])
	{
   		int i;

	   	for (i=0; i<3; ++i)
			y[i] = x[i]*d;

	   	return;
	}

/******************************************************************
bij_uij(bij, rms, vib, amat) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
eigen(w, valu, vect) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
inver(a,b):  Matrix b is set to the inverse of matrix a.
******************************************************************/

double inver (double a[3][3], double b[3][3])
	{
		int idata1[3] = {1,2,0};
		int idata2[3] = {2,0,1};
		int id1, id2;

		double x[4][3];

		int j=0,i;

		for (i=0;i<3;++i)
			{
				id1 = idata1[i];
				id2 = idata2[i];
				cross(a[id1], a[id2], x[i]);
			}

		for (i=0;i<3;++i)
			{
				x[3][2] = fdot(a[i], x[i]);
				if (fabs(x[3][2])<1.0e-20) x[3][2]=1.0e-20;
				x[3][2] = 1./x[3][2];
				conmul (x[i],x[3][2],x[i]);
			}

		trnspz(x,b);

		return 0;
	}

/******************************************************************
trnspz(a,b):  Places the transpose of matrix a in matrix b.
******************************************************************/

void trnspz(double a[3][3], double b[3][3])
	{
    	int i, j;
		double t[3][3];

    	for (i=0;i<3;++i)
			for (j=0;j<3;++j)
	    		t[i][j]=a[j][i];

    	for (i=0;i<3;++i)
			for (j=0;j<3;++j)
	    		b[i][j]=t[i][j];

    	return;
	}

/******************************************************************
f = eigen1(a,v,b) and f = eigen2(a, v, b) NOT INCLUDED because they're not relevant
******************************************************************/

/******************************************************************
inrtia(x,f,a) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
cpymat(a,b):  Equates matrix b to matrix a. (b=a).
******************************************************************/

void cpymat(double a[3][3], double b[3][3])
	{
    	int i,j;

    	for (i=0;i<3;++i)
			for (j=0;j<3;++j)
	    		b[i][j] = a[i][j];

    	return;
	}

/******************************************************************
cpyvec(x,y):  Equates vector y to vector x. (y=x).
******************************************************************/

void cpyvec(double x[3], double y[3])
	{
		int i;

		for (i=0;i<3;++i)
			y[i] = x[i];

		return;
	}

/******************************************************************
matmul(a,b,c) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
vecmul(x,a,y):  Multiplies the matrix a by the vector x to produce
				vector y. (y=ax).
******************************************************************/
// assuming a[row index][col index] and y=a*x, x and y are column vectors
void vecmul(double x[3], double a[3][3], double y[3])
	{
		int i;
		double t[3];

		for (i=0;i<3;++i)
			t[i]=a[i][0]*x[0]+a[i][1]*x[1]+a[i][2]*x[2];

		for (i=0;i<3;++i) y[i]=t[i]; // XXX: waste of spacetime

		return;
	}

/******************************************************************
a = magnitude(x): a = magnitude of vector x.
******************************************************************/

double magnitude(double x[3])
	{
		double mg;

		mg = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

		return mg;
	}

/******************************************************************
vecdif(x,y,z): z = x - y.  x,y,z are vectors.
******************************************************************/

void vecdif(double x[3], double y[3], double z[3])
	{
		int i;

		for (i=0;i<3;++i)
			z[i] = x[i] - y[i];

		return;
	}

/******************************************************************
a=alngth(x,y) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
f = torsin(a,b,c,d) NOT INCLUDED because it's not relevant
******************************************************************/

/******************************************************************
transpose(a): in-place transposition of matrix a
******************************************************************/

void transpose(double a[3][3])
	{
		double b[3][3];
		int i, j;

		for (i=0;i<3;++i)
			for (j=0;j<3;++j)
				b[i][j] = a[j][i];

		for (i=0;i<3;++i)
			for (j=0;j<3;++j)
				a[i][j] = b[i][j];

      return;
	}

// multiplies vector a by a constant, stores in b
void fconmul(double* a, double factor, double* b, double a_len)
	{
		for (int i = 0; i < a_len; i++)
			b[i] = factor * a[i];
	}
