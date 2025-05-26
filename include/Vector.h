#ifndef VECTOR_H
#define VECTOR_H

// prototypes for the functions contained in vector.cpp

void make_rotmat(double rangle1, double a[3][3], int axis);
void matmul(double a[3][3], double b[3][3], double c[3][3]);
void vecmul(double x[3], double a[3][3], double y[3]);
void vecdif(double x[3], double y[3], double z[3]);
double magnitude(double x[3]);
void unit(double x[3], double y[3]);
void fvecsum(double x[3], double y[3], double z[3]);
void conmul(double x[3], double d, double y[3]);
void cpyvec(double x[3], double y[3]);
void cpymat(double a[3][3], double b[3][3]);
double vangle (double a[3], double b[3], double c[3]);
void identity2(double a[3][3]);
void cross(double x[3], double y[3], double z[3]);
void rotmata(double v[3], double angle, double amat[3][3]);
void normto(double a[3], double b[3], double c[3]);
void subrot(double v[3], double axis[3], double sinang, double cosang, double row[3]);
double fdot(double x[3], double y[3]);
double alngth(double x[3], double y[3]); 
double torsin (double a[3], double b[3], double c[3], double d[3]);
void transpose(double a[3][3]);
void inrtia(double x[3], double wt, double ten[3][3]);
double eigen1(double a[3][3], double v[3], double b[3][3]);
double eigen2(double a[3][3], double v[3], double b[3][3]); //does this actually get used
double inver (double a[3][3], double b[3][3]);
void unitry(double a[3][3], double b[3][3]);
void trnspz(double a[3][3], double b[3][3]);
void bij_uij(double bij[6], double rms[3], double vib[3][3],double amat[3][3]);//does this get used
void eigen(double w[3][3],double valu[3],double vect[3][3]);
double amax1(double a,double b, double c);
void colunit(double v[3],double a[3][3],int col);
double vmv(double a[3][3],int col1, double w[3][3], double b[3][3], int col2);
void axeqb(double a1[3][3], double v1[3], double v2[3], int col);
double amax2(double a, double b);
void zerovec(double x[3]);
void vecnorm(double x[3]);

void fconmul(double* a, double factor, double* b, double a_len);
int idot(int x[3], int y[3]);
void ivecsum(int x[3], int y[3], int z[3]);

#endif // VECTOR_H
