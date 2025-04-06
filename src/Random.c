#include "stdafx.h"
//#include "Mesosim Resources.h"
#include "Defs.h"
#include "Geometry.h"
#include "Prototypes.h"
#include "Vector.h"
#include "Global_Externs.h"
#include "Simulation_Global_Externs.h"
// TODO: remove most of these

//this file is the exact same as Random Number.cpp except for the obvious change in header

/******************************************************************************/
/******************************************************************************/

double drandj(long *idum)
	{
		int j;
		long k;
		double temp;

		k = *idum/IQ1;
		*idum = IA1*(*idum-k*IQ1) - k*IR1;
		if (*idum < 0) *idum += IM1;
		k = idum2/IQ2;
		idum2 = IA2*(idum2-k*IQ2)-k*IR2;
		if (idum2 < 0) idum2 += IM2;
		j = iy/NDIV;
		iy = iv[j] - idum2;
		iv[j] = *idum;
		if (iy < 1) iy += IMM1;
		temp = AM*iy;
		if (temp > RNMX) return RNMX;
		else return temp;
	}

/******************************************************************************/
/******************************************************************************/

void srandj(long *idum)
	{
		int j;
		long k;

		if (*idum < 0) *idum = -(*idum);

		if (*idum < 1) *idum = 1;
		idum2 = *idum;

		for (j = NTAB+7; j>=0; j--)
			{
				k = (*idum)/IQ1;
				*idum = IA1*((*idum)-k*IQ1) - k*IR1;
				if (*idum < 0) *idum += IM1;
				if (j < NTAB) iv[j] = *idum;
			}
		iy = iv[0];

	}


