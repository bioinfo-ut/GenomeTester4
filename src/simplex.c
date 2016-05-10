#define __KHAYYAM_SIMPLEX_C__

//
// Simplex parameter fitting
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simplex.h"

float
downhill_simplex (int NDIM, float MX[], float MdX[], float EMax, int nruns, int niterations, float (*func) (int, const float[], void *), void *data)
{
	double MP[26][25];          /* Main matrix of simplex vertices         */
	double Pb[25];              /* The point Pb.                           */
	double Pr[25];              /* The point Pr.                           */
	double Prr[25];             /* The point Prr or Prr'.                  */
	double Y[26];               /* Results Y[i]=FUNK(MP[i])                */
	double RTol;                /* Real toleration                         */
	double Al, Bt, Gm, Ypr, Yprr, Yavr, Rmp;
	double xa, xb, xc, xd, lMin;
	int i, j, MPTS;
	int ITR0, MITR, ITR1, NITR;
	int iLo, iHi, iNHi;			/* Lowest point, Highest point, Next highest point */
	double YLo, YHi, DLo, DHi;

	MPTS = NDIM + 1;
	Al = 1.0;
	Bt = 0.5;
	Gm = 2.0;
	Rmp = MPTS;
	MITR = niterations;
	NITR = nruns;

	Y[0] = func (NDIM, MX, data);
	for (ITR1 = 0; ITR1 < NITR; ITR1++) {
		ITR0 = 0;
		srand (0);
		for (i = 0; i < NDIM; i++) {
			for (j = 0; j < MPTS; j++) MP[j][i] = MX[i];
			MP[i][i] += MdX[i] * (0.9 + 0.2 * rand () / RAND_MAX) / (5 * ITR1 + 1);
			// MP[i][i] += MdX[i];
			// MdX[i] /= 2;
		}
		/* Rotarion of the matrix MP - not implememnted */
		for (j = 0; j <= NDIM; j++) {
			for (i = 0; i < NDIM; i++) MX[i] = (float) MP[j][i];
			Y[j] = func (NDIM, MX, data);
		}
		YLo = Y[0];
		YHi = Y[0];
		while (ITR0 < MITR) {
			Yavr = 0;
			iLo = 0;
			if (Y[0] > Y[1]) {
				iHi = 0;
				iNHi = 1;
			} else {
				iHi = 1;
				iNHi = 0;
			}
			for (i = 0; i < MPTS; i++) {
				Yavr += Y[i];
				if (Y[i] < Y[iLo]) iLo = i;
				if (Y[i] > Y[iHi]) {
					iNHi = iHi;
					iHi = i;
				} else if (Y[i] > Y[iNHi]) {
					if (i != iHi) iNHi = i;
				}
			}
			Yavr /= Rmp;
			RTol = 2.0 * fabs (Y[iHi] - Y[iLo]) / (fabs (Y[iHi]) + fabs (Y[iLo]));

			DLo = fabs (YLo - Y[iLo]);
			DHi = fabs (YHi - Y[iHi]);
			YLo = Y[iLo];
			YHi = Y[iHi];
			// fprintf (stderr, "%d %d YLo %g YHi %g DLo %g DHi %g\n", ITR1, ITR0, YLo, YHi, DLo, DHi);
			// if (YLo <= EMax) break;

			ITR0++;
			// printf(">%3d%4d RTol=%e Ymin=%f Ymax=%f Yavr=%f\n",ITR1,ITR0,RTol,Y[iLo],Y[iHi],Yavr);
			for (j = 0; j < NDIM; j++) Pb[j] = 0.0;
			for (i = 0; i < MPTS; i++) {
				if (i != iHi) {
					for (j = 0; j < NDIM; j++) Pb[j] += MP[i][j];
				}
			}
			for (j = 0; j < NDIM; j++) {
				Pb[j] /= NDIM;
				Pr[j] = (1.0 + Al)*Pb[j] - Al*MP[iHi][j];
			}
			for (j = 0; j < NDIM; j++) MX[j] = (float) Pr[j];
			Ypr = func (NDIM, MX, data);
			if (Ypr <= Y[iLo]) {
				for (j = 0; j < NDIM; j++) Prr[j] = Gm*Pr[j] + (1.0 - Gm)*Pb[j];
				for (j = 0; j < NDIM; j++) MX[j] = (float) Prr[j];
				Yprr = func (NDIM, MX, data);
				if (Ypr > Yprr) {
					for (j = 0; j < NDIM; j++) MP[iHi][j] = Prr[j];
					Y[iHi] = Yprr;
				} else {
					for (j = 0; j < NDIM; j++) MP[iHi][j] = Pr[j];
					Y[iHi] = Ypr;
				}
			} else {
				if (Ypr >= Y[iNHi]) {
					if (Ypr < Y[iHi]) {
						for (j = 0; j < NDIM; j++) MP[iHi][j] = Pr[j];
						Y[iHi] = Ypr;
					}
					for (j = 0; j < NDIM; j++) Prr[j] = Bt*MP[iHi][j] + (1.0 - Bt)*Pb[j];
					for (j = 0; j < NDIM; j++) MX[j] = (float) Prr[j];
					Yprr = func (NDIM, MX, data);
					if (Yprr < Y[iHi]) {
						for (j = 0; j < NDIM; j++) MP[iHi][j] = Prr[j];
						Y[iHi] = Yprr;
					} else {
						/*      for(i=0;i<MPTS;i++) if(i!=iLo) {
								  for(j=0;j<NDIM;j++) {Pr[j]=0.5*(MP[i][j]+MP[iLo][j]);MP[i][j]=Pr[j];}
								  for(j=0;j<NDIM;j++) *(MX[j])=Pr[j];Y[i]=func();}}}                     */
						for (j = 0; j < NDIM; j++) Pr[j] = 0.5*(MP[iHi][j] + MP[iLo][j]);
						for (j = 0; j < NDIM; j++) MX[j] = (float) Pr[j];
						Ypr = func (NDIM, MX, data);
						if (Ypr < Y[iHi]) {
							for (j = 0; j < NDIM; j++) MP[iHi][j] = Pr[j];
							Y[iHi] = Ypr;
						} else {
							for (j = 0; j < NDIM; j++) Prr[j] = -MP[iHi][j] + 2.0*MP[iLo][j];
							for (j = 0; j < NDIM; j++) MX[j] = (float) Prr[j];
							Yprr = func (NDIM, MX, data);
							if (Yprr < Y[iHi]) {
								for (j = 0; j < NDIM; j++) MP[iHi][j] = Prr[j];
								Y[iHi] = Yprr;
							} else {
								xa = 3 * Y[iHi] - 8 * Ypr + 6 * Y[iLo] - Yprr;
								xb = Y[iHi] - 2 * Y[iLo] + Yprr;
								xc = -0.5*Y[iHi] + 8 * Ypr / 3 - 2 * Y[iLo] + Yprr / 6;
								xd = xb*xb - 4 * xa*xc;
								if (xd > 0) {
									lMin = 0.5*(-xb - sqrt (xd)) / xa;
									for (j = 0; j < NDIM; j++) Pr[j] = lMin*MP[iHi][j] + (1 - lMin)*MP[iLo][j];
									for (j = 0; j < NDIM; j++) MX[j] = (float) Pr[j];
									Ypr = func (NDIM, MX, data);
								}
								if (Ypr < Y[iHi]) {
									for (j = 0; j < NDIM; j++) MP[iHi][j] = Pr[j];
									Y[iHi] = Ypr;
								} else {
									for (j = 0; j < NDIM; j++) MP[iHi][j] = MP[iLo][j];
									Y[iHi] = Y[iLo];
								}
							}
						}
					}
				} else {
					for (j = 0; j < NDIM; j++) MP[iHi][j] = Pr[j];
					Y[iHi] = Ypr;
				}
			}
		}
		iLo = 0;
		for (i = 1; i < MPTS; i++) if (Y[i] < Y[iLo]) iLo = i;
		for (i = 0; i < NDIM; i++) MX[i] = (float) MP[iLo][i];
		// fprintf (stderr, "%d %d YLo %g YHi %g\n", ITR1, ITR0, YLo, YHi);
		// for (i = 0; i < NDIM; i++) MdX[i] /= 4;
	}
	return (float) Y[iLo];
}
