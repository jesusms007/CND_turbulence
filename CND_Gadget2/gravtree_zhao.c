#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;

  /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
#endif
	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[i].Type == 0)
		      GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;
		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;
		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);

      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] *= All.G;


  /* Finally, the following factor allows a computation of a cosmological simulation
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif


//External potential!!!!!!!!!!

#ifdef GRAVITY_GALAXY /////SALAS MOD
  double r=0,v=0;
  double x=0,y=0,z=0,vx=0,vy=0,vz=0;
  double phi=0;
  double r0=100.0,rho0=0.0004,a=0.25,b22=0.1,b20=0.3; //for the bar potential ... Salas mod Oct 25, 2018 == chnged units
 // double omega=63.; //bar pattern speed from Kim et al. 2011, 63 km s^-1 kpc^-1
  double mult=1;
  double barx=0,bary=0,barz=0;
  double cenx=0,ceny=0,cenz=0;
  double corx=0,cory=0,corz=0;
  ///double fric=0,Msat=2e-4,Mgal=1.26e-2,chi=0,pi=3.141592,rho=0;
  double i1=0,k0=0,k1=0;//Bessel functions
//  double gamma=0.577216;//Euler's constant
  for(i = 0; i < NumPart; i++)
    {
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	x = P[i].Pos[0];
    y = P[i].Pos[1];
    z = P[i].Pos[2];
	vx = P[i].Vel[0];
    vy = P[i].Vel[1];
    vz = P[i].Vel[2];
	phi = atan2(y,x);
          if (i == 100 && ThisTask == 0)
          {  printf("--Calculating GC gravity Zhao+ 1991--\n"); ///SALAS MOD OCT 25, 2018
          }
	r = sqrt(x*x + y*y + z*z);
	if(r > 0)
	  {
	    r = sqrt(r*r + 0.05*0.05);

	    if(r>1.0)
          {
	    mult = -6.2832*pow((r/r0),a)*rho0*r0*r0*All.G/((3.+a)*(-2.+a)*(1.+a)*pow(r,4));

	    barx = mult*(-6.*b20*z*z*x-x*b20*y*y*a-x*b20*y*y*a*a-4.*x*b20*z*z*a+2.*x*b20*z*z*a*a+12.*x*x*x+12.*z*z*x-2.*x*y*y*a
		-2.*x*y*y*a*a-2.*x*z*z*a-2.*x*z*z*a*a-b20*x*x*x*a-b20*x*x*x*a*a-2.*x*x*x*a-2.*x*x*x*a*a+12.*x*y*y
		+12.*b22*z*z*x*cos(2.*phi)+12.*b22*sin(2.*phi)*y*x*x
		+12.*b22*sin(2.*phi)*y*z*z+6.*b22*cos(2.*phi)*x*x*x*a
		+6.*b22*cos(2.*phi)*x*x*x*a*a+12.*b22*sin(2.*phi)*y*y*y*a
		+12.*b22*sin(2.*phi)*y*y*y+6.*x*b22*cos(2.*phi)*y*y*a
		+6.*x*b22*cos(2.*phi)*y*y*a*a+12.*b22*z*z*x*cos(2.*phi)*a
		+12.*b22*sin(2.*phi)*y*x*x*a+12.*b22*sin(2.*phi)*y*z*z*a);

	    bary = mult*(-6.*b20*z*z*y-y*b20*x*x*a-y*b20*x*x*a*a-4.*y*b20*z*z*a+2.*y*b20*z*z*a*a+12.*y*y*y+12.*z*z*y-2.*y*x*x*a
		-2.*y*x*x*a*a-2.*y*z*z*a-2.*y*z*z*a*a-b20*y*y*y*a-b20*y*y*y*a*a-2.*y*y*y*a-2.*y*y*y*a*a+12.*y*x*x
		-12.*b22*sin(2.*phi)*x*x*x+6.*y*b22*cos(2.*phi)*x*x*a
		+6.*y*b22*cos(2.*phi)*x*x*a*a+12.*b22*z*z*y*cos(2.*phi)*a
		-12.*b22*sin(2.*phi)*x*y*y*a-12.*b22*sin(2.*phi)*x*z*z*a
		+12.*b22*z*z*y*cos(2.*phi)-12.*b22*sin(2.*phi)*x*y*y
		-12.*b22*sin(2.*phi)*x*z*z+6.*b22*cos(2.*phi)*y*y*y*a
		+6.*b22*cos(2.*phi)*y*y*y*a*a-12.*b22*sin(2.*phi)*x*x*x*a);

	    barz = mult*z*(-6.*b22*cos(2.*phi)*x*x*a-6.*b22*cos(2.*phi)*y*y*a
		+6.*b22*cos(2.*phi)*x*x*a*a+6.*b22*cos(2.*phi)*y*y*a*a+12.*x*x+12.*y*y+12.*z*z
		-12.*b22*cos(2.*phi)*x*x-12.*b22*cos(2.*phi)*y*y+5.*b20*x*x*a+5.*b20*y*y*a
		+2.*b20*z*z*a-b20*x*x*a*a-b20*y*y*a*a+2.*b20*z*z*a*a+6.*b20*x*x+6.*b20*y*y-2.*x*x*a-2.*y*y*a-2.*z*z*a-2.*x*x*a*a
		-2.*y*y*a*a-2.*z*z*a*a);

	    corx = 2.*All.OmegaBar*vy;
        cory = -2.*All.OmegaBar*vx;
        corz = 0;

	    cenx = All.OmegaBar*All.OmegaBar*x;
        ceny = All.OmegaBar*All.OmegaBar*y;
        cenz = All.OmegaBar*All.OmegaBar*z;

         #ifdef ROTATING_BAR
        if (i == 100 && ThisTask == 0)
        {  printf("---Adding bar rotation = %g km/(s kpc)\n", All.OmegaBar * kpc_to_km / All.UnitTime_in_s);
        }

	    P[i].GravAccel[0] += -barx+corx+cenx;
	    P[i].GravAccel[1] += -bary+cory+ceny;
	    P[i].GravAccel[2] += -barz+corz+cenz;
          #endif
          }
	    //BH contribution 4e6 Msun
	    P[i].GravAccel[0] += -All.G*4.0*x/r/r/r;
	    P[i].GravAccel[1] += -All.G*4.0*y/r/r/r; // mass unit = 1e6 msun
	    P[i].GravAccel[2] += -All.G*4.0*z/r/r/r;

	    } ///end of if r(>0)

	  } ///end of if(P[i].Ti_endstep == All.Ti_Current)

      }
#endif



  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
        All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
        All.SofteningTable[0] = All.SofteningGas;

      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
        All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
        All.SofteningTable[1] = All.SofteningHalo;

      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
        All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
        All.SofteningTable[2] = All.SofteningDisk;

      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
        All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
        All.SofteningTable[3] = All.SofteningBulge;

      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
        All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
        All.SofteningTable[4] = All.SofteningStars;

      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
        All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
        All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
