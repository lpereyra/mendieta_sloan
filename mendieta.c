#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesloan.h"
#include "timer.h"
#include "iden.h"
#include "colores.h"
#include "grid.h"

<<<<<<< HEAD
#define CVEL  299792.458
#define GCONST 6.6726E-11 /*M³s¯²kg¯¹*/
#define VUNIT 1000.
#define RUNIT 3.0856E+22
#define MSUN 1.989E+30  /*Masa del Sol en Kg*/

double mediana(double *x, int n, int n2);
int fsort(const void *a, const void *b);
=======
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
void Write_Groups(double *fof);
void Write_Sloan(double *fof);

int main(int argc, char **argv)
{
  int    i;
  double start,end;

  TIMER(start);
  
  init_variables(argc,argv);

  omp_set_nested(1);
  
  readsloan();

  // IMPORTANTE
  // Cambia origen de coordenadas
  change_positions(cp.npart);

  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmabmin, rmapmin);
  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmablim, rmaplim);
  fflush(stdout);

  for(i=0;i<nfrac;i++)
  {
    fprintf(stdout,"\nBegins Identification : Step %d of %d \n",i+1,nfrac);
   
    #ifdef LEN_FOF_MERCHAN
      iden.d0 = 3.0/(4.0*M_PI*(fof[i]+1)*(0.4*log(10.0)*flfia*rintlim));
    #else
      iden.d0 = 3.0*cp.vol/(4.0*M_PI*(cp.npart*(fof[i]+1.0))); // cbrt(vol/(npart*(fof+1)))
    #endif
   
    iden.d0 = cbrt(iden.d0);
    iden.step = i;
    iden.nobj = cp.npart;

    fprintf(stdout,"Linking length = %f \n",iden.d0);

    identification();

    Write_Groups(fof);
    
    free(Temp.head);
    free(Temp.npgrup);
    free(Temp.ll);
  }

  /************* TERMINO LA IDENTIFICACION ***************/
  Write_Sloan(fof);

  free(P);
  free(gal);
  grid_free();
  free(fof);
  for(i=0;i<NTABLA;i++)
    free(rinttabla[i]);
  free(rinttabla);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

void Write_Groups(double *fof)
{
  int i,j,k,id,npar,gn,save_sub;
  float dx,dy,dz,xc,yc,zc;
  char filename[200];
  FILE *pfout, *pfcentros, *pfcentros_ascii;

<<<<<<< HEAD
  FILE *pfgrupos;
  int   indmagmax;
  double *alfag, *deltag, *velg, *alfagx, *deltagx, *velgx;
  double magmax, magabi;
  double alfamed,deltamed,velmed;
  double alfamad,deltamad,velmad;
  double sum1alfa, sum1delta, sum1vel;
  double sum2alfa, sum2delta, sum2vel;
  double tmp, sum1sigma, sum2sigma, gi, wi;
  double rv, rvt, a1, a2, d1, d2, titaij, masa, sigma;
  double sigmaerr=30.0;
  double sqrtpi = sqrt(M_PI);

=======
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
  //////////////////////////////////////////////////////
  if(iden.step!=0)
  {
    #ifdef LIM_VOL
      sprintf(filename,"%.2f_%.2f_%.2f_fof.bin",zcut,fof[1],fof[0]);
    #else
      sprintf(filename,"%.2f_%.2f_fof.bin",fof[1],fof[0]);
    #endif
    pfout=fopen(filename,"w");
    i = iden.ngrupos-1;
    fwrite(&i,sizeof(int),1,pfout);

    #ifdef LIM_VOL
      sprintf(filename,"%.2f_%.2f_%.2f_centros.bin",zcut,fof[1],fof[0]);
    #else
      sprintf(filename,"%.2f_%.2f_centros.bin",fof[1],fof[0]);
    #endif
    pfcentros=fopen(filename,"w");
    fwrite(&i,sizeof(int),1,pfcentros);

    #ifdef LIM_VOL
      sprintf(filename,"%.2f_%.2f_%.2f_centros.dat",zcut,fof[1],fof[0]);
    #else
      sprintf(filename,"%.2f_%.2f_centros.dat",fof[1],fof[0]);
    #endif
    pfcentros_ascii=fopen(filename,"w");

<<<<<<< HEAD
    #ifdef LIM_VOL
      sprintf(filename,"%.2f_%.2f_%.2f_grupos.dat",zcut,fof[1],fof[0]);
    #else
      sprintf(filename,"%.2f_%.2f_grupos.dat",fof[1],fof[0]);
    #endif
    pfgrupos=fopen(filename,"w");

=======
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
  }else{
    #ifdef LIM_VOL
      sprintf(filename,"%.2f_%.2f_%.2f_fof.bin",zcut,fof[0],fof[1]);
    #else
      sprintf(filename,"%.2f_%.2f_fof.bin",fof[0],fof[1]);
    #endif
    pfout=fopen(filename,"w");
    i = iden.ngrupos-1;
    fwrite(&i,sizeof(int),1,pfout);
<<<<<<< HEAD
  }
=======
 }
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
  //////////////////////////////////////////////////////

  npar = gn = 0;

  for(i=1;i<iden.ngrupos;i++)
  {

    j = 0;
    id = k = Temp.head[i];
    if(iden.step!=0)
    {
<<<<<<< HEAD

      alfag   = (double *) malloc(Temp.npgrup[i]*sizeof(double));
      deltag  = (double *) malloc(Temp.npgrup[i]*sizeof(double));
      velg    = (double *) malloc(Temp.npgrup[i]*sizeof(double));
      alfagx  = (double *) malloc(Temp.npgrup[i]*sizeof(double));
      deltagx = (double *) malloc(Temp.npgrup[i]*sizeof(double));
      velgx   = (double *) malloc(Temp.npgrup[i]*sizeof(double));

      magmax = 0.0;

      xc = yc = zc = 0.0;
      save_sub = P[k].sub;
      fwrite(&save_sub,sizeof(int),1,pfout);

=======
      xc = yc = zc = 0.0;
      save_sub = P[k].sub;
      fwrite(&save_sub,sizeof(int),1,pfout);
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
    }
    fwrite(&i,sizeof(int),1,pfout);
    fwrite(&Temp.npgrup[i],sizeof(int),1,pfout);

    while(k != -1)
    {
      if(iden.step!=0)
      {
<<<<<<< HEAD

        alfag[j]  = gal[k].alfa           ;
        deltag[j] = gal[k].delta          ;
        velg[j]   = cp.h0*red2dis(gal[k].red);
        magabi    = gal[k].m[2]-25.0-5.0*log10(3000.0*gal[k].red);

        if(magabi<magmax)
        {
          magmax    = magabi;
          indmagmax = k;
        }

        //// cuidado con el orden {pos[i]-centro} en este caso
=======
        // cuidado con el orden {pos[i]-centro} en este caso
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
        dx = P[k].Pos[0] - P[id].Pos[0];
        dy = P[k].Pos[1] - P[id].Pos[1];
        dz = P[k].Pos[2] - P[id].Pos[2];

        #ifdef PERIODIC
        dx = dx > cp.lbox*0.5 ? dx-cp.lbox : dx;
        dy = dy > cp.lbox*0.5 ? dy-cp.lbox : dy;
        dz = dz > cp.lbox*0.5 ? dz-cp.lbox : dz;
  
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif

        xc += dx;
        yc += dy;
        zc += dz;
<<<<<<< HEAD


      }else{  

        P[k].sub = P[k].gr;

=======
      }else{   
        P[k].sub = P[k].gr;
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
      }
      //fwrite(&P[k].id,sizeof(int),1,pfout);
      fwrite(&k,sizeof(int),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    #ifdef DEBUG
    assert(j == Temp.npgrup[i]);
    #endif

    if(iden.step!=0)
    {
      xc /= (float)Temp.npgrup[i];
      yc /= (float)Temp.npgrup[i];
      zc /= (float)Temp.npgrup[i];

      xc += P[id].Pos[0];
      yc += P[id].Pos[1];
      zc += P[id].Pos[2];
      
      xc += pmin[0];
      yc += pmin[1];
      zc += pmin[2];
           
      #ifdef PERIODIC
      xc = xc<0 ? cp.lbox+(float)fmod(xc,cp.lbox) : (float)fmod(xc,cp.lbox);
      yc = yc<0 ? cp.lbox+(float)fmod(yc,cp.lbox) : (float)fmod(yc,cp.lbox);
      zc = zc<0 ? cp.lbox+(float)fmod(zc,cp.lbox) : (float)fmod(zc,cp.lbox);
      #endif

      fwrite(&save_sub,sizeof(int),1,pfcentros);
      fwrite(&i,sizeof(int),1,pfcentros);
      fwrite(&xc,sizeof(float),1,pfcentros);
      fwrite(&yc,sizeof(float),1,pfcentros);
      fwrite(&zc,sizeof(float),1,pfcentros);
      fwrite(&Temp.npgrup[i],sizeof(int),1,pfcentros);
      fprintf(pfcentros_ascii,"%d %d %f %f %f %d\n",save_sub,i,xc,yc,zc,Temp.npgrup[i]);
<<<<<<< HEAD
      
      qsort(alfag,Temp.npgrup[i],sizeof(double),fsort);
      qsort(deltag,Temp.npgrup[i],sizeof(double),fsort);
      qsort(velg,Temp.npgrup[i],sizeof(double),fsort);
      j = Temp.npgrup[i]/2;
    
      ////////// MEDIANAS /////////////
      if(Temp.npgrup[i]%2 == 0)
      {
        alfamed  = (alfag[j] + alfag[j - 1])*0.5;
        deltamed = (deltag[j] + deltag[j - 1])*0.5;
        velmed   = (velg[j] + velg[j - 1])*0.5;
      }else{
        alfamed  = alfag[j];
        deltamed = deltag[j];
        velmed   = velg[j]; 
      }

      /// CALCULA LAS MAD ///
      for(k=0; k<Temp.npgrup[i]; k++)
      {
        alfagx[k]  = fabs(alfag[k]-alfamed)   ;
        deltagx[k] = fabs(deltag[k]-deltamed) ;
        velgx[k]   = fabs(velg[k]-velmed)     ;
      }

      alfamad  = mediana(alfagx,Temp.npgrup[i],j);
      deltamad = mediana(deltagx,Temp.npgrup[i],j);
      velmad   = mediana(velgx,Temp.npgrup[i],j); 

      /// DEFINICION DE LOS MUI (reuso las gx) ///
      for(k=0;k<Temp.npgrup[i];k++)
      {
        if(alfamad>0)
          alfagx[k]  = (alfag[k]-alfamed)/(6.0*alfamad);
        else     
          alfagx[k]  = 1.0;

        if(deltamad>0)
          deltagx[k] = (deltag[k]-deltamed)/(6.0*deltamad);
        else
          deltagx[k] = 1.0 ;

        if(velmad>0)
          velgx[k]   = (velg[k]-velmed)/(6.0*velmad);
        else
          velgx[k] = 1.0 ;
      }

      //// CALCULO DE LAS SUMAS PARA LA MEDIA ////
      sum1alfa  = sum2alfa  = 0.0;
      sum1delta = sum2delta = 0.0;
      sum1vel   = sum2vel   = 0.0;

      for(k=0;k<Temp.npgrup[i];k++)
      {
        if(fabs(alfagx[k]) < 1.0)
        {
          tmp       = pow(1.0-alfagx[k]*alfagx[k],2) ;
          sum1alfa += (alfag[k]-alfamed)*tmp         ;
          sum2alfa += tmp                            ;
        }

        if(fabs(deltagx[k]) < 1.0)
        {
          tmp        = pow(1.0-deltagx[k]*deltagx[k],2) ;
          sum1delta += (deltag[k]-deltamed)*tmp         ;
          sum2delta += tmp                              ;
        }

        if(fabs(velgx[k]) < 1.0)
        {
          tmp      = pow(1.0-velgx[k]*velgx[k],2) ;
          sum1vel += (velg[k]-velmed)*tmp         ;
          sum2vel += tmp                          ;
        }
      }

      //// CALCULO DE LAS MEDIAS BIWEIGHTED (reuso las medianas) ////
      if(sum2alfa>0)
        alfamed  += sum1alfa/sum2alfa   ;
      if(sum2delta>0)
        deltamed += sum1delta/sum2delta ;
      if(sum2vel>0)
        velmed   += sum1vel/sum2vel     ;

      /// CALCULO DE LAS DISPERSIONES ////
      if(Temp.npgrup[i]<15)
      {
        //// GAPPER ////
        sum1sigma = 0.0;
        for(k=0;k<(Temp.npgrup[k]-1);k++)
        {
          gi = velg[k+1]-velg[k];
          wi = ((double)k+1.0)*((double)Temp.npgrup[k]-(double)k-1.0);
          sum1sigma+=wi*gi;
        }

        sigma=(sqrtpi/((double)Temp.npgrup[i]*((double)Temp.npgrup[i]-1.0)))*sum1sigma;

      }else{

        //// BIWEIGHTED ////
        sum1sigma = sum2sigma = 0.0;

        for(k=0;k<Temp.npgrup[i];k++)
        {
          if(fabs(velgx[k]) < 1.0 )
          {
            tmp        = 1.0-velgx[k]*velgx[k];
            sum1sigma += pow((velg[k]-velmed),2)*pow(tmp,4);
            sum2sigma += tmp*(1.0-5.0*velgx[k]*velgx[k]);
          }
        }

        sigma = sqrt((double)Temp.npgrup[i])*sqrt(sum1sigma)/fabs(sum2sigma);
      }

      sigma = sigma*sigma*(double)Temp.npgrup[i]/((double)Temp.npgrup[i]-1.0)-sigmaerr*sigmaerr;
      sigma = sigma>0 ? sqrt(sigma) : 0.0;

      //// RADIO VIRIAL Y MASA ////
      rv=0.0;
      for(k=0;k<Temp.npgrup[i]-1;k++)
      {
        a1 = alfag[k];
        d1 = deltag[k];
        gi = sin(d1);
        wi = cos(d1);

        for(j=k+1;j<Temp.npgrup[i];j++)
        {
          a2 = alfag[j];
          d2 = deltag[j];

          titaij = gi*sin(d2)+wi*cos(d2)*cos(a1-a2);
          
          if(fabs(titaij)>1.0) titaij = 1.0;
          titaij = acos(titaij);
          rv += 1.0/(2.0*(velmed/cp.h0)*sin(titaij/2.0));
          //rv+=1.0/(titaij*velmed/h0);
        }
      }

      rv = ((double)(Temp.npgrup[i]*Temp.npgrup[i]))/rv;

      tmp = sigma*VUNIT;
      rvt = rv*3.1415926/2.0*RUNIT;
      masa = (3.0*tmp*tmp)*rvt/GCONST/MSUN;

      fprintf(pfgrupos,"%9.4f %9.4f %11.4f %9.4f %9.4f %6d %7d %9.4e %9.4f %7d\n",
      alfamed,deltamed,velmed,sigma,rv,Temp.npgrup[i],i,masa,magmax,indmagmax);
      //if(masa<=0.0 && Temp.npgrup[i] > 3 ) printf("%d %d %e %e\n",i,Temp.npgrup[i],sigma,rv);
      //fflush(stdout);
  
      free(alfag);
      free(deltag);
      free(velg);
      free(alfagx);
      free(deltagx);
      free(velgx);
    }  

    npar+=Temp.npgrup[i];
=======
    }  

    npar+=j;
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
    gn++;
  }

  assert(gn == (iden.ngrupos-1));
  fclose(pfout);

  if(iden.step!=0)
  {
    fclose(pfcentros);
    fclose(pfcentros_ascii);
  }

  fprintf(stdout,"num de grupos %d num de particulas en grupos %d\n",gn,npar);
  fflush(stdout);

  return;
}

void Write_Sloan(double *fof)
{
  int i;
<<<<<<< HEAD
  char filename[200], filename_ascii[200];
  FILE *pfout, *pfout_ascii;
=======
  char filename[200];
  FILE *pfout;
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
  int **npgrupo;
  struct galsloang galg;

  npgrupo = (int **) calloc(cp.npart,sizeof(int *));
  for(i=0;i<cp.npart;i++)
    npgrupo[i] = (int *) calloc(nfrac,sizeof(int));

  for(i=0;i<cp.npart;i++)
  {
    npgrupo[P[i].sub][0]++;
    npgrupo[P[i].gr][1]++;
  }

  #ifdef LIM_VOL
<<<<<<< HEAD
    sprintf(filename,"cat_%.2f_%.2f_%.2f.bin",zcut,fof[0],fof[1]);
    sprintf(filename_ascii,"cat_%.2f_%.2f_%.2f.dat",zcut,fof[0],fof[1]);
  #else
    sprintf(filename,"cat_%.2f_%.2f.bin",fof[0],fof[1]);
    sprintf(filename_ascii,"cat_%.2f_%.2f.dat",fof[0],fof[1]);
  #endif
  pfout=fopen(filename,"w");
  pfout_ascii=fopen(filename_ascii,"w");
=======
    sprintf(filename,"cat_%.2f_%.2f_%.2f.dat",zcut,fof[0],fof[1]);
  #else
    sprintf(filename,"cat_%.2f_%.2f.dat",fof[0],fof[1]);
  #endif
  pfout=fopen(filename,"w");
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b

  fwrite(&cp.npart,sizeof(int),1,pfout);
  for(i=0;i<cp.npart;i++)
  {
    galg.gal       = gal[i];
    galg.gal.alfa  /= PI180;
    galg.gal.delta /= PI180;
    galg.grupo[0]   = P[i].sub;
    galg.grupo[1]   = P[i].gr;
    galg.npgrupo[0] = npgrupo[P[i].sub][0];
    galg.npgrupo[1] = npgrupo[P[i].gr][1];
    fwrite(&galg,sizeof(struct galsloang),1,pfout);
<<<<<<< HEAD
    
    fprintf(pfout_ascii,\
    "%ld %ld %ld \
%f %f %f %f \
%f %f %f %f %f \
%f %f %f %f %f \
%f %f %f %f %f \
%f %f %f %f %f \
%f %f %f %f %f \
%f %f %f \
%f %f %f \
%f %f %f \
%f %f %f \
%f %f %f %f %f \
%d %d %d %d\n",
    galg.gal.id, galg.gal.targetid, galg.gal.specid, \
    galg.gal.alfa, galg.gal.delta, galg.gal.red, galg.gal.rederr, \
    galg.gal.m[0], galg.gal.m[1], galg.gal.m[2], galg.gal.m[3], galg.gal.m[4], \
    galg.gal.merr[0], galg.gal.merr[1], galg.gal.merr[2], galg.gal.merr[3], galg.gal.merr[4], \
    galg.gal.modelm[0], galg.gal.modelm[1], galg.gal.modelm[2], galg.gal.modelm[3], galg.gal.modelm[4],
    galg.gal.modelmerr[0], galg.gal.modelmerr[1], galg.gal.modelmerr[2], galg.gal.modelmerr[3], galg.gal.modelmerr[4], \
    galg.gal.ext[0], galg.gal.ext[1], galg.gal.ext[2], galg.gal.ext[3], galg.gal.ext[4], \
    galg.gal.petrosg[0], galg.gal.petrosg[1], galg.gal.petrosg[2], \
    galg.gal.petrosgerr[0], galg.gal.petrosgerr[1], galg.gal.petrosgerr[2], \
    galg.gal.petrosr[0], galg.gal.petrosr[1], galg.gal.petrosr[2], \
    galg.gal.petrosrerr[0], galg.gal.petrosrerr[1], galg.gal.petrosrerr[2], \
    galg.gal.k[0], galg.gal.k[1], galg.gal.k[2], galg.gal.k[3], galg.gal.k[4], \
    galg.grupo[0], galg.grupo[1], galg.npgrupo[0], galg.npgrupo[1]);

  } 

  fclose(pfout);
  fclose(pfout_ascii);

=======
  } 

>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
  for(i=0;i<cp.npart;i++)
    free(npgrupo[i]);
  free(npgrupo);

  return;
}
<<<<<<< HEAD

double mediana(double *x, int n, int n2)
{
  double *xx;
  double m;

  xx = malloc(n*sizeof(double));
  memcpy(xx,x,n*sizeof(double));
  qsort(xx,n,sizeof(double),fsort);

  if((n%2)==0)  
    m = 0.5*(xx[n2]+xx[n2-1]);
  else
    m = (xx[n2]);

  free(xx);

  return m;
}

int fsort(const void *a, const void *b)
{
  if((double *) a < (double *) b)
    return -1;

  if((double *) a > (double *) b)
    return +1;

  return 0;
}
=======
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b
