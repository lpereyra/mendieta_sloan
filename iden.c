#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "iden.h"
#include "bitmask.h"
#include "leesloan.h"

void identification(void)
{
  int i,j,l,tid,nvec;
  int ngrid_old = grid.ngrid;
  int *test;
  double r0 = v0;
  int *list;
  int *listnew;
  int *c;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif
 
  #ifdef LOCK
  omp_lock_t *lock; 
  lock = (omp_lock_t *) malloc(cp.npart*sizeof(omp_lock_t));
  for(i=0;i<cp.npart;i++) omp_init_lock(&(lock[i]));
  #endif

  if(iden.d0>r0)
  {
    r0=iden.d0;
    fprintf(stdout,"d0 > v0, %f > %f \n",iden.d0,v0);
  }else{
    fprintf(stdout,"v0 > d0, %f > %f \n",v0,iden.d0); 
  }

  printf("lado    %f\n",cp.lbox);
  printf("r0      %f\n",r0);
  printf("rmaplim %f\n",rmaplim);
  printf("dismax  %f\n",cp.dlummax);
  printf("rintlim %f\n",rintlim);
  

  #ifdef LIM_VOL
    grid.ngrid = (int)(cp.lbox/(r0*pow(intfl(MAGMENOSINF,rmaplim-25.0-5.0*log10(cp.dlummax))/rintlim ,-1./3.)));
  #else
    grid.ngrid = (int)(cp.lbox/(r0*pow(intfl(rmapmin-25.0-5.0*log10(cp.dlummax),rmaplim-25.0-5.0*log10(cp.dlummax))/rintlim,-1./3.)));
  #endif

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  
  grid.nobj = iden.nobj;
  grid.step = iden.step;

  if(ngrid_old != grid.ngrid || iden.step==1)
  {
    grid_free();
    grid_init();
    grid_build();
  }

  test = (int *) calloc(cp.npart/32 + 1,sizeof(int));
  list = (int *) malloc(cp.npart*sizeof(int));
  c    = (int *) malloc(NTHREADS*sizeof(int));

  fprintf(stdout,"Comienza identificacion.....%d\n",iden.step);

  printf("\nRunning on %d threads\n",NTHREADS);

  #pragma omp parallel default(none) private(tid,i) \
  shared(P,iden,cp,test,stdout)   
  {
    tid = omp_get_thread_num(); 
    
    for(i = tid*floor((float)cp.npart/NTHREADS);
    i<(tid==NTHREADS-1 ? cp.npart : (tid+1)*floor((float)cp.npart/NTHREADS));
    i++)
    {
      P[i].gr = i;
      if(iden.step!=0 && P[i].sub == 0)
        SetBit(test,i);     
    }
  }

  for(i=0;i<cp.npart;i++)
  {

    if(TestBit(test,i)) continue ;  // Salta a la siguiente

    fprintf(stdout,"centro %d redshift %f\r",i,gal[i].red);
    fflush(stdout);

    nvec = 0;
    SetBit(test,i);
    #ifdef LOCK
      busv(i,&nvec,list,test,lock);
    #else
      busv(i,&nvec,list,test);
    #endif

    if(nvec > 0)
    {

      do
      {

        #ifdef LOCK
          #pragma omp parallel default(none) private(tid,j,l,listnew) \
          shared(iden,cp,test,c,list,nvec,lock,stdout)   
        #else
          #pragma omp parallel default(none) private(tid,j,l,listnew) \
          shared(iden,cp,test,c,list,nvec,stdout)   
        #endif
        {
          tid = omp_get_thread_num(); 
          c[tid] = 0;
          listnew = (int *) malloc(cp.npart*sizeof(int));

          for(j = tid*floor((float)nvec/NTHREADS);
          j<(tid==NTHREADS-1 ? nvec : (tid+1)*floor((float)nvec/NTHREADS));
          j++)
          {

            if(TestBit(test,list[j])) continue ;  // Salta a la siguiente

            #ifdef LOCK
              
              omp_set_lock(&(lock[list[j]]));
              if(TestBit(test,list[j]))
              {               
                omp_unset_lock(&(lock[list[j]]));
                continue;  // Salta a la siguiente
              }else{
                SetBit(test,list[j]);
                omp_unset_lock(&(lock[list[j]]));
              }

              busv(list[j],&c[tid],listnew,test,lock);

            #else

              SetBit(test,list[j]);
              busv(list[j],&c[tid],listnew,test);

            #endif
          }    

          #pragma omp barrier

          l = 0;
          for(j=0;j<tid;j++)
            l += c[j];

          #pragma omp barrier

          for(j=0;j<c[tid];j++)
          {
            assert(listnew[j]>=0 && listnew[j]<cp.npart);
            assert((j+l)>=0 && (j+l)<cp.npart);
            list[j+l] = listnew[j];
          }

          free(listnew);

        }
        
        nvec = 0;
        for(j=0;j<NTHREADS;j++)
          nvec += c[j];

      }while( nvec != 0 );

    }

  }

  #ifdef LOCK
  for(i=0;i<cp.npart;i++) 
    omp_destroy_lock(&(lock[i]));
  free(lock);
  #endif

  linkedlist(test);
  free(test);
  free(list);
  free(c);

  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);

  return;

}

int Raiz(int i)
{

 if(i != P[i].gr)
   P[i].gr = Raiz(P[i].gr);

 return P[i].gr;

}

#ifdef LOCK
void Unir(int u, int v, omp_lock_t *lock)
#else
void Unir(int u, int v)
#endif
{
  int z;

  while(P[u].gr != P[v].gr)
  { 
      if(P[u].gr < P[v].gr)
      {
        #ifdef LOCK
          if(u == P[u].gr)
          {
            omp_set_lock(&(lock[u]));
            z = 0;
            if(u == P[u].gr)
            {
                P[u].gr = P[v].gr;  
                z = 1;
            }
            omp_unset_lock(&(lock[u]));
            if(z==1) break;             
          }
        #else
          if(u == P[u].gr)
          {
            P[u].gr = P[v].gr;  
            break;             
          }
        #endif
          
          z = P[u].gr;   
          P[u].gr = P[v].gr;
          u = z;

      }else{

        #ifdef LOCK
          if(v == P[v].gr)
          {
            omp_set_lock(&(lock[v]));
            z = 0;
            if(v == P[v].gr)
            {
                P[v].gr = P[u].gr;   
                z = 1;
            }
            omp_unset_lock(&(lock[v]));
            if(z == 1) break;            
          }
        #else
          if(v == P[v].gr)
          {
              P[v].gr = P[u].gr;   
              break;             
          }
        #endif

          z = P[v].gr;   
          P[v].gr = P[u].gr;   
          v = z;

      }
  }

  return;
}

#ifdef LOCK
  void busv(int ic, int *nvec, int *list, int *test, omp_lock_t *lock)
#else
  void busv(int ic, int *nvec, int *list, int *test)
#endif
{

  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  int  itabla,imintabla;
  float dl,dl1,dl2,d12;
  float rscale,vl,v12;
  float rm12,rmmin12;
  float rint12,rint121,rint122;
  float d1,d2,tita12;
  int i;
  float lbox,fac;
  float tridif,numera,alfa12;
  long ngrid;
  #ifdef PERIODIC
  float lbox2;
  #endif

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (float)ngrid/lbox;
  #ifdef PERIODIC
  lbox2 = lbox/2.0;
  #endif

  ixc  = (int)(P[ic].Pos[0]*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(P[ic].Pos[1]*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(P[ic].Pos[2]*fac);
  izci = izc - 1;
  izcf = izc + 1;

  d1   = P[ic].Dis;                // distancia cosmologica
  dl1  = d1*(1.0+gal[ic].red);     // distancia luminosidad

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif

      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {
          
          if(P[i].sub != P[ic].sub)
          {
            i = grid.ll[i];
            continue;
          }
 
          if(TestBit(test,i))
          {
            i = grid.ll[i];
            continue;
          }      

          // DISTANCIAS DE LA VECINA
          d2   = P[i].Dis;                    // distancia cosmologica
          dl2  = d2*(1.0+gal[i].red);         // distancia luminosidad 

          rm12    = rmaplim-25.0-5.0*log10((dl1+dl2)*0.5);
          rmmin12 = rmapmin-25.0-5.0*log10((dl1+dl2)*0.5);

          // CALCULA LA INTEGRAL DEL FACTOR DE ESCALA
          itabla=(int)((rm12-mt2min)/dmt2);
          imintabla=(int)((rmmin12-mt1min)/dmt1);

          if(itabla<0 || itabla >= NTABLA )
          {
            fprintf(stdout,"\nXXXX itabla  = %d\n ",itabla);
            exit(-33);
          }

          //assert(itabla>=0 && itabla<NTABLA);

          if(imintabla<0 || imintabla >= NTABLA )
          {
            fprintf(stdout,"\nXXXX imintabla= %d\n ",imintabla);
            exit(-33);
          }

          //assert(imintabla>=0 && imintabla<NTABLA);

          rint12  = (rinttabla[imintabla][itabla+1]-rinttabla[imintabla][itabla])/dmt2;
          rint121 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla][itabla];

          rint12  = (rinttabla[imintabla+1][itabla+1]-rinttabla[imintabla+1][itabla])/dmt2;
          rint122 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla+1][itabla];

          rint12 = (rint122-rint121)/dmt1;
          rint12 = rint12*(rmmin12-imintabla*dmt1-mt1min)+rint121;

          rscale = pow(rint12/rintlim,-0.33333);
          dl = iden.d0*rscale;
          vl = v0*rscale;

          alfa12 = gal[i].alfa-gal[ic].alfa;
          tita12 = sin(gal[i].delta)*sin(gal[ic].delta)+cos(gal[i].delta)*cos(gal[ic].delta)*cos(alfa12);

          numera = (cos(gal[ic].delta)*sin(alfa12));
          numera *= numera;

          tridif = cos(gal[i].delta)*sin(gal[ic].delta)-sin(gal[i].delta)*cos(gal[ic].delta)*cos(alfa12);
          tridif *= tridif;

          numera += tridif;
          numera = sqrt(numera);

          tita12 = atan2(numera,tita12);

          d12 = sin(tita12/2.0)*(d1+d2);
          v12 = fabs(d1-d2);

          //tita12=sin(gal[i].delta)*sin(gal[ic].delta)+cos(gal[i].delta)*cos(gal[ic].delta)*cos(gal[i].alfa-gal[ic].alfa);

          //if(fabs(tita12)>1.0){ printf("tita12 = %f\n",tita12); tita12=1.0;}

          //tita12=acos(tita12);
          //d12 = sin(tita12/2.0)*(d1+d2);
          //v12 = fabs(d1-d2);

          if((d12<dl) && (v12<vl))
          {
            #ifdef LOCK
              Unir(ic,i,lock);
            #else
              Unir(ic,i);
            #endif
            list[*nvec] = i;
            *nvec = *nvec + 1;
          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

}

#ifdef LOCK
void busv_rec(int ic, int *test, omp_lock_t *lock)
#else
void busv_rec(int ic, int *test)
#endif
{

  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  int  itabla,imintabla;
  float dl,dl1,dl2,d12;
  float rscale,vl,v12;
  float rm12,rmmin12;
  float rint12,rint121,rint122;
  float d1,d2,tita12;
  float tridif,numera,alfa12;
  int i;
  float lbox,fac;
  long ngrid;
  #ifdef PERIODIC
  float lbox2;
  #endif

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (float)ngrid/lbox;
  #ifdef PERIODIC
  lbox2 = lbox/2.0;
  #endif

  ixc  = (int)(P[ic].Pos[0]*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(P[ic].Pos[1]*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(P[ic].Pos[2]*fac);
  izci = izc - 1;
  izcf = izc + 1;
  d1   = P[ic].Dis;                // distancia cosmologica
  dl1  = d1*(1.0+gal[ic].red);     // distancia luminosidad

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
     if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif

      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {

          if(TestBit(test,i))
          {
            i = grid.ll[i];
            continue;
          }      

          // DISTANCIAS DE LA VECINA
          d2   = P[i].Dis;                    // distancia cosmologica
          dl2  = d2*(1.0+gal[i].red);         // distancia luminosidad 

          rm12 = rmaplim-25.0-5.0*log10((dl1+dl2)/2.);
          rmmin12 = rmapmin-25.0-5.0*log10((dl1+dl2)/2.);

          // CALCULA LA INTEGRAL DEL FACTOR DE ESCALA
          itabla=(int)((rm12-mt2min)/dmt2);
          imintabla=(int)((rmmin12-mt1min)/dmt1);

          if(itabla<0 || itabla >= NTABLA )
          {
            printf("XXXX itabla= %d\n ",itabla);
            exit(-33);
          }

          if(imintabla<0 || imintabla >= NTABLA )
          {
            printf("XXXX imintabla= %d\n ",imintabla);
            exit(-33);
          }

          rint12  = (rinttabla[imintabla][itabla+1]-rinttabla[imintabla][itabla])/dmt2;
          rint121 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla][itabla];

          rint12  = (rinttabla[imintabla+1][itabla+1]-rinttabla[imintabla+1][itabla])/dmt2;
          rint122 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla+1][itabla];

          rint12 = (rint122-rint121)/dmt1;
          rint12 = rint12*(rmmin12-imintabla*dmt1-mt1min)+rint121;

          rscale = pow(rint12/rintlim,-0.33333);
          dl = iden.d0*rscale;
          vl = v0*rscale;

          alfa12 = gal[i].alfa-gal[ic].alfa;
          tita12 = sin(gal[i].delta)*sin(gal[ic].delta)+cos(gal[i].delta)*cos(gal[ic].delta)*cos(alfa12);

          numera = (cos(gal[i].delta)*sin(alfa12));
          numera *= numera;

          tridif = cos(gal[ic].delta)*sin(gal[i].delta)-sin(gal[ic].delta)*cos(gal[i].delta)*cos(alfa12);
          tridif *= tridif;

          numera += tridif;
          numera = sqrt(numera);

          tita12 = atan2(numera,tita12);

          tita12=acos(tita12);
          d12 = sin(tita12/2.0)*(d1+d2);
          v12 = fabs(d1-d2);

          if((d12<dl) && (v12<vl))
          {
             #ifdef LOCK
               Unir(ic,i,lock);
             #else
               Unir(ic,i);
             #endif

             if(i>ic)
             {
              SetBit(test,i);
              #ifdef LOCK
                busv_rec(i,test,lock);
              #else
                busv_rec(i,test);
              #endif
             }
          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

}

void linkedlist(int *test)
{
  int i,g;

  Temp.ll = (int *) calloc(cp.npart,sizeof(int));

  iden.ngrupos = 0;
  for(i=0;i<cp.npart;i++)
  {
    P[i].gr = Raiz(i);
    if(TestBit(test,P[i].gr))
    {
      Temp.ll[P[i].gr]++;
      if(Temp.ll[P[i].gr]>=NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[P[i].gr] = iden.ngrupos;
        ClearBit(test,P[i].gr);
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (int *) malloc(iden.ngrupos*sizeof(int));
  Temp.npgrup = (unsigned int *) malloc(iden.ngrupos*sizeof(unsigned int));

  for(i=0;i<cp.npart;i++)
  { 
    if(TestBit(test,P[i].gr))
      P[i].gr = 0;
    else
      P[i].gr = Temp.ll[P[i].gr];

    if(i<iden.ngrupos)
    {
      Temp.head[i] = -1;
      Temp.npgrup[i] = 0;
    }
  }

  for(i=0;i<cp.npart;i++)
  {
    g = P[i].gr;

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }
}
