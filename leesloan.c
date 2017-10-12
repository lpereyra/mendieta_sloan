#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesloan.h"

void readsloan()
{

  FILE *pf;
  char filename[200];
  int i,j;
  float disred, dlum; 
  double mt1, mt2;
  double dismin, dismax;
  #ifdef LIM_VOL
  float mag_cut,mtest;
  #endif

  cp.omegam  = 0.3                       ;  /* OMEGA MATERIA                              */
  cp.omegal  = 0.7                       ;  /* OMEGA LAMBDA                               */
  cp.omegak  = 1.0-cp.omegam-cp.omegal   ;  /* OMEGA CURVATURA                            */
  cp.h0      = 100.                      ;  /* ESTO DEJA TODO EN UNIDADES DE H^-1         */
  cp.dlummax = -1.E26                    ;
  mt1min     = 1.E26                     ;
  mt1max     = -1.E26                    ;
  mt2min     = 1.E26                     ;
  mt2max     = -1.E26                    ;
  dismin     = 1.E26                     ;
  dismax     = -1.E26                    ;

  RED("Read Sloan...\n");

  sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");

  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&cp.npart, sizeof(int),1,pf);

  /* ALOCATACION Y LECTURA */
  gal = (struct galsloan *) malloc(cp.npart*sizeof(struct galsloan));
  P   = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));

  #ifdef LIM_VOL
    mag_cut = rmaplim-25.0-5.0*log10(red2dis(zcut)*(1.0+zcut));
    fprintf(stdout,"%f MABS CUT a %f redshift cut\n",mag_cut,zcut);
    fflush(stdout);
  #endif

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  i=0;
  for(j=0;j<cp.npart;j++)
  {
    fread(&gal[i].id,sizeof(long int),1,pf);
    fread(&gal[i].targetid,sizeof(long int),1,pf);
    fread(&gal[i].specid,sizeof(long int),1,pf);
    fread(&gal[i].alfa,sizeof(float),1,pf);
    fread(&gal[i].delta,sizeof(float),1,pf);
    fread(&gal[i].red,sizeof(float),1,pf);
    fread(&gal[i].rederr,sizeof(float),1,pf);
    fread(&gal[i].m[0],sizeof(float),1,pf);
    fread(&gal[i].m[1],sizeof(float),1,pf);
    fread(&gal[i].m[2],sizeof(float),1,pf);
    fread(&gal[i].m[3],sizeof(float),1,pf);
    fread(&gal[i].m[4],sizeof(float),1,pf);
    fread(&gal[i].merr[0],sizeof(float),1,pf);
    fread(&gal[i].merr[1],sizeof(float),1,pf);
    fread(&gal[i].merr[2],sizeof(float),1,pf);
    fread(&gal[i].merr[3],sizeof(float),1,pf);
    fread(&gal[i].merr[4],sizeof(float),1,pf);
    fread(&gal[i].modelm[0],sizeof(float),1,pf);
    fread(&gal[i].modelm[1],sizeof(float),1,pf);
    fread(&gal[i].modelm[2],sizeof(float),1,pf);
    fread(&gal[i].modelm[3],sizeof(float),1,pf);
    fread(&gal[i].modelm[4],sizeof(float),1,pf);
    fread(&gal[i].modelmerr[0],sizeof(float),1,pf);
    fread(&gal[i].modelmerr[1],sizeof(float),1,pf);
    fread(&gal[i].modelmerr[2],sizeof(float),1,pf);
    fread(&gal[i].modelmerr[3],sizeof(float),1,pf);
    fread(&gal[i].modelmerr[4],sizeof(float),1,pf);
    fread(&gal[i].ext[0],sizeof(float),1,pf);
    fread(&gal[i].ext[1],sizeof(float),1,pf);
    fread(&gal[i].ext[2],sizeof(float),1,pf);
    fread(&gal[i].ext[3],sizeof(float),1,pf);
    fread(&gal[i].ext[4],sizeof(float),1,pf);
    fread(&gal[i].petrosg[0],sizeof(float),1,pf);
    fread(&gal[i].petrosg[1],sizeof(float),1,pf);
    fread(&gal[i].petrosg[2],sizeof(float),1,pf);
    fread(&gal[i].petrosgerr[0],sizeof(float),1,pf);
    fread(&gal[i].petrosgerr[1],sizeof(float),1,pf);
    fread(&gal[i].petrosgerr[2],sizeof(float),1,pf);
    fread(&gal[i].petrosr[0],sizeof(float),1,pf);
    fread(&gal[i].petrosr[1],sizeof(float),1,pf);
    fread(&gal[i].petrosr[2],sizeof(float),1,pf);
    fread(&gal[i].petrosrerr[0],sizeof(float),1,pf);
    fread(&gal[i].petrosrerr[1],sizeof(float),1,pf);
    fread(&gal[i].petrosrerr[2],sizeof(float),1,pf);
    fread(&gal[i].k[0],sizeof(float),1,pf);
    fread(&gal[i].k[1],sizeof(float),1,pf);
    fread(&gal[i].k[2],sizeof(float),1,pf);
    fread(&gal[i].k[3],sizeof(float),1,pf);
    fread(&gal[i].k[4],sizeof(float),1,pf);

    if(gal[i].red<=redmin || gal[i].red >=redmax || gal[i].m[2] > rmaplim || gal[i].m[2] < rmapmin) continue;
    
    #ifdef LIM_VOL
    if(gal[i].red>zcut) continue;
    #endif

    if(gal[i].alfa<100. || gal[i].alfa>300) continue; // corta los bigotes

    if((gal[i].alfa>250.) && (gal[i].alfa<270) && (gal[i].delta>50.) && (gal[i].delta<70.)) continue; // corta una esquina

    gal[i].alfa  *= PI180;
    gal[i].delta *= PI180;

    disred = red2dis(gal[i].red);   /*EN MPC*/
    dlum = disred*(1.0+gal[i].red); 

    #ifdef LIM_VOL
    mtest = (gal[i].m[2]-gal[i].k[2])-25.0-5.0*log10(dlum);
    if(mtest>mag_cut) continue;
    #endif

    if(dlum>cp.dlummax) cp.dlummax = dlum; /* DISTANCIA LUMINOSIDAD MAXIMA */

    mt2  = rmaplim-25.0-5.0*log10(dlum);
    mt1  = rmapmin-25.0-5.0*log10(dlum) ;

    if(mt2>mt2max) mt2max = mt2; /*magnitud absoluta limite maxima*/
    if(mt2<mt2min) mt2min = mt2; /*magnitud absoluta limite minima*/

    if(mt1>mt1max) mt1max = mt1; /*magnitud absoluta minima maxima*/
    if(mt1<mt1min) mt1min = mt1; /*magnitud absoluta minima minima*/

    P[i].Pos[0] =  disred*cos(gal[i].delta)*cos(gal[i].alfa) ;
    P[i].Pos[1] =  disred*cos(gal[i].delta)*sin(gal[i].alfa) ;
    P[i].Pos[2] =  disred*sin(gal[i].delta)                  ;
    P[i].Dis    =  disred;

    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];

    if(P[i].Dis > dismax) dismax = P[i].Dis;
    if(P[i].Dis < dismin) dismin = P[i].Dis;

    i++;
  }

  fclose(pf);

  if(i!=cp.npart)
  {
    cp.npart = i;
    // REALLOCATEA //
    gal = (struct galsloan *) realloc(gal,cp.npart*sizeof(struct galsloan));
    P   = (struct particle_data *) realloc(P,cp.npart*sizeof(struct particle_data));
  }

  cp.vol = (2.210)*(pow(dismax,3.)-pow(dismin,3.))/3.0;

  fprintf(stdout,"Volumen aprox %.8e\n",cp.vol);
  fprintf(stdout,"Num Total %d\n",cp.npart);
  RED("End Read Sloan\n");
  fflush(stdout);
  
  #ifdef LIM_VOL
  rmablim = rmaplim-25.0-5.0*log10(red2dis(zcut)*(1.0+zcut)); // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = MAGMENOSINF;                                      // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  #else
  rmablim = rmaplim-25.0-5.0*log10(vfid);  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = rmapmin-25.0-5.0*log10(vfid);  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  #endif

  rintlim = intfl(rmabmin,rmablim);        // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
  fprintf(stdout, "INTEGRAL LIMITES\n");
  fprintf(stdout, "%f MABS correspondiente a %f MAPA con una velocidad fiducial %f\n",rmabmin,rmapmin,vfid);
  fprintf(stdout, "%f MABS correspondiente a %f MAPA con una velocidad fiducial %f\n",rmablim,rmaplim,vfid);
  fflush(stdout);

  mt2min *= (1.0-1.0e-2*mt2min/fabs(mt2min)); // agranda y achica los limites 
  mt2max *= (1.0+1.0e-2*mt2max/fabs(mt2max));
  mt1min *= (1.0-1.0e-2*mt1min/fabs(mt1min));
  mt1max *= (1.0+1.0e-2*mt1max/fabs(mt1max));

  // TABLA DE LA INTEGRAL DEL FACTOR DE ESCALA //
  dmt1 = (mt1max-mt1min)/(int)NTABLA;
  dmt2 = (mt2max-mt2min)/(int)NTABLA;

  printf("mt1min mt1max %f %f\n", mt1min, mt1max);
  printf("mt2min mt2max %f %f\n", mt2min, mt2max);

  fprintf(stdout,"Crea la tabla\n");

  rinttabla = (double **) malloc(NTABLA*sizeof(double *));
  for(i=0;i<NTABLA;i++)
    rinttabla[i] = (double *) malloc(NTABLA*sizeof(double));

  #ifdef NTHREADS
    #pragma omp parallel for default(none) \
    num_threads(NTHREADS) private(i,j,mt1,mt2) \
    shared(mt1min,dmt1,mt2min,dmt2,rinttabla)   
  #endif
  for(i=0;i<NTABLA;i++ )
  {
    mt1 = mt1min + (double)i*dmt1;
    for(j=0;j<NTABLA;j++ )
    {
      mt2 = mt2min + (double)j*dmt2;
      rinttabla[i][j] = intfl(mt1,mt2);
    }
  }

  fprintf(stdout,"Termina la tabla\n");

}

double funlum(double x, void *p)
{
  struct paramfl *par=(struct paramfl *)p ;
  double ma=(par->ma)                     ;
  double alfa=(par->alfa)                 ;
  double t                                ;

  t=pow(10.0,0.4*(ma-x)) ;
  t=pow(t,1.0+alfa)*exp(-t) ;

  return(t) ;
}

double intfl(double x1, double x2)
{
  double resultado, error;
  size_t neval;
  struct paramfl pfl;
  gsl_function F; 

  pfl.ma=flma;
  pfl.alfa=flalfa;
  pfl.fia=flfia;

  F.function = &funlum;
  F.params = &pfl;

  if(x1<MAGMENOSINF)x1=MAGMENOSINF;
  gsl_integration_qng(&F,x1,x2,1.0e-7,1.0e-7,&resultado,&error,&neval);
  return(resultado);
}

double f(double z, void *p) /* FUNCION A INTEGRAR PARA LA DISTANCIA EN FUNCION DE Z*/
{
  struct paramcos *par=(struct paramcos *)p ;
  double om=(par->omegam);
  double ol=(par->omegal);
  double ok=(par->omegak);
  double q;
  q=pow(1.+z,3.)*om+pow(1.+z,2.)*ok+ol;
  return(1.0/sqrt(q));
}

double red2dis(double z)
{
  double resultado, error;
  size_t neval;
  struct paramcos pcos;
  gsl_function F;

  pcos.omegam=cp.omegam;
  pcos.omegal=cp.omegal;
  pcos.omegak=cp.omegak;

  F.function=&f;
  F.params=&pcos;

  gsl_integration_qng(&F,0.0,z,1.0e-7,1.0e-7,&resultado,&error,&neval);
 
  return(CVEL/cp.h0*resultado);
}

void change_positions(int n)
{
  int ip, idim;
  #ifdef PRINT_XYZ
  FILE *pfout;
  char filename[200];
  #endif
  RED("Inicio Change Positions\n");

  #ifdef PRINT_XYZ

  #ifdef LIM_VOL
    sprintf(filename,"sloan_%.2f_xyz.bin",zcut);
  #else
    sprintf(filename,"sloan_xyz.bin");
  #endif
  pfout=fopen(filename,"w");
  fwrite(&n,sizeof(int),1,pfout);

  #endif

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip = 0; ip < n; ip++)
  {
    #ifdef PRINT_XYZ
    fwrite(&P[ip].Pos[0],sizeof(float),1,pfout);
    fwrite(&P[ip].Pos[1],sizeof(float),1,pfout);
    fwrite(&P[ip].Pos[2],sizeof(float),1,pfout);
    #endif
    P[ip].sub = 0;
    for(idim = 0; idim < 3; idim++)
      P[ip].Pos[idim] -= pmin[idim];
  }

  #ifdef PRINT_XYZ
  fclose(pfout);
  #endif

  cp.lbox = pmax[0] - pmin[0];
  for(idim = 1; idim < 3; idim++)
    if(cp.lbox < (pmax[idim] - pmin[idim])) cp.lbox = (pmax[idim] - pmin[idim]);

  cp.lbox *= 1.001;
  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
}

