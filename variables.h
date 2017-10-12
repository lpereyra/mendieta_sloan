#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NPARTMIN
  #define NPARTMIN 4 
#endif

#ifndef NTABLA
  #define NTABLA 1000
#endif

/* Posiciones, velocidades y energias de las partículas */
struct particle_data 
{
  float          Pos[3];
  float          Dis;
  int            sub;
  int            gr;
} *P;

int     nfrac               ;
double  *fof                ;
float   pmin[3]             ;
float   pmax[3]             ;
double  rmaplim             ;  // MAGNITUD APARENTE LIMITE DEL CATALOGO    
double  rmapmin             ;  // MAGNITUD APARENTE MINIMA DEL CATALOGO    
double  redmax              ;  // REDSHIT MAXIMO                           
double  redmin              ;  // REDSHIT MINIMO                           
double  vfid                ;  // VELOCIDAD FIDUCIAL                       
double  flfia               ;  // AMPLITUD DE LA FL                        
double  flma                ;  // MAGNITUD CARACTERISTICA DE LA FL         
double  flalfa              ;  // PENDIENTE EN EL EXTREMO DEBIL DE LA FL   
double  v0                  ;  // EN MPC H^-1                              
double  rmablim             ;  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
double  rmabmin             ;  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
double  rintlim             ;  // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
#ifdef LIM_VOL
double  zcut;
#endif

// AUXILIARES TABLA //
double mt1min                        ;  
double mt1max                        ;
double mt2min                        ;  
double mt2max                        ;
double dmt1                          ;
double dmt2                          ;
double **rinttabla                   ;

void init_variables(int argc, char **argv);

#endif
