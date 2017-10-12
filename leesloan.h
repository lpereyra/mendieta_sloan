#ifndef LEESLOAN_H
#define LEESLOAN_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define PI180 0.01745329251994329577
#define CVEL  299792.458
#define MAGMENOSINF -26.0

struct SnapST
{
  int nfiles;
  char root[200], name[200]; 
} snap;

struct paramfl
{
  double ma   ;
  double alfa ;
  double fia  ;
};

struct paramcos
{
  double omegam   ;
  double omegal   ;
  double omegak   ;
};

struct galsloan
{
  long int id,targetid,specid;
  float alfa,delta,red,rederr;
  float m[5],merr[5],modelm[5],modelmerr[5],ext[5];
  float petrosg[3],petrosgerr[3],petrosr[3],petrosrerr[3];
  float k[5];
} *gal;

struct galsloang
{
  struct galsloan gal;
  int grupo[2], npgrupo[2];
} *gr;

void readsloan();
double funlum(double x, void *p);
double intfl(double x1, double x2);
double f(double z, void *p);
double red2dis(double z);
void change_positions(int n);

#endif
