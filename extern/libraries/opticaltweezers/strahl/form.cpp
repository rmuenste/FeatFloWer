/***************************************************************************
                          form.cpp  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "form.h"
#include "misc.h"
#include "matrix.h" 
Form::Form(){
  r0=1;
  H=unity();
  R=unity();
  sf=1.0;
  Ealpha=0;
  Ebeta=0;
  Egamma=0; 
}

Form::Form (const Vector<double> &P,
        complex<double>  n,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex,
        const Vector<double> &Ey,
        const Vector<double> &Ez,
        const int type
       )
{ 
 Ealpha=0.0;
 Ebeta=0.0;
 Egamma=0.0;

 inelactive=true;
 sf=1.0;
 r0=1.0;
 this->type=type;
 this->P=P;
 this->alpha=alpha;
 this->n=n;
 ninel=n;
 e[0]=Ex;
 e[1]=Ey;
 e[2]=Ez;
 R(0,0)=Ex[0]; R(0,1)=Ex[1]; R(0,2)=Ex[2];
 R(1,0)=Ey[0]; R(1,1)=Ey[1]; R(1,2)=Ey[2];
 R(2,0)=Ez[0]; R(2,1)=Ez[1]; R(2,2)=Ez[2];

 H(0,0)=Ex[0]; H(0,1)=Ey[0]; H(0,2)=Ez[0];
 H(1,0)=Ex[1]; H(1,1)=Ey[1]; H(1,2)=Ez[1];
 H(2,0)=Ex[2]; H(2,1)=Ey[2]; H(2,2)=Ez[2];
}

Form::Form (const Form &F)
{
 inelactive=true; 
 sf=1.0; 
 type=F.type;
 P=F.P;
 alpha=F.alpha;
 Ealpha=F.Ealpha;
 Ebeta=F.Ebeta;
 Egamma=F.Egamma;
 H=F.H;
 R=F.R;
 por=F.por;
 pul=F.pul;
 r0=F.r0;
 type=F.type;
 n=F.n;
 ninel=F.ninel;
 for (int i=0; i<3; i++)
 e[i]=F.e[i];
 setMatrix (Ealpha,Ebeta,Egamma);
}


/*Form & Form::operator = (Form &f)
{
 if (this == &f) return *this;
 H=f.H;
 P=f.P;
 R=f.R;
 alpha=f.alpha;
 e=f.e;
 n=f.n;
 por=f.por;
 pul=f.pul;
 //r=f.r;
 r0=f.r0;
 type=f.type;
 return *this;
} */

void Form::setMatrix (Matrix<double> H)
{
	this->H=H;
	R=invert(H);
	// R=invert(H);
 e[0]=H*ex;
 e[1]=H*ey;
 e[2]=H*ez;
}

void Form::setMatrix (double alpha, double beta, double gamma)
{  
 double ca,cb,cg;
 double sa,sb,sg;
 ca=cos(alpha); cb=cos(beta); cg=cos(gamma);
 sa=sin(alpha); sb=sin(beta); sg=sin(gamma);

 H(0,0)= cg*cb;   H(0,1)= sa*sb*cg+ca*sg; H(0,2)=-ca*sb*cg+sa*sg;
 H(1,0)=-cb*sg;   H(1,1)=-sa*sb*sg+ca*cg; H(1,2)= ca*sb*sg+cg*sa;
 H(2,0)=    sb;   H(2,1)=-sa*cb;          H(2,2)= ca*cb;

 R(0,0)=cb*cg;             R(0,1)=-cb*sg;          R(0,2)= sb;
 R(1,0)=ca*sg+sa*sb*cg;    R(1,1)=ca*cg-sa*sb*sg;  R(1,2)=-sa*cb;
 R(2,0)=sa*sg-ca*sb*cg;    R(2,1)=sa*cg+ca*sb*sg;  R(2,2)=ca*cb;
 
 // R=invert(H);
 e[0]=H*ex;
 e[1]=H*ey;
 e[2]=H*ez;
 Ealpha=alpha;
 Ebeta=beta;
 Egamma=gamma;
}

/*
void Form::initQuad()
{
 initInc(this);
}*/


/*void Form::setr0(double r0)
{
 setR0(this,r0);
}*/

Matrix<double> computeInertia(Form *F)
{
    switch (F->type)
	{
        case ELLIPSOID : return ((FormEllipsoid *)F)->computeInertia();
		case SURFACE   : return ((surface *)F)->computeInertia();
	}
}