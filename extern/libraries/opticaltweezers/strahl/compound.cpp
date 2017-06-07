/***************************************************************************
 *   Copyright (C) 2005 by Thomas Weigel                                   *
 *   weigel@lat.ruhr-uni-bochum.de                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "compound.h"
#include "misc.h"
extern long int deletedItems, createdItems;

Compound::Compound( const Form& f) : Form (f)
{
 stack=TStack(AnzEl);
 type=COMPOUND;
 lastHit=-1;
 initQuad();
 }

Compound::Compound(const Compound& C) : Form(C)
{
 AnzEl=C.AnzEl;
 copyFormList( Element,C.Element,AnzEl);
 stack=C.stack;
 lastHit=C.lastHit;
 initQuad();
}

Compound& Compound::operator = (Compound &f)
{
 if (this == &f) return *this;
 H=f.H;
 P=f.P;
 R=f.R;
 alpha=f.alpha;
 for (int i=0; i<3; i++) 
 e[i]=f.e[i];
 n=f.n;
 por=f.por;
 pul=f.pul;
 for (int l=0; l<AnzEl; l++)
   deleteInc( Element[l]);
 free(Element);
 copyList(f.Element,f.AnzEl);
 type=f.type;
 Ealpha=f.Ealpha;
 Ebeta=f.Ebeta;
 Egamma=f.Egamma;
 stack=f.stack;  
 return *this;
}

Compound& Compound::operator = (Compound f)
{
 if (this == &f) return *this;
 H=f.H;
 P=f.P;
 R=f.R;
 alpha=f.alpha;
 for (int l=0; l<AnzEl; l++)
   deleteInc( Element[l]);
 free(Element);
 for (int i=0; i<3; i++) 
 e[i]=f.e[i];
 n=f.n;
 por=f.por;
 pul=f.pul;
 copyList(f.Element,f.AnzEl);
 type=f.type;
 Ealpha=f.Ealpha;
 Ebeta=f.Ebeta;
 Egamma=f.Egamma;
 stack=f.stack;
 return *this;
}

Compound::Compound( const Vector<double> &P,
             Form **F,
             int AnzEl,
             complex<double> n,
             double r0,
             const Matrix<complex<double> > alpha,
             const Vector<double> &Ex,
             const Vector<double> &Ey, 
             const Vector<double> &Ez
             ) : Form (P, n, alpha, Ex, Ey, Ez, COMPOUND)
{
 copyFormList( Element,F,AnzEl);
 Ealpha=0;
 Ebeta=0;
 Egamma=0;
  trafo (Ex,Ey,Ez,H,R);
 stack=TStack(AnzEl);
 this->AnzEl=AnzEl;
  this->r0=r0;
 this->P=P;
 type=COMPOUND;
 lastHit=-1; 
 initQuad();
}


Compound::Compound()
 : Form()
{
 AnzEl=0;
 Element=0;
 type=COMPOUND;
 lastHit=-1;
 r0=1;
}


Compound::~Compound()
{
 if (AnzEl>0)
 { 
  for (int l=0;l<AnzEl; l++) 
   deleteInc (Element[l]);
  free(Element);
  AnzEl=0; 
 }
 stack.clear();
}




/*!
    \fn Compound::addElement(Form *El)
 */

void Compound::addList(Form **El, int n)
{
 for (int i=0; i<n; i++)
  addElement( El[i]);
 initQuad();
}

void Compound::copyList(Form **El, int n)
{
 if (AnzEl>0) { for (int l=0; l<AnzEl; l++) deleteInc(Element[l]); free(Element); }
 copyFormList( Element,El,n);
 AnzEl=n;

 /*for (int i=0; i<n; i++)
  addElement( El[i]);*/
}

 
void Compound::addElement(Form *El)
{
  if (AnzEl==0) Element=(Form **) malloc (sizeof (Form *));
  else Element=(Form **) realloc (Element,sizeof (Form *) * (AnzEl+1));
  Element [AnzEl]=El;
  AnzEl++;
  initStack(AnzEl);   
  initQuad();
}



/*!
    \fn Compound::deleteElement(int Item)
 */
void Compound::deleteElement(int Item)
{
    if ((Item<=AnzEl) && (Item>0) && (AnzEl>0))
   {
    if (Item < AnzEl-1) AnzEl--;
    {
     for (int i=Item+1; i<AnzEl; i++)
       Element[i-1]=Element[i];
    }
    AnzEl--;
   }
initStack(AnzEl); 
initQuad();
}



/*!
    \fn Compound::scale (double sf)
 */
void Compound::scale (double sf)
{
    for (int i=0; i<AnzEl; i++)
      Element[i]->scale(sf);
	initQuad();
}



/*!
    \fn Compound::isInside (const Vector<double> &p)
 */
bool Compound::isInside (const Vector<double> &p)
{
	bool found=false;
    for (int i=0; i<AnzEl; i++)
    found=found || Element[i]->isInside(p);
    return found;
}



/*!
    \fn Compound::Volume()
 */
double Compound::Volume()
{
    double Erg=0;
    for (int i=0; i<AnzEl; i++)
      Erg+=Element[i]->Volume();
    return Erg;
}



bool Compound::next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
 bool found;
 int Item;
 Vector<double> k,p=Ps-P;
 p=H*p;
 k=H*K;
 if (!inside) // Strahl kommt von auﬂen
 {
  if (find_next(p,k,pout,Item)) 
  {  
   pout=R*pout+P;
   stack.clear();
   lastHit=Item;
   stack.newItem(Item);  
//    cout << "Item=" << Item << "  p=" << p << endl;       
   return true;
  }
  else 
  {
   pout=Ps;
   return false;
  }
 }
 else // Strahl ist Innen
 {
  if (stackSize==0) stack.newItem(lastHit);
  do 
  {
   found=find_next(p,k,pout,Item); // erst mal den Schnittpunkt mit dem n‰chsten Element suchen
   if (found) // Schnittpunkt gefunden -> im Stack nachschauen ggf. einf¸gen
    {
     p=pout;
     stack.checkItem (Item);     
    }  
  } while((stack.numItems()>0) && found); 
  lastHit=Item;
   pout=R*pout+P;
  return true;
 }
}



/*!
    \fn Compound::next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside=-1)
 */
/*bool Compound::next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
  struct TStack *hs;
  bool found;
  int Item,Item_new;  
  Vector<double> p_new;
  Item=-1;
   Vector<double> k,p=Ps-P;
 p=H*p;
 k=H*K;
 pout=Vector<double> (-1000,-1000,-1000);
 if (!inside)  // kommt von ausserhalb
 { 
  if (find_next(p,k,pout,Item)) // Schnittpunkt mit einem der Elemente gefunden
   {
    lastHit=Item;
    if (stack!=0) { cout << "stack gibts ja schon !! next:" << stack->next << endl; delete stack; }
    stack=new struct TStack;    // Stack auf das gefundene Element setzen
    stack->prev=0;
    stack->next=0; 
    stack->id=Item;
    stackSize++; 
    return true;
   } // if (find_next...)
   else return false;
 } // if (!inside)
 else // im Einschluss 
 {
   if (stackSize==0)  
   {
	stack=new struct TStack;
    stack->prev=0;
    stack->next=0; 
    stack->id=lastHit;
    stackSize++; 
   }
   pout=p; 
   do 
   {     
     found=find_next(pout,k,p_new,Item_new);   // n‰chsten Schnittpunkt mit einem Element innerhalb des Compounds finden     
     lastHit=Item_new;
     if (found && (stackSize>0))  // Es wurde was gefunden 
     {
      hs=stack;    
     while ((hs->next!=0) && (hs->id!=Item_new) )  // Erst mal schauen, ob der Einschluss schon mal getroffen wurde
     {
      hs=hs->next;
     }
// Der Einschluss wurde schon einmal getroffen 
// -> wir sind jetzt ausserhalb dieses Einschlusses 
// -> Eintrag aus Stack lˆ    /// @todo implement meschen	

     if ((hs->id==Item_new) && (stack!=0)) 
     {
      if ((hs->next==0) && (hs->prev==0)) {  cout << "hs=" << hs << "   stack=" << stack << endl; lastHit=hs->id; delete hs; pout=p_new; stackSize=0; stack=0; return true;}      // letzter Einschluss im Stack soll gelˆscht werden -> fertig !
	  if (hs->prev==0) // erstes Element in der Liste wird gelˆscht
	  {
	   hs->next->prev=0; 
	   stack=hs->next;
	   delete hs;
	  } // if (hs->prev==0)
	  else 
	  {	   
	   if (hs->next==0) // letztes Element in der Liste wird gelˆscht
	   {
	    hs->prev->next=0;
	    delete hs; 
	   }
	   else // ein Element in der Mitte wird gelˆscht
	   {
 		hs->prev->next=hs->next;
 		hs->next->prev=hs->prev;
		delete hs;
       }
	   stackSize--;
      }
	 } 
     else // Ein neues Element wurde gefunden ==> wird vorne zugef¸gt
     {
      hs=new struct TStack;
	  hs->prev=0;
	  hs->next=stack;
	  hs->id=Item_new;          
	  stack->prev=hs;
	  stack=hs;
	  stackSize++;
     }
      pout=p_new;
      Item=Item_new; 
     }
   } while (found && stack!=0);    
   lastHit=Item;
   pout=p_new;
   if (stack!=0) cout << "Das ist aber komisch ?!  stack=" << stack << "   stackSize=" << stackSize << endl;
   stack=0;
   stackSize=0;
   return true;
  }  // if 
  return true; 
}*/

Vector<double> Compound::norm (const Vector<double> &r)
{
 Vector<double> p=r-P;
 p=H*p; 
 if ((lastHit<=-1) || (lastHit>AnzEl)) return zero;
 return R*Element[lastHit]->norm(p); 
}

/*!
    \fn Compound::find_next(const Vector<double> &p, const Vector<double> &k, int &Item)
 */
bool Compound::find_next(const Vector<double> &p, const Vector<double> &k, Vector<double> &P, int &Item)
/* 
 Finde einen Schnittpunkt des Strahls mit dem n‰chsten Element (egal
ob der Strahl von Innen oder von Aussen kommt) 
*/
{
 double d,dh;
 Vector<double> Ph;
 bool found,foundh;
 d=-1;
 P=p; 
 found=false;
 for (int i=0; i<AnzEl; i++) 
 {
  foundh=Element[i]->next(p,k,Ph);
  dh=abs(Ph-p);
  if ( ((d<0) || (d>dh)) && foundh ) 
  {
   Item=i;
   d=dh;   
   P=Ph;
   found=true;
  }
 }
 return found; 
}

 void Compound::binWrite (ofstream &os)
 {
  int i;
  P.binWrite(os);
  H.binWrite(os);
  R.binWrite(os);
  os.write ((char *) &n, (char) sizeof (n)); 
  alpha.binWrite(os);
  pul.binWrite (os);
  por.binWrite (os);
  for (i=0; i<3; i++)
    e[i].binWrite(os);
  os.write ((char *) &Ealpha,(char) sizeof (Ealpha));
  os.write ((char *) &Ebeta,(char) sizeof (Ebeta));
  os.write ((char *) &Egamma,(char) sizeof (Egamma));
  os.write ((char *) &sf,(char) sizeof (sf));
  os.write ((char *) &r0,(char) sizeof (r0));
  os.write ((char *) &AnzEl,(char) sizeof (AnzEl));
  for (i=0; i<AnzEl; i++)
  {
   binWriteInc( os,Element[i]);   
  }
 }

 void Compound::binRead(ifstream &os)
{
  int i;
  P.binRead(os);
  H.binRead(os);
  R.binRead(os);
  os.read ((char *) &n, (char) sizeof (n)); 
  alpha.binRead(os);
  pul.binRead (os);
  por.binRead (os);
  for (i=0; i<3; i++)
    e[i].binRead(os);
  os.read ((char *) &Ealpha,(char) sizeof (Ealpha));
  os.read ((char *) &Ebeta,(char) sizeof (Ebeta));
  os.read ((char *) &Egamma,(char) sizeof (Egamma));
  os.read ((char *) &sf,(char) sizeof (sf));
  os.read ((char *) &r0,(char) sizeof (r0));
  os.read ((char *) &AnzEl,(char) sizeof (AnzEl));
  Element=(Form **) malloc (sizeof(Form *) * AnzEl);
  for (i=0; i<AnzEl; i++)
  {   
   binReadInc(os,Element[i],true);   
  }
 }

void Compound::setr0(double r0)
{
  P=P/this->r0*r0;
  for (int i=0; i<AnzEl; i++)
    Element[i]->setr0(r0);
 this->r0=r0;
  initQuad();
}

void Compound::initQuad()
{
  double xmax,xmin,ymin,ymax,zmax,zmin;	 
  xmin=r0;
  ymin=r0;
  zmin=r0;
  xmax=-r0; 
  ymax=-r0;
  zmax=-r0;
  cout << "COMPOUND::INIT" << endl;
	
  for (int i=0; i<AnzEl; i++)
  {
   Element[i]->initQuad();
   if (Element[i]->pul[0]<xmin) xmin=Element[i]->pul[0];
   if (Element[i]->pul[1]<ymin) ymin=Element[i]->pul[1];
   if (Element[i]->pul[2]<zmin) zmin=Element[i]->pul[2];

   if (Element[i]->por[0]>xmax) xmax=Element[i]->por[0];
   if (Element[i]->por[1]>ymax) ymax=Element[i]->por[1];
   if (Element[i]->por[2]>zmax) zmax=Element[i]->por[2];
  }
  pul=Vector<double>(xmin,ymin,zmin);
  por=Vector<double>(xmax,ymax,zmax);
  cout << "pul=" << pul << "   por=" << por << endl;
}

/*!
    \fn Compound::clearStack()
 */

void TStack::newItem (int Item)
{
 /*if (Anz==0) 
  ID=(int *) malloc (sizeof(int));
 else 
  ID=(int *) realloc (ID,sizeof (int) * (Anz+1));*/
 if (Anz<0) Anz=0;
 if (Anz<MAXANZSTACKITEMS)
 ID[Anz]=Item;
 Anz++;
}

void TStack::deleteItem (int i)
{
 if ((i<Anz-1) &&(i>=0) ) 
  for (int l=i; l<Anz-2; l++)
     ID[l]=ID[l+1];
  Anz--;
}

void TStack::clear()
{
 
 if (maxAnz>0) { /*delete[] ID;*/  }
 Anz=0;
 maxAnz=0;
}

bool TStack::checkItem(int Item)
{
 bool found=false;
 for (int l=0; (l<Anz) && !found; l++)
  if (ID[l]==Item) { deleteItem(l); found=true; }
 if (!found) newItem( Item);
 return found;
}

TStack& TStack::operator = (TStack stack)
{ 
 if (this==&stack) return *this;
 // clear();
 Anz=stack.numItems();
 maxAnz=stack.maxAnz; 
// ID=(int *) malloc (sizeof(int) * Anz);
 
 for (int l=0; l<Anz; l++)
  ID[l]=stack.Id(l);
 return *this;
}

/*!
    \fn Compound::setP(Vector<double> r)
 */
void Compound::setP(Vector<double> r)
{
    Vector<double> dP=r-P;
    for (int i=0; i<AnzEl; i++)
    Element[i]->setP(Element[i]->P+dP);
    P=r;
}


/*!
    \fn Compound::setP(double x, double y, double z)
 */
void Compound::setP(double x, double y, double z)
{
  Vector<double> r(x,y,z);
  setP(r);
}
