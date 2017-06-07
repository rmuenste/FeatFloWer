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
#ifndef COMPOUND_H
#define COMPOUND_H

#include "form.h"

/*struct TStack
{
 int id;
 struct TStack *next;
 struct TStack *prev;
};*/

#define MAXANZSTACKITEMS 200
class TStack 
{
 public :
  TStack () {maxAnz=0; Anz=0; }
  TStack (int Anz) { /*ID=new int[Anz]; */ maxAnz=Anz; this->Anz=0;}
  TStack& operator = (TStack stack);
  void reset () { Anz=0; } 
  void newItem(int Item);
  void deleteItem(int i);
  bool checkItem(int Item);
  int numItems() {return Anz;}
  int Id(int Item) {return ID[Item]; }
  void clear();
  ~TStack() {clear();}
 //protected :
  int Anz,maxAnz;
  int ID[MAXANZSTACKITEMS]; 
};
 

/**
Zusammenschluss verschiedener Form-Objekte zu einem einzigen Objekt.


	@author Thomas Weigel <weigel@lat.ruhr-uni-bochum.de>
*/
class Compound : public Form
{
public:
    Compound(); 
    Compound(const Form&);
    Compound(const Compound&); 
	Compound& operator = (Compound &f);
	Compound& operator = (Compound f);
    Compound(
             const Vector<double> &P,
             Form **F,
             int AnzEl,
             complex<double> n,
             double r0=1.0,
             const Matrix<complex<double> > alpha=CUNITY,
             const Vector<double> &Ex=ex,
             const Vector<double> &Ey=ey,
             const Vector<double> &Ez=ez
             );
    ~Compound();
    Form **getElement () {return Element; }
	Form *getElement (int i) {return Element[i]; }
    int getAnzElements() {return AnzEl;}
    void addList(Form **El, int n);
    void copyList(Form **El, int n); 
    void addElement(Form *El);
    void deleteElement(int Item);
    void scale (double sf);
    bool isInside (const Vector<double> &p);
    double Volume();
    bool next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside=-1);
    bool find_next(const Vector<double> &p, const Vector<double> &k, Vector<double> &P, int &Item);
   Vector<double> norm (const Vector<double> &P);
    void binWrite (ofstream &os);
    void binRead(ifstream &os);
    void setr0(double r0);   
    void initStack(int AnzEl) {stack.clear(); stack=TStack(AnzEl);}    
	void initQuad();
    void setP(Vector<double> r);
    void setP(double x, double y, double z);
//private:
  Form **Element;
  int AnzEl;
  int lastHit;
  TStack stack;
  int stackSize;
protected:    
};

#endif
