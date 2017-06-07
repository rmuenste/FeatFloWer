#include "types.h"
#include <iostream>
#include <stdlib.h>

xfig_polyline::xfig_polyline()
{
 object_code=XFIG_POLYLINE;
 pen_style=0;
 pen_color=XFIG_COLOR_DEFAULT;
 area_fill=-1;
 forward_arrow=0;
 backward_arrow=0;
 npoints=0;
 x=0;
 y=0;
}

void xfig_polyline::addPoint (int Px, int Py)
{
 if (npoints==0) 
 {
  x=(int *) malloc (sizeof(int));
  y=(int *) malloc (sizeof(int));
 }
 else 
 {
  x=(int *) realloc (x, sizeof(int) * (npoints+1));
  y=(int *) realloc (y, sizeof(int) * (npoints+1));
  npoints++;
 }
 x[npoints-1]=Px;
 y[npoints-1]=Py;
}


using namespace std;
istream& operator >> (istream& is, xfig_color_pseudo_object &cpo )
{
 is >> cpo.object_code >> cpo.color_number >> cpo.rgb_values;
 return is;
}

istream& operator >> (istream& is, xfig_polyline &p)
{
is >> p.object_code;
is >> p.sub_type;
is >> p.line_style;
is >> p.thickness;
is >> p.pen_color;
is >> p.fill_color;
is >> p.depth;
is >> p.pen_style;
is >> p.area_fill;
is >> p.style_val;
is >> p.join_style;
is >> p.cap_style;
is >> p.radius;
is >> p.forward_arrow;
is >> p.backward_arrow;
is >> p.npoints; 
 if (p.forward_arrow) is >> p.fwa_def;
if (p.backward_arrow) is >> p.bwa_def;

 // if (p.npoints<0) { delete[] p.x; delete[] p.y; }
 p.x=new int[p.npoints];
 p.y=new int[p.npoints];
 for (int i=0; i<p.npoints; i++)
  is >> p.x[i] >> p.y[i];
 return is; 
}

bool spyComment (istream &is, bool jump)
{
 char buffer[256];
 int pos=is.tellg();
 is >> buffer;
 is.seekg(pos);
 char c=buffer[0];
 if (c=='#') 
 {
  if (jump)  
  {
   is.getline(buffer,254);
  }
  return true;
 }
 return false;
}

istream& operator >> (istream& is, xfig_header &h)
{
 is >> h.orientation;
 is >> h.justification;
 is >> h.units;
 is >> h.papersize;
 is >> h.magnification;
 is >> h.multiple_page;
 is >> h.transparent_color;
 spyComment(is);
 is >> h.resolution;
 is >> h.coord_system;
 return is;
}



ostream& operator << (ostream& os, xfig_color_pseudo_object cpo )
{
 os << cpo.object_code << "  " << cpo.color_number << "   " << cpo.rgb_values; 
 return os;
}

istream& operator >> (istream& is, xfig_arrow &a)
{
 is >> a.arrow_type;
 is >> a.arrow_style;
 is >> a.arrow_thickness;
 is >> a.arrow_width;
 is >> a.arrow_height;
 return is;
}

ostream& operator << (ostream& os, xfig_arrow a)
{
 os << a.arrow_type << "  ";
 os << a.arrow_style << "  ";
 os << a.arrow_thickness << "  ";
 os << a.arrow_width << "  ";
 os << a.arrow_height; 
 return os;
}

ostream& operator << (ostream& os, xfig_polyline p)
{
 char tab=9;
os << p.object_code << "   ";
os << p.sub_type << "   ";
os << p.line_style << "   ";
os << p.thickness << "   ";
os << p.pen_color << "   ";
os << p.fill_color << "   ";
os << p.depth << "   ";
os << p.pen_style << "   ";
os << p.area_fill << "   ";
os << fixed << p.style_val << "   ";
os << p.join_style << "   ";
os << p.cap_style << "   ";
os << p.radius << "   ";
os << p.forward_arrow << "   ";
os << p.backward_arrow << "   ";
os << p.npoints << endl;
if (p.forward_arrow) os << tab << p.fwa_def << endl;
if (p.backward_arrow) os << tab << p.bwa_def << endl;
os << tab; 
 
 for (int i=0; i<p.npoints; i++)
   os << p.x[i] << "   " << p.y[i] << "   ";
   os << endl;
 return os;
}

ostream& operator << (ostream& os, xfig_header h)
{
 os << h.orientation << endl;
 os << h.justification << endl;
 os << h.units << endl;
 os << h.papersize << endl;
 os << h.magnification << endl; 
 os << h.multiple_page << endl;
 os << h.transparent_color << endl;
 os << h.resolution << "  ";
 os << h.coord_system;
 return os;
}

istream& operator >> (istream& is, xfig_compound &c)
{
 char buffer[255];
 int pos,counter;
 char dummy[255];
 is >> c.object_code >> c.upperleft_corner_x >> c.upperleft_corner_y >> c.lowerright_corner_x >> c.lowerright_corner_y;
is.getline(buffer,254);
 do 
 {
  spyComment(is);
//  is.getline (buffer,254);
//  cout << "buffer=" << buffer << endl;
  switch (spyObjectType(is))
  {
   case XFIG_POLYLINE : 
                       if (c.npolylines>0)
     c.p=(xfig_polyline *) realloc (c.p,sizeof(xfig_polyline)*(c.npolylines+1)); 
                       else
     c.p=(xfig_polyline *) malloc (sizeof (xfig_polyline));
			is >> c.p[c.npolylines];
			c.npolylines++; 
//			cout << c.p[c.npolylines]<< endl << "---------------" << endl;
		       break;
  }
  pos=is.tellg();
   is >> dummy;
  is.seekg(pos);
  }
 while ((atoi(dummy)!=-6) && (!is.eof())); 
 return is;
}

ostream& operator << (ostream& os, xfig_compound c)
{
 os << c.object_code;
 os << "  " << c.upperleft_corner_x << "  " << c.upperleft_corner_y << "   ";
 os << "  " << c.lowerright_corner_x << "  " << c.lowerright_corner_y << "   ";
 os << endl;
 for (int i=0; i<c.npolylines; i++)
  os << c.p[i];
 os << "-6" << endl; 
 return os; 
}

int spyObjectType (istream &is)
{
 int pos=is.tellg();
 char c[255];
 is >> c;
 is.seekg(pos);
 return atoi(c);
}
