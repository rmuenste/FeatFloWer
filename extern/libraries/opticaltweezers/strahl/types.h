#include <string>
#include <istream>

#ifndef TYPES_H
#define TYPES_H 

#define XFIG_COLOR_PSEUDO_OBJECT 0
#define XFIG_ELLIPSE             1
#define XFIG_POLYLINE            2
#define XFIG_SPLINE              3
#define XFIG_TEXT                4
#define XFIG_ARC                 5
#define XFIG_COMPOUND            6

#define XFIG_COLOR_DEFAULT      -1
#define XFIG_COLOR_BLACK         0
#define XFIG_COLOR_BLUE          1
#define XFIG_COLOR_GREEN         2
#define XFIG_COLOR_CYAN          3
#define XFIG_COLOR_RED           4
#define XFIG_COLOR_MAGENTA       5
#define XFIG_COLOR_YELLOW        6
#define XFIG_COLOR_WHITE         7

using namespace std;

typedef struct 
{
 string orientation;
 string justification;
 string units;
 string papersize;
 
 float magnification; 
 string multiple_page;
 int transparent_color;
 int resolution;
 int coord_system;
} xfig_header;

typedef struct
{
 int object_code;
 int color_number;
 string rgb_values;
} xfig_color_pseudo_object;

typedef struct
{
 int arrow_type;           //   (enumeration type)
 int arrow_style;          //   (enumeration type)
 float arrow_thickness;    //     (1/80 inch)
 float arrow_width;        //     (Fig units)
 float arrow_height;       //     (Fig units)
} xfig_arrow;

class xfig_polyline
{
 void addPoint (int Px, int Py);
 public :
 xfig_polyline();
 int object_code;         //  (always 2)
 int sub_type;            /* (1: polyline
                              2: box
                              3: polygon
                              4: arc-box)
                              5: imported-picture bounding-box)*/
 int line_style;         //   (enumeration type)
 int thickness;          //   (1/80 inch)
 int pen_color;          //   (enumeration type, pen color)
 int fill_color;         //   (enumeration type, fill color)
 int depth;              //   (enumeration type)
 int pen_style;          //   (pen style, not used)
 int area_fill;          //   (enumeration type, -1 = no fill)
 float style_val;        //   (1/80 inch)
 int join_style;         //   (enumeration type)
 int cap_style;          //   (enumeration type, only used for POLYLINE)
 int radius;             //   (1/80 inch, radius of arc-boxes)
 int forward_arrow;      //   (0: off, 1: on)
 int backward_arrow;     //   (0: off, 1: on)
 int npoints;            //   (number of points in line)
 xfig_arrow fwa_def, bwa_def; 
 int *x;
 int *y;
};

typedef struct
{
 int object_code;
 int upperleft_corner_x;
 int upperleft_corner_y;
 int lowerright_corner_x;
 int lowerright_corner_y;
 int npolylines,ncpo;
 xfig_polyline *p;
} xfig_compound;

typedef struct
{ 
  xfig_header header;
  int ncpo, npoly, ncompounds;
  xfig_color_pseudo_object *cpo;
  xfig_polyline *polyline;
  xfig_compound *comp;
} xfig_file;


bool spyComment (istream &is, const bool jump=true);
int spyObjectType (istream &is);
istream& operator >> (istream& is, xfig_color_pseudo_object &cpo );
istream& operator >> (istream& is, xfig_polyline &p);
istream& operator >> (istream& is, xfig_header &h);
istream& operator >> (istream& is, xfig_compound &h);
istream& operator >> (istream& is, xfig_arrow &a);


ostream& operator << (ostream& os, xfig_color_pseudo_object cpo );
ostream& operator << (ostream& os, xfig_polyline p);
ostream& operator << (ostream& os, xfig_header h);
ostream& operator << (ostream& os, xfig_compound h);
ostream& operator << (ostream& os, xfig_arrow a);

#endif 
