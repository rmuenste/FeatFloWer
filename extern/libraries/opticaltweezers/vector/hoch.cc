# include <iostream> 

double hoch (double x, int y)
{
int i;
double z;
 z = 1.0 ; 
 if (y == 0) return 1.0;

 if(y > 0)
 { 
  for (i = 1; i<=y; i++)
       {
	   z = z * x;  
      }
 }

 if (y < 0)
 { 
  for (i = 1; i<= -y ;i++)
  {
   z = ( 1/ x) * z;
  }
  } 

return z;
}

/*
main (int argc, char **argv)
{
    double a, v;
    int b;
 cout <<" gib 2 Zahlen";
 cin >>a >>b;
v = hoch (a ,b);
 cout << v  <<endl ;
     }
*/









