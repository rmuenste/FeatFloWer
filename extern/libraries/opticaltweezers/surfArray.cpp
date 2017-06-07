#include "surfArray.h"


surfArray::surfArray()
{
	nObjs = 0;
	N = 0;
	type = 0;
}

surfArray::surfArray(const surfArray &SA)
{
	nObjs = SA.nObjs;
	data = new Vector<double>*[SA.nObjs];
	
	P = new Vector<double>*[nObjs];
	N = new int[nObjs];
	for (int i = 0; i < nObjs; i++)
	{
		N[i] = SA.N[i];
		data[i] = new Vector<double>[N[i]];
 		P[i] = new Vector<double>[N[i]];
		for (int j = 0; j < N[i]; j++)
		{
			data[i][j] = SA.data[i][j];
			P[i][j] = SA.P[i][j];
		}
	}
}

surfArray::surfArray(int nObj, Form **Objs)
{
	this->N = new int[nObj];
	data = new Vector<double> *[nObj];
	P = new Vector<double> *[nObj];
	this->nObjs = nObj;

	for (int i = 0; i < nObj; i++)
	{
		switch (Objs[i]->type)
		{
		   case SURFACE:  
						  surface *s = (surface *)Objs[i];
			              N[i] = s->anzp;
						  data[i] = new Vector<double>[N[i]];
			          
			              P[i] = new Vector<double>[N[i]];
						  for (int j = 0; j < N[i]; j++)
						  {							  
							  P[i][j] = (s->S[j][0] + s->S[j][1] + s->S[j][2]) / 3.0;
						  }
			   break;
		}
	}
}

surfArray& surfArray::operator =(const surfArray& SA)
{
	nObjs = SA.nObjs;
	data = new Vector<double> *[nObjs];
	P = new Vector<double>*[nObjs];
	N = new int[nObjs];
	for (int i = 0; i < nObjs; i++)
	{
		N[i] = SA.N[i];
		data[i] = new Vector<double>[N[i]];
		P[i] = new Vector<double>[N[i]];
		for (int j = 0; j < N[i]; j++)
		{
			data[i][j] = SA.data[i][j];
			P[i][j] = SA.P[i][j];
		}
	}
	return *this;
}

void surfArray::addValue(int objIndex, surface *S, int i, Vector<double> v)
{
	// cout << "addvalue: objindex=" << objIndex << "    i=" << i << endl;
	data[objIndex][i] += v;
}

void surfArray::genEllipsoidGitter(int currIndex, FormEllipsoid *E)
{
 
}

void surfArray::clearData()
{
	if (nObjs > 0)
	{
		for (int i = 0; i < nObjs; i++)
		{
			delete[] P[i];
			delete[] data[i];
		}
		delete[] data;
		delete[] P;
	}
}

surfArray::~surfArray()
{
}

