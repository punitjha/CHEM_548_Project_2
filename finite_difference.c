#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <armadillo>
#include <>



using namespace std;

//extern "C" void dstevx_(const char *JOBZ,const char *RANGE, int *N, double *D, double *E,const  double *VL, const  double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *IWORK, int *IFAIL, int *INFO);       
/* ***********************************************************************************************************/	

int main()
{
	int points=200, states=7; //no. of mesh points and no. of states to be obtained
	int i,n,d=1;                    //variables to be used in looping
	int width;
	int xxx, m, iwork[5*points], ifail[points], info;
	double potent[points], pmin, diag[points], subdiag[points-1], energies[points], eigenvectors[points][states], scale_graph, abstol, work[5*points], dlamch;
	double vl=0,vu=1;
//	char jobz='v',range='i';	

	double *diag = new double[points];

	cout<<"Enter the widht of the box"<<endl;
	cin>>width;

	xxx=width/points-1; //caluculating the distane between the mesh points
	cout<<"Particle in a box potential"<<endl;
	for(i=1;i<=points;i++)
	{
		potent[i]=0;
	}
	/* we are now creating the diagonal and the subdiagonal matrix elements */

	for (i=i;i<points;i++ )
	{
		diag[i]=potent[i]+(2/xxx^2)/2;
		subdiag[i]=-(1/xxx^2)/2;

	}
	// we are now calling the Lapack Package written in FORTRAN to diagonalize the matrix and find the eigenvalues and the eigenvectors
	abstol=2.0*dlamch;
	
//	 dstevx_(&jobz , &range, &points, diag, & *subdiag , &vl, &vu, &d, &states, &abstol, &m, & *energies, & **eigenvectors, &points, & *work, & *iwork, & *ifail, &info);
	if(info == 0)
	cout<<"diagnolization failed"<<endl;
	cout<<"the eigenvalues are :"<<endl;
 	for (int i=0; i<=states; i++ )
	{
		cout<<energies[i]<<endl;
	}



}
