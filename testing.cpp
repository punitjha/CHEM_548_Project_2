#include <armadillo>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>

using namespace std;
int main()
{
	arma::mat M(3,3, arma::fill::randu);
	arma::mat B(3,3, arma::fill::randu);
	M.print();
	cout<<endl;
	B.print();
	cout<<endl;
	arma::mat C = M.t()*B;
	C.print();
	arma::cx_mat time_evolve=arma::zeros<arma::cx_mat>(10,10);
	float t=0;
	for(int i=0; i<10; ++i)
		{
			for(int a=0; a < 10; ++a)
					{
						time_evolve(i,a)=0;
						cout<<time_evolve(i,a)<<endl; 
					}	
		}
	arma::cx_fmat corr = arma::zeros<arma::cx_fmat> (26,10);		
	arma::vec nos = arma::linspace(1,5,6);
	arma::vec add = nos*10;
	complex <double> i = (0,1);
	arma::cx_vec no=nos*i;
	no.print();
	add.print("This is add vector");
	cout<<arma::dot(nos,nos);
	arma::vec co= arma::zeros<arma::vec> (100000);
	for (int i=0; i<co.n_rows; i++)
	{
		co(i)=sin((2*3.14159*i)/180)+sin((10*2*3.14159*i)/180);
	}
	co.save("sine.txt",arma::raw_ascii);
	arma::cx_vec corf= arma::fft(co);
	fstream myfile99("FT3.txt",fstream::out| fstream::trunc);
    	for (int j=0; j<corf.n_elem; ++j)
		{
			if (myfile99.is_open())
				{
					myfile99<<abs(corf(j))<<endl;
				}
		}

		

}
