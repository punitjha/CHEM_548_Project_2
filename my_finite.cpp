#include <armadillo> // the Armadillo linear algebra package of c++
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include <complex> //for using complex numbers in the program
using namespace std;
int main()
{
	float points=400;
	float width=10, n;
	float  meshx;

	//cout<<"Enter the number of points in the mesh"<<endl;
	//cin>>points;
	//cout<<"Enter the width of the box"<<endl;
	//cin>>width;
	meshx=width/(points-1);	
	arma::fvec potent1= arma::zeros<arma::fvec>(points);
	arma::fvec potent2= arma::zeros<arma::fvec>(points);
	arma::fvec diagonal1= arma::zeros<arma::fvec>(points);
	arma::fvec subdiag1= arma::zeros<arma::fvec>(points-1);
	arma::fvec diagonal2= arma::zeros<arma::fvec>(points);
	arma::fvec subdiag2= arma::zeros<arma::fvec>(points-1);
	for( float i=0; i<potent1.n_elem; ++i )
		{
			potent1[i]=300*(pow((1-exp(-1.2*(i-(points/3))*meshx)),2.0));
			potent2[i]= 295*(pow((1-exp(-0.8*(i-(points/2.56))*meshx)),2.0))+184;
		}
	potent1.save("potent1.txt",arma::raw_ascii);
	potent2.save("potent2.txt",arma::raw_ascii);

	for ( int i=0; i<diagonal1.n_elem; ++i)
		{
			diagonal1[i]=potent1[i]+(2.0/(meshx*meshx))/2.0;//notice that the power operator in c++ is not ^ or ** you have to use a function pow ()
			diagonal2[i]=potent2[i]+(2.0/(meshx*meshx))/2.0;
		}
	for ( int i=0; i<subdiag1.n_elem; ++i)
		{
			subdiag1[i]=(-1.0/(meshx*meshx))/2;
			subdiag2[i]=(-1.0/(meshx*meshx))/2;
		}
	arma::fmat Matrix1(points, points, arma::fill::zeros);
	Matrix1.diag()=diagonal1;
	Matrix1.diag(-1)=subdiag1;
	Matrix1.diag(1)=subdiag1;
	arma::fmat Matrix2(points, points, arma::fill::zeros);
	Matrix2.diag()=diagonal2;
	Matrix2.diag(-1)=subdiag2;
	Matrix2.diag(1)=subdiag2;
	
	arma::fvec eigenvalues1;
	arma::fmat eigenvectors1;
	arma::eig_sym(eigenvalues1, eigenvectors1, Matrix1);
	arma::fvec eigenvalues2;
	arma::fmat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Matrix2);
		
	//	eigenvectors.print("Eigenvectors");
	//	arma::uvec sort1 = arma::sort_index(eigenvalues);
	//	eigenvalues = eigenvalues(sort1);
	//	eigenvectors.print("Eigenvectors");
	eigenvalues1.save("eigen1.txt",arma::raw_ascii);
	eigenvectors1.save("eigenvectors1.txt",arma::raw_ascii);
	eigenvalues2.save("eigen2.txt",arma::raw_ascii);
	eigenvectors2.save("eigenvectors2.txt",arma::raw_ascii);
	



//	calculating the Frank-Condon factors
	arma::fmat frank = eigenvectors1.t()*eigenvectors2;
	arma::fvec frank2 = arma::zeros<arma::fvec>(frank.n_cols);

	for(int i=0; i< frank.n_cols; i++)
		{
			frank2(i)=frank(0,i)*frank(0,i);  //it showed segmentation fault if frank2 was not initialized as matrix full of zeros
		}
	
	frank2.save("frank.txt",arma::raw_ascii);


//	the lines below are for the second part of the homework assigment

	arma::fmat square=eigenvectors2.t()*eigenvectors2;
	arma::fvec new_wave = eigenvectors1.col(0);
	arma::fvec ddd=square.diag();
	arma::fvec newwave1 = ddd%new_wave;
	newwave1.save("newwave.txt", arma::raw_ascii);

	

//	the third part of the assingment regarding the time evolution of the wavepacket formed in the excited state
	complex<float> I(0.0,1.0);
	arma::cx_mat time_evolve=arma::zeros<arma::cx_mat>(400,400);
	float t=0;
		for(int a=0; a <time_evolve.n_cols; ++a)
		{
			for(int i=0; i<time_evolve.n_rows; ++i)
					{
						complex<float> coef =exp(-I*t*eigenvalues2(i));
						time_evolve(i,a)=newwave1[i]*coef;
					}	
			
			t+=0.005;	
		}
	fstream myfile("time.txt",fstream::out |fstream::trunc );
	for (int i=0; i<time_evolve.n_rows; i++)
		{
			for (int j=0; j<time_evolve.n_cols; ++j)
				{
					if (myfile.is_open())
						{
							myfile<<std::real(time_evolve(i,j))<<' '<<std::imag(time_evolve(i,j))<<' ';
						}
				}
				
			myfile<<endl;
		}
//	this part of the program correspond to the 4th question of the project
	
//	arma::cx_fmat corr=arma::zeros<arma::cx_fmat>(400,400);
//	t=0;
//	for(int c=0; c <corr.n_cols; ++c)
//		{
//			for(int row=0; row <corr.n_rows; ++row)
//				{
//					complex<float> coef =exp(-I*t*eigenvalues2(row));
//					corr(row,c)=frank2(row)*coef; 
//				}
//			t+=0.006;	
//		}
//	fstream myfile1("correlation.txt",fstream::out |fstream::trunc );
//		for(int row=0;  row<corr.n_rows; ++row)
//			{
//				for(int col=0; col <corr.n_cols; ++col)
//					{
//						if (myfile1.is_open())
//							{
//								myfile1<<std::real(corr(row,col))<<' '<<std::imag(corr(row,col))<<' ';
//							}
//					}
//				myfile1<<endl;
//			}
	arma::cx_fvec corr_store=arma::zeros<arma::cx_fvec>(400);
	int time_steps = 16384;
	arma::cx_fvec corr=arma::zeros<arma::cx_fvec>(time_steps);
	arma::fvec time=arma::zeros<arma::fvec>(time_steps);
	int loop=0;
	//t=0;
	float final_t = 100;
	float initial_t = 0;
	float dt = (final_t - initial_t) / (time_steps - 1);
	while(loop < time_steps)
	{	
		t = initial_t + dt * loop;
		time(loop)=t;
		for(int row=0; row <corr_store.n_rows; ++row)
		{
			complex<float> coef =exp(-I*t*eigenvalues2(row));
			corr_store(row)=frank2(row)*coef;
		}
		corr(loop)=arma::sum(corr_store);
		loop++;
	}
	time.save("time.txt",arma::raw_ascii);
	fstream myfile1("correlation.txt",fstream::out |fstream::trunc );
	for(int row=0; row <corr.n_rows; ++row)
	{
		if (myfile1.is_open())
		{
			myfile1<<std::real(corr(row))<<' '<<std::imag(corr(row))<<' ';
			myfile1<<endl;

		}
	}

//	carrying out the Fourier transform of the correlation function.

	arma::cx_fvec four= arma::fft(corr);
	fstream myfile2("FT1.txt",fstream::out |fstream::trunc );
	for (int j=0; j<four.n_elem; ++j)
	{
		if (myfile2.is_open())
		{
			myfile2<<std::real(four(j))<<' '<<std::imag(four(j))<<endl;
		}
	}
	four.save("FT.txt",arma::raw_ascii);
	corr.save("correlation1.txt",arma::raw_ascii);
	time_evolve.save("time_evolve.txt",arma::raw_ascii);
	
// 	Question no. 6 of the project -- evolving the wavefunction using the finite difference algorithym 
	

	
}
	
	
				
