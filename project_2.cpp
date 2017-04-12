// to compile this program rum g++  project_2.cpp -01 -larmadillo llapack -lblas  -Wall 
// the -Wall attribute is used to show all the errors you may possible get in your score code
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
	double points=400;
	double width=10, n;
	double  meshx;
	meshx=width/(points-1);	
	arma::vec potent1= arma::zeros<arma::vec>(points);
	arma::vec potent2= arma::zeros<arma::vec>(points);
	arma::vec diagonal1= arma::zeros<arma::vec>(points);
	arma::vec subdiag1= arma::zeros<arma::vec>(points-1);
	arma::vec diagonal2= arma::zeros<arma::vec>(points);
	arma::vec subdiag2= arma::zeros<arma::vec>(points-1);
	for( double i=0; i<potent1.n_elem; ++i )
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
	arma::mat Matrix1(points, points, arma::fill::zeros);
	Matrix1.diag()=diagonal1;
	Matrix1.diag(-1)=subdiag1;
	Matrix1.diag(1)=subdiag1;
	arma::mat Matrix2(points, points, arma::fill::zeros);
	Matrix2.diag()=diagonal2;
	Matrix2.diag(-1)=subdiag2;
	Matrix2.diag(1)=subdiag2;
	arma::vec eigenvalues1;
	arma::mat eigenvectors1;
	arma::eig_sym(eigenvalues1, eigenvectors1, Matrix1);
	arma::vec eigenvalues2;
	arma::mat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Matrix2);
	eigenvalues1.save("eigen1.txt",arma::raw_ascii);
	eigenvectors1.save("eigenvectors1.txt",arma::raw_ascii);
	eigenvalues2.save("eigen2.txt",arma::raw_ascii);
	eigenvectors2.save("eigenvectors2.txt",arma::raw_ascii);




//	calculating the Frank-Condon factors
	arma::mat frank = eigenvectors1.t()*eigenvectors2;
	arma::vec frank2 = arma::zeros<arma::vec>(frank.n_cols);

	for(int i=0; i< frank.n_cols; i++)
	{
		frank2(i)=frank(0,i)*frank(0,i);  //it showed segmentation fault if frank2 was not initialized as matrix full of zeros
	}
	frank2.save("frank.txt",arma::raw_ascii);



//	the lines below are for the second part of the homework assigment
	arma::mat square=eigenvectors2.t()*eigenvectors2;
	arma::vec new_wave = eigenvectors1.col(0);
	arma::vec ddd=square.diag();
	arma::vec newwave1 = ddd%new_wave;
	newwave1.save("newwave.txt", arma::raw_ascii);




//	the third part of the assingment regarding the time evolution of the wavepacket formed in the excited state
	complex<double> I(0.0,1.0);
	arma::cx_mat time_evolve=arma::zeros<arma::cx_mat>(400,400);
	double t=0;
	for(int a=0; a <time_evolve.n_cols; ++a)
	{
		for(int i=0; i<time_evolve.n_rows; ++i)
		{
			complex<double> coef =exp(-I*t*eigenvalues2(i));
			time_evolve(i,a)=newwave1[i]*coef;
		}	
		t+=0.005;	
	}
	fstream myfile("time1.txt",fstream::out |fstream::trunc );
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





//	arma::cx_vec corr_store=arma::zeros<arma::cx_vec>(400);
	int time_steps = 100000;
	arma::cx_vec corr=arma::zeros<arma::cx_vec>(time_steps);
	arma::vec time=arma::zeros<arma::vec>(time_steps);
	int loop=0;
	//t=0;
	double final_t = 100;
	double initial_t = 0;
	double dt = (final_t - initial_t) / (time_steps - 1);
	while(loop < time_steps)
	{	
		t = initial_t + dt * loop;
		time(loop)=t;
		for(int row=0; row <corr_store.n_rows; ++row)
		{
			complex<double> coef =exp(-I*t*eigenvalues2(row));
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
	arma::cx_vec four= arma::fft(corr);
	fstream myfile2("FT1.txt",fstream::out |fstream::trunc );
	for (int j=0; j<four.n_elem; ++j)
	{
		if (myfile2.is_open())
		{
			myfile2<<abs(four(j))<<endl;
		}
	}
	four.save("FT.txt",arma::raw_ascii);
	corr.save("correlation1.txt",arma::raw_ascii);
	time_evolve.save("time_evolve.txt",arma::raw_ascii);



// 	Question no. 6 of the project -- evolving the wavefunction using the finite difference algorithym 
	
	double it=0; //initial time
	double ft=100; // final time
	double ddt=10000; // number of steps 
	double increment = (ft-it)/(ddt-1); // step size
	arma::cx_mat wavepro = arma::zeros<arma::cx_mat> (newwave1.n_rows,ddt);
	arma::cx_vec store= newwave1*I; 
	arma::vec time2 = arma::zeros<arma::vec> (ddt); //vector to store the time as it increases
	double tim=it; // tim is the time looping variable
	for (int col=0; col<wavepro.n_cols; col++)
	{
		store=((Matrix2*store)*(tim/I))+store; // propagates the wave by the Hamiltonian matrix Matrix2
		tim+=increment; //increments the time
		time2(col)=it+tim*col; // storing the time in the array/vec
		store = store/sqrt(arma::dot(store,store));
		store = arma::normalise(store);
		for (int row=0; row<wavepro.n_rows; row++)
		{
			wavepro(row,col)=store(row);	
		}
	}
	time2.save("time2.txt",arma::raw_ascii);
	wavepro.save("wavepro.txt",arma::raw_ascii);	
	fstream myfile77("wavepro1.txt",fstream::out |fstream::trunc);
	for (int i=0; i<wavepro.n_rows; i++)
	{
		for (int j=0; j<wavepro.n_cols; ++j)
		{
			if (myfile.is_open())
			{
					myfile77<<std::real(wavepro(i,j))<<' '<<std::imag(wavepro(i,j))<<' ';
			}
		}
		myfile77<<endl;
	}

//	calculating the correlation function for wave propagated by the finite difference method
	arma::cx_vec corr2 = arma::zeros<arma::cx_vec> (wavepro.n_cols);
	for(int i=0; i<wavepro.n_cols; i++)
	{
		complex <double> sum = (0,0);
		for (int j=0; j<wavepro.n_rows; j++)
		{
		sum=sum+wavepro(j,i)*newwave1(j);
		}
		corr2(i)=sum;
	}
	corr2.print();
	fstream myfile88("corr2.txt", fstream::out | fstream::trunc);
	for(int i =0; i<corr2.n_rows; i++)
	{
		if (myfile88.is_open())
		{
			myfile88<<std::real(corr2(i))<<' '<<std::imag(corr2(i))<<endl;
		}
	}

// Doing the FT of the correlation function above
	
	arma::cx_vec four2= arma::fft(corr2);
	fstream myfile99("FT2.txt",fstream::out |fstream::trunc );
	for (int j=0; j<four2.n_elem; ++j)
	{
		if (myfile99.is_open())
		{
			myfile99<<abs(four2(j))<<endl;
		}
	}
}
