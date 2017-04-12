//Lapack Test Code

#include <iostream>
#include <vector>



using namespace std;


/*Lets say we want to solve the following linear system:

x + y = 2

x - y = 0

The solution is (x,y) = (1,1)*/

//this program solves a system of linier equations using the Lapack library
//the commmanf to compile and link the C++ program to the libraries can be found here : http://www.hoffman2.idre.ucla.edu/lapack/#cpp


// the line below declares each LAPACK function that is used as being extern
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

int main()
{

    char trans = 'N';
    int dim = 2;    
    int nrhs = 1;
    int LDA = dim;
    int LDB = dim;
    int info;

    vector<double> a, b;

    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(-1);

    b.push_back(2);
    b.push_back(0);

    int ipiv[3];
    
    dgetrs_(&trans, &dim, &nrhs, & *a.begin(), &LDA, ipiv, & *b.begin(), &LDB, &info);


    std::cout << "solution is:";    
    std::cout << "[" << b[0] << ", " << b[1] << ", " << "]" << std::endl;
    std::cout << "Info = " << info << std::endl; 

    return(0);
}













