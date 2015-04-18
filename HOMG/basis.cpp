#include "StdAfx.h"
#include "basis.h"


static matrix<double> polynomial(std::vector<double>& x,double alpha,double beta,int N)
{
  /*
	function [P] = basis::polynomial(x,alpha,beta,N)
    Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
             (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
    Note   : They are normalized to be orthonormal.

	x: denotes an array of points (float or double based on the required precision)
	alpha, beta: denotes the constants in the Jacobi polynomial
	N: denotes the order of the Jacobi polynomial

 */ 
	const int dim_1=N+1;
	const int dim_2=x.size();
	matrix<double> PL=zero_matrix<double>(dim_1,dim_2);
	matrix<double> P1=zero_matrix<double>(dim_1,dim_2);
	matrix<double> P2=zero_matrix<double>(dim_2,1);
	
	double gamma0=(pow(2,(alpha+beta+1))/(alpha+beta+1))*tgamma(alpha+1)*((tgamma(beta+1))/(tgamma(alpha+beta+1)));
	
	for(int i=0;i<dim_2;i++)
		PL(1,i)=1.0/sqrt(gamma0);

	if(N==0)
	{
		for(int i=0;i<dim_1;i++)
		 for(int j=0;j<dim_2;j++)
			 P1(j,i)=PL(i,j);

		return P1;
	}

	double gamma1=(alpha+1)*((beta+1)/(alpha+beta+3))*gamma0;
	for(int i=0;i<dim_2;i++)
		PL(2,i)=((alpha+beta+2)*x[i]/2 + (alpha-beta)/2)/sqrt(gamma1);

	if(N==1)
	{
		for(int i=0;i<dim_2;i++)
			P2(i,0)=PL(N+1,i);

		return P2;
	}
	
	//Repeat value in recurrence
	double aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

	for(int i=0;i<dim_1-2;i++)
	{
		double h1=2*i+alpha+beta;
		double anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
		double bnew = -(pow(alpha,2)-pow(beta,2))/h1/(h1+2);
		
		for(int k=0;k<dim_2;k++)
			PL(i+2,k)=1/anew*( -aold*PL(i,k) + (x[k]-bnew)*PL(i+1,k));
		
		aold=anew;
	}

	for(int i=0;i<dim_2;i++)
		P2(i,0)=PL(N+1,i);


}


static matrix<double> gradient(std::vector<double>& r,double alpha, double beta, int N)
{
	/*
	 function [dP] = basis::gradient(r, alpha, beta, N);
      Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
                at points r for order N and returns dP[1:length(r))]
	*/

	matrix<double> dP=zero_matrix<double>(r.size(),1);
	matrix <double> poly;
	if(N!=0)
	{
		poly=polynomial(r,alpha+1,beta+1,N-1);
		for(int i=0;i<dP.size1();i++)
			dP(i,0)=sqrt(N*(N+alpha+beta+1))*poly(i,1);
	}

	return dP;


}


static void gauss(double alpha,double beta,int N, /*out parameter list ==>*/ matrix<double>&x, matrix<double>& w)
{
	/*
	  function [x,w] = basis::gauss(alpha,beta,N)
             Purpose: Compute the N'th order Gauss quadrature points, x,
                      and weights, w, associated with the Jacobi
                      polynomial, of type (alpha,beta) > -1 ( <> -0.5).
	*/
	
	// clear the output parameters. 
	x=zero_matrix<double>(N+1,1);
	w=zero_matrix<double>(N+1,1);
	if(N==0)
	{   
		double temp=(alpha-beta)/(alpha+beta+2);
		x(1,0)=temp;
		w(1,0)=2;
		return;	
	}

	//% Form symmetric matrix from recurrence.
	matrix<double> h1=zero_matrix<double>(1,N+1);
	matrix<double> J=zero_matrix<double>(N+1,N+1);
	matrix<double> J1=zero_matrix<double>(N+1,N+1);
	for(int i=0;i<h1.size2();i++)
		h1(0,i)=2*i+alpha+beta;

	for(int i=0;i<J.size1();i++)
	{
		for(int j=0;j<J.size2();j++)
		{
			if(i==j)
			{
				J(i,j)=-1/2*(pow(alpha,2)-pow(beta,2))/(h1(0,i)+2)/h1(0,i);
			}
			if(i+1==j)
			{
				J(i+1,j)=(2/(h1(0,j)+2))*sqrt((float)j)*(((j)+alpha+beta)*((j)+alpha)*((j)+beta)/(h1(0,j)+1)/(h1(0,j)+3));
			}

		}
	}
	
	if(alpha+beta<10*std::numeric_limits<double>::epsilon())
	{
		J(0,0)=0.0;	
	}

	for(int i=0;i<J.size1();i++)
		for(int j=0;j<J.size2();j++)
			J1(i,j)=J(i,j)+J(j,i);

	J=J1;

	// Need to find eigen values and eigen vectors of J and creat x and w values. 


}

static void gll(double alpha,double beta,int N, matrix<double>&x, matrix<double>& w)
{
	/*
	 function [x] = basis::gll (alpha,beta,N)
             Purpose: Compute the N'th order Gauss Lobatto quadrature
                      points, x, associated with the Jacobi polynomial,
                      of type (alpha,beta) > -1 ( <> -0.5).
	*/

	x=zero_matrix<double>(N+1,1);
	w=zero_matrix<double>(N+1,1);
	if(N==1)
	{
		x(0,0)=-1.0;
		x(1,0)=1.0;

		w(0,0)=1.0;
		w(1,0)=1.0;

		return;
	}

	std::vector<double> xint;
	std::vector<double> wq;

	gauss(alpha+1,beta+1,N-2,xint,wq);

	
	for(int i=0;i<x.size1();i++)
	{
		if(i==0){
			x(i,0)=-1;
			continue;
		}

		x(i,0)=xint[i];

		if(i==N){
			x(i,0)=1;
			continue;
		}

	}
		
	w=polynomial(x,alpha,beta,N);

	

}