#ifndef _PROBLEM_H
#define _PROBLEM_H

//#include "cec09.h"

//// Toolkit includes. //////////////////////////////////////////////////////

#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"
#define PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;




// *********************** CEC 2009 ************************************


void UF1(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
	}

	void UF2(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 0) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
	}

	void UF3(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
	}

	void UF4(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0*sum1 / (double)count1;
		f[1] = 1.0 - x[0]*x[0]	+ 2.0*sum2 / (double)count2;
	}

	void UF5(double *x, double *f, const unsigned int nx)
	{
		 unsigned int j, count1, count2;
                double sum1, sum2, yj, hj, N, E;
                sum1   = sum2   = 0.0;
                count1 = count2 = 0;
                N = 10.0; E = 0.1;
                for(j = 2; j <= nx; j++)
                {
                        yj = x[j-1]-sin(6.0*M_PI*x[0]+j*M_PI/nx);
                        hj = 2.0*yj*yj - cos(4.0*M_PI*yj) + 1.0;
                        if (j % 2 == 0)
                        {
                                sum2  += hj;
                                count2++;
                        }
                        else
                        {
                                sum1  += hj;
                                count1++;
                        }
                }
                hj = (0.5/N + E)*fabs(sin(2.0*N*M_PI*x[0]));
                f[0] = x[0]           + hj + 2.0*sum1 / (double)count1;
                f[1] = 1.0 - x[0] + hj + 2.0*sum2 / (double)count2;

	}

	void UF6(double *x, double *f, const unsigned int nx)
	{
  		unsigned int j, count1, count2;
                double sum1, sum2, prod1, prod2, yj, hj, pj, N, E;
                N = 2.0; E = 0.1;
                sum1   = sum2   = 0.0;
                count1 = count2 = 0;
                prod1  = prod2  = 1.0;
                for(j = 2; j <= nx; j++)
                {
                        yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
                        pj = cos(20.0*yj*PI/sqrt(j+0.0));
                        if (j % 2 == 0)
                        {
                                sum2  += yj*yj;
                                prod2 *= pj;
                                count2++;
                        }
                        else
                        {
                                sum1  += yj*yj;
                                prod1 *= pj;
                                count1++;
                        }
                }
                hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*x[0]);
                if(hj<0.0) hj = 0.0;
                f[0] = x[0]           + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
                f[1] = 1.0 - x[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;

	}

	void UF7(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x[0],0.2);
		f[0] = yj	    + 2.0*sum1 / (double)count1;
		f[1] = 1.0 - yj + 2.0*sum2 / (double)count2;
	}

	void UF8(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
	}

	void UF9(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, E;
		
		E = 0.1;
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		yj = (0.5+E)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
		if(yj<0.0) yj = 0.0;
		f[0] = 0.5*(yj + 2*x[0])*x[1]		+ 2.0*sum1 / (double)count1;
		f[1] = 0.5*(yj - 2*x[0] + 2.0)*x[1] + 2.0*sum2 / (double)count2;
		f[2] = 1.0 - x[1]                   + 2.0*sum3 / (double)count3;
	}

	void UF10(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
			if(j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += hj;
				count2++;
			}
			else
			{
				sum3  += hj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
	}




void CEC09_F1(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	UF1(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F2(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	UF2(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F3(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	UF3(&(*(X.begin())), &(*(F.begin())), X.size());
}
void CEC09_F4(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	UF4(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	UF5(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F6(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	UF6(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F7(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	UF7(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F8(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	UF8(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F9(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	UF9(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F10(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	UF10(&(*(XX.begin())), &(*(F.begin())), X.size());
}


///////////////////////////////////////wfg problems...
void wfg1(std::vector<double> &F, std::vector<double> &X)
{
	 //---- Get the function name.

	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG1( X, k, M );
}
void wfg2(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG2( X, k, M );
}
void wfg3(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG3( X, k, M );
}
void wfg4(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG4( X, k, M );
}
void wfg5(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG5( X, k, M );
}
void wfg6(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG6( X, k, M );
}
void wfg7(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG7( X, k, M );
}
void wfg8(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG8( X, k, M );
}
void wfg9(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG9( X, k, M );
}
void dtlz1(std::vector<double> &F, std::vector<double> &X)
{

  int k = X.size()- F.size() + 1;

  double g = 0.0 ;

  for (int i = X.size() - k; i < X.size(); i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100 * (k + g);
  for (int i = 0; i < F.size(); i++)
    F[i] = (1.0 + g) * 0.5;

  for (int i = 0; i < F.size(); i++){
    for (int j = 0; j < F.size() - (i + 1); j++)
      F[i] *= X[j];
      if (i != 0){
        int aux = F.size() - (i + 1);
        F[i] *= 1 - X[aux];
      } //if
  }//for

}
void dtlz2(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
   double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } //if
  } // for
}
void dtlz3(std::vector<double> &F, std::vector<double> &X)
{
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  int k = numberOfVariables_ - numberOfObjectives_ + 1;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100.0 * (k + g);
  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } // if
  } //for

}
void dtlz4(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++) {
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(pow(X[j],alpha)*(M_PI/2.0));
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(pow(X[aux],alpha)*(M_PI/2.0));
      } //if
  } // for

}
void dtlz5(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz6(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += pow(X[i],0.1);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz7(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += X[i] ;

  g = 1 + (9.0 * g)/k ;


  for (int i = 0; i < numberOfObjectives_ - 1; i++)
    F[i] = X[i] ;

  double h = 0.0 ;
  for (int i = 0; i < numberOfObjectives_ - 1; i++){
    h+=(F[i]/(1.0+g))*(1 + sin(3.0*M_PI*F[i])) ;
  } //for

  h = numberOfObjectives_ - h ;

  F[numberOfObjectives_ - 1] = (1+g)*h ;	
}

#endif
