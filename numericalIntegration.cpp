/**
 * @file   numericalIntegration.cpp
 * @author HU Guanghui <gary@Vivian>
 * @date   Fri Mar 11 08:43:45 2022
 * 
 * @brief Numerical integration, including Trapezoidal rule, Midpoint
 * rule, Simpson rule.
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <omp.h>

/// Exact value of the integral
const double value = 0.7468241328124270253994674361318530053544996868126063290276544989;

/** 
 * The function exp(-x * x), whose antiderivative is given using error
 * function erf(x).
 * 
 * @param x input value
 * 
 * @return function value
 */
double f(double& x)
{
  return exp(- x * x);
}

/** 
 * main function
 * 
 * @param argc 
 * @param argv 
 * 
 * @return 
 */
int main(int argc, char* argv[])
{

  
  double b_pnt = atof(argv[1]);
  double e_pnt = atof(argv[2]);
  int n_interval = atoi(argv[3]);

  double h = (e_pnt - b_pnt)/n_interval;

  std::vector<double> f_vals(n_interval + 1);
#pragma omp parallel for
  for(int i = 0;i < n_interval + 1;++ i){
    double pnt = b_pnt + i * h;
    f_vals[i] = f(pnt);
  }

  std::vector<double> f_mid_vals(n_interval);
#pragma omp parallel for  
  for(int i = 0;i < n_interval;++ i){
    double pnt = b_pnt + (i + 0.5) * h;
    f_mid_vals[i] = f(pnt);
  }
  

//@{  Trapezoidal rule
  double start_time = omp_get_wtime();  
  double Trape = 0.;
#pragma omp parallel for reduction (+:Trape)  
  for(int i = 1;i < n_interval;++ i){
    Trape += f_vals[i];
  }
  std::cout << std::fixed << std::setprecision(15)
	    << "Trapezoidal rule gives: "
	    << h * (0.5 * f(b_pnt) + Trape + 0.5 * f(e_pnt))
	    << ", with error: "
	    << std::scientific << std::setprecision(4)
	    << fabs(h * (0.5 * f(b_pnt) + Trape + 0.5 * f(e_pnt)) - value)
	    << ", with CPU seconds: "
	    << omp_get_wtime() - start_time
	    << std::endl;
//@}

//@{  Midpoint rule
  start_time = omp_get_wtime();    
  double midRule = 0.;
#pragma omp parallel for reduction (+:midRule)    
  for(int i = 0;i < n_interval;++ i){
    midRule += f_mid_vals[i];
  }
  std::cout << std::fixed << std::setprecision(15)
	    << "   Midpoint rule gives: "
	    << midRule * h
	    << ", with error: "
	    << std::scientific << std::setprecision(4)
	    << fabs(midRule * h - value)    
	    << ", with CPU seconds: "
	    << omp_get_wtime() - start_time
	    << std::endl;
//@}

//@{ Simpson rule
  start_time = omp_get_wtime();      
  double Simpson = 0., Simpson_e = 0., Simpson_o = 0.;
#pragma omp parallel for reduction (+:Simpson_e)    
  for(int i = 1;i < n_interval;i = i + 2){
    Simpson_e += f_vals[i];
  }
#pragma omp parallel for reduction (+:Simpson_o)  
  for(int i = 2;i < n_interval;i = i + 2){
    Simpson_o += f_vals[i];
  }
  std::cout << std::fixed << std::setprecision(15)
	    << "    Simpson rule gives: "
	    << 1./3. * h * (f(b_pnt) + 4. * Simpson_e + 2. * Simpson_o + f(e_pnt))
	    << ", with error: "
	    << std::scientific << std::setprecision(4)
	    << fabs(1./3. * h * (f(b_pnt) + 4. * Simpson_e + 2. * Simpson_o + f(e_pnt)) - value)    
	    << ", with CPU seconds: "
	    << omp_get_wtime() - start_time
	    << std::endl;
//@}
  return 0;
}

/**
 * end of file
 * 
 */
