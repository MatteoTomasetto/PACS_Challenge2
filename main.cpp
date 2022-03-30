#include <iostream>
#include <cmath>
#include <memory>
#include <string>
#include <limits>

#include "ZeroFun.hpp"
#include "GetPot"


// The function for which we want the zero
double myfun(const double & x)
{
	return 0.5 - std::exp(M_PI * x);
}


// Derivative of the function for which we want the zero (needed for Newton method)
double mydfun(const double & x)
{
	return - M_PI * std::exp(M_PI * x);
}


// Compute the zero of a function with the method in input
int main(int argc, char **argv)
{
	// Read in input from command line the datafile name and the method to use
	GetPot cl(argc, argv);
	
	const std::string filename = cl.follow("data", 2, "-f", "--file");
	const std::string method_name = cl("method", "Bisection"); 

	// Read constant parameter in input from datafile
	GetPot datafile(filename.c_str());
	const std::string section = "ZeroFun/";
	const std::string subsection{method_name + "/"};
	
	const SolverTraits::Real		tol = datafile((section + "tol").data(), 1.e-5);
	const SolverTraits::Uint		maxIt = datafile((section + "maxIt").data(), 200);
	const SolverTraits::Real		tola = datafile((section + subsection + "tola").data(), 1.e-10);
	const SolverTraits::InputType	x = datafile((section + subsection + "x").data(), 0.0);
	const SolverTraits::InputType	a = datafile((section + subsection + "a").data(), 0.0);
	const SolverTraits::InputType	b = datafile((section + subsection + "b").data(), 1.0);
	const SolverTraits::Interval	interval = std::make_pair(a, b);
	const SolverTraits::InputType	h = datafile((section + subsection + "h").data(), 1.e-2); 
	const SolverTraits::InputType	sol_ex = datafile((section + "sol_ex").data(),
													   std::numeric_limits<SolverTraits::InputType>::quiet_NaN()); 
	const SolverTraits::Uint		maxIt_interval = datafile((section + subsection + "maxIt_interval").data(), 200);
	const SolverTraits::InputType	h_interval = datafile((section + subsection + "h_interval").data(), 0.1);
	
	// Save the parameters in a struct
	Parameters param = {myfun, tol, maxIt, tola, interval, x, h, maxIt_interval, h_interval, mydfun};
	
	// Initialize the solver factory
	SolverFactory solver(param);
		
	// Solve the problem
	std::unique_ptr<SolverBase> my_ptr = solver(method_name);

	std::cout << "Finding the zero with " << method_name << " method" << std::endl;
	SolverBase::SolverOutput res = my_ptr -> solve();
		
	if(res.second)
	{
		std::cout << "The zero is " << res.first << std::endl;
		if(!std::isnan(sol_ex))
		{
			std::cout << "The exact zero is " << sol_ex << std::endl;
			std::cout << "Approximation error " << std::abs(sol_ex - res.first) << std::endl;
		}
	}
	else
		std::cout << "Zero not found! Try to change the parameters or the initial values" << std::endl;
	
	return 0;
}
