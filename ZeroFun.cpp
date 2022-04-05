#include <iostream>
#include <cmath>
#include <limits>
#include "ZeroFun.hpp"


/*!
 * This function tries to find an interval that brackets the zero of a
 * function f. It does so by sampling the value of f at points
 * generated starting from a given point
 *
 * f --> The function.
 * x1 --> initial point
 * h --> initial increment for the sampling
 * maxIt --> maximum number of iterations
 * It retruns a pair with the bracketing points and a bool which is true if number
 * of iterations not exceeded (bracket found)
 */

std::pair<SolverTraits::Interval, bool>
SolverBaseInterval::bracketInterval(InputType x1)
{
	
	constexpr SolverTraits::InputType	expandFactor = 1.5;
	double								direction = 1.0;
	SolverTraits::InputType				x2 = x1 + h_interval;
	SolverTraits::OutputType			y1 = f(x1);
	SolverTraits::OutputType			y2 = f(x2);
	unsigned int 						iter{0u};

	// Get the initial decrement direction
	while((y1 * y2 > 0) && (iter < maxIt_interval))
	{ 
		++iter;
		if(std::abs(y2) > std::abs(y1))
		{
			std::swap(y1, y2);
			std::swap(x1, x2);
		}

		direction = (x2 > x1) ? 1.0 : -1.0;
		x1 = x2;
		y1 = y2;
		x2 += direction * h_interval;
		y2 = f(x2);
		h_interval *= expandFactor;
	}
	
	// Swap to have the endpoints in the correct order
	if(x1 > x2)
		std::swap(x1, x2);
	
	if (iter < maxIt_interval)
		std::cout << "Bracket interval found: [" << x1 << ", " << x2 << "]" << std::endl;
		 
	return std::make_pair(Interval{x1, x2}, iter < maxIt_interval);
}


/*!
 * This function check if the interval in input brackets a zero of the function f stored 
 * in the calss. If not, it tries to change the interval with bracketInterval() 
 *
 * It retruns a pair with the bracketing interval and a bool which is true if it found the bracket interval
 */
 
std::pair<SolverTraits::Interval, bool>
SolverBaseInterval::CheckInterval()
{
	InputType		a{interval.first};
	InputType		b{interval.second};
	OutputType		ya = f(a);
	OutputType		yb = f(b);

	if(ya * yb > 0)
	{
		std::cout << "Function must change sign at the two end values...trying to find a proper interval" << std::endl;
		
		std::pair<Interval, bool> newinterval1 = bracketInterval(a);
		if (newinterval1.second == false)
		{
			std::pair<Interval, bool> newinterval2 = bracketInterval(b);
			if (newinterval2.second == false)
			{
				std::cout << "ERROR: unable to find a proper interval that brackets the zero" << std::endl;
				return std::make_pair(this -> interval, false);
			}
			
			else
				return std::make_pair(newinterval2.first, true);
		}

		else
			return std::make_pair(newinterval1.first, true);
	
	}
	
	else
		return std::make_pair(this -> interval, true);
}

/*!
 * Computes the zero of a function with the method of the bisection
 *
 * f --> The function
 * interval.fist --> First end of initial interval
 * interval.second --> Second end of initial interval
 * tol --> Tolerance
 * It returns the approximation of the zero of f and a status (true if converging)
 *
 */

SolverTraits::SolverOutput
Bisection::solve()
{
	std::pair<Interval, bool> check_interval = CheckInterval();
	
	if(check_interval.second == false)
		return std::make_pair(std::numeric_limits<InputType>::quiet_NaN(), false);
	
	interval = check_interval.first;
	InputType		a{interval.first};
	InputType		b{interval.second};
	OutputType		ya = f(a);
	OutputType		yb = f(b);
	
	if((ya * yb) == 0.0)
	{
		if(ya == 0.)
			return std::make_pair(a, true);
		else if(yb == 0.)
			return std::make_pair(b, true);
    };
	
	InputType		delta = b - a;
	unsigned int	iter{0u};
	OutputType		yc{ya};
	InputType		c{a};

	while(std::abs(delta) > 2 * tol && iter < maxIt)
	{
		++iter;
		c = (a + b) / 2.;
		yc = f(c);
		if(yc * ya < 0.0)
		{
			yb = yc;
			b = c;
		}
		else
		{
			ya = yc;
			a = c;
		}
		delta = b - a;
	}
	return std::make_pair((a + b) / 2., (iter < maxIt));
}


/*!
 * Computes the zero of a scalar function with the method of the Regula Falsi
 * Stop when the residual is below tolerance;
 *
 * f --> The function
 * interval.first --> First end of initial interval
 * interval.second --> Second end of initial interval
 * tol --> Tolerance (relative)
 * tola --> Tolerance (absolute)
 * It returns the approximation of the zero of f
 */

SolverTraits::SolverOutput
RegulaFalsi::solve()
{
	std::pair<Interval, bool> check_interval = CheckInterval();
	
	if(check_interval.second == false)
		return std::make_pair(std::numeric_limits<InputType>::quiet_NaN(), false);
	
	interval = check_interval.first;
	InputType				a{interval.first};
	InputType				b{interval.second};
	OutputType				ya = f(a);
	OutputType				yb = f(b);
	
	if((ya * yb) == 0.0)
	{
		if(ya == 0.)
			return std::make_pair(a, true);
		else if(yb == 0.)
			return std::make_pair(b, true);
    };
	
	InputType				delta = b - a;
	OutputType				resid0 = std::max(std::abs(ya), std::abs(yb));
	OutputType				yc{ya};
	InputType				c{a};
    OutputType				incr = std::numeric_limits<double>::max();
	constexpr OutputType	small = 10.0 * std::numeric_limits<double>::epsilon();
	unsigned int			iter{0u};
	
	while(std::abs(yc) > tol * resid0 + tola && incr > small && iter < maxIt)
	{
		++iter;
		double incra = -ya / (yb - ya);
		double incrb = 1. - incra;
		double incr = std::min(incra, incrb);

		if(std::max(incra, incrb) >= 1.0 || incr <= 0)
		{
			std::cout << "ERROR: Chord is failing" << std::endl;
			return std::make_pair(std::numeric_limits<InputType>::quiet_NaN(), false);
		}
		
		c = a + incra * delta;
		yc = f(c);
		
		if(yc * ya < 0.0)
		{
			yb = yc;
			b = c;
		}

		else
		{
			ya = yc;
			a = c;
		}
		
		delta = b - a;
    }
  return std::make_pair(c, (iter < maxIt));
}


/*!
 * Brent type search
 * If converging, it finds a zero with error below the given tolerance.
 *
 * f --> The functino 
 * interval.first --> The first end of a bracketing interval
 * interval.second --> The second end of a bracketing interval
 * tol --> Tolerance
 * maxIter --> Max number of iteration.
 * It returns the found approximated zero and a status flag (true if converged).
 *
 */

SolverTraits::SolverOutput
Brent::solve()
{
	std::pair<Interval, bool> check_interval = CheckInterval();
	
	if(check_interval.second == false)
		return std::make_pair(std::numeric_limits<InputType>::quiet_NaN(), false);
	
	interval = check_interval.first;
	InputType		a{interval.first};
	InputType		b{interval.second};
	OutputType		ya = f(a);
	OutputType		yb = f(b);

	if((ya * yb) == 0.0)
	{
		if(ya == 0.)
			return std::make_pair(a, true);
		else if(yb == 0.)
			return std::make_pair(b, true);
    };

	if(std::abs(ya) < std::abs(yb))
	{
		std::swap(a, b);
		std::swap(ya, yb);
	}

	InputType		c = a;
	InputType		d = c;
	OutputType		yc = ya;
	bool			mflag{true};
	InputType		s = b;
	OutputType		ys = yb;
	unsigned int 	iter{0u};
	
	do
	{
		if(ya != yc and yb != yc)
		{
			OutputType yab = ya - yb;
			OutputType yac = ya - yc;
			OutputType ycb = yc - yb;

			// Inverse quadratic interpolation
			s = a * ya * yc / (yab * yac) + b * ya * yc / (yab * ycb) - c * ya * yb / (yac * ycb);
		}

		else
		{
			// Secant
			s = b - yb * (b - a) / (yb - ya);
		}

		if(((s - 3 * (a + b) / 4) * (s - b) >= 0) or
		(mflag and (std::abs(s - b) >= 0.5 * std::abs(b - c))) or
		(!mflag and (std::abs(s - b) >= 0.5 * std::abs(c - d))) or
		(mflag and (std::abs(b - c) < tol)) or
		(!mflag and (std::abs(c - d) < tol)))
		{
			mflag = true;
			s = 0.5 * (a + b); // Back to bisection step
		}
		
		else
			mflag = false;

		ys = f(s);
		d = c;
		c = b;
		yc = yb;

		if(ya * ys < 0)
		{
			b = s;
			yb = ys;
		}

		else
		{
			a = s;
			ya = ys;
		}

		if(std::abs(ya) < std::abs(yb))
		{
			std::swap(a, b);
			std::swap(ya, yb);
		}
	} while(ys != 0. && std::abs(b - a) > tol && iter < maxIt);

return std::make_pair(s, (iter < maxIt));
}

/*!
 * Computes the zero of a scalar function with the method of the secant
 * It stops when |f(solution)| <= tol|f(initial_solution)| + tola
 *
 * f --> The function
 * x --> First point for computation of derivatives
 * y --> Second point for computing derivatives
 * tol --> relative tolerance
 * tola --> absolute tolerance
 * maxIt --> maximum number of iterations
 * It returns the approximation of the zero of f and a status (false if not
 * converging)
 *
 */
 
SolverTraits::SolverOutput
Secant::solve()
{
	InputType		a{interval.first};
	InputType		b{interval.second};
	OutputType		ya = f(a);
	OutputType		resid = std::abs(ya);
	InputType		c{a};
	unsigned int	iter{0u};
	Real			check = tol * resid + tola;
	bool			goOn = resid > check;
	
	while(goOn && iter < maxIt)
	{
		++iter;
		OutputType yb = f(b);
		c = a - ya * (b - a) / (yb - ya);
		OutputType yc = f(c);
		resid = std::abs(yc);
		goOn = resid > check;
		a = c;
		ya = yc;
	}

	return std::make_pair(c, (iter < maxIt));
}


/*!
 * Computes the zero of a scalar function with the Newton
 * It stops when |f(solution)| <= tol|f(initial_solution)| + tola
 *
 * f --> The function
 * df --> the derivative (it can be approximated through differences)
 * x --> Initial point
 * tol --> relative tolerance
 * tola --> absolute tolerance
 * maxIt --> maximum number of iterations
 * It returns the approximation of the zero of f and a status (false if not
 * converging)
 *
 */

SolverTraits::SolverOutput
Newton::solve()
{  
	InputType		a{x};
	OutputType		ya = f(a);
	OutputType		resid = std::abs(ya);
	unsigned int	iter{0u};
	Real			check = tol * resid + tola;
	bool			goOn = resid > check;
	
	while(goOn && iter < maxIt)
	{
		++iter;
		a += - ya/df(a);
		ya = f(a);
		resid = std::abs(ya);
		goOn = resid > check;
	}

	return std::make_pair(a, (iter < maxIt));
}



/*!
 * Computes the zero of a scalar function with the Newton
 * It stops when |f(solution)| <= tol|f(initial_solution)| + tola
 *
 * f --> The function
 * df --> the derivative (it can be approximated through differences)
 * x --> Initial point
 * tol --> relative tolerance
 * tola --> absolute tolerance
 * maxIt --> maximum number of iterations
 * It returns the approximation of the zero of f and a status (false if not
 * converging)
 *
 */

SolverTraits::SolverOutput
QuasiNewton::solve()
{  
	InputType		a{x};
	OutputType		ya = f(a);
	OutputType		resid = std::abs(ya);
	unsigned int	iter{0u};
	Real			check = tol * resid + tola;
	bool			goOn = resid > check;
		
	while(goOn && iter < maxIt)
	{
		++iter;
		a += - ya/df(a);
		ya = f(a);
		resid = std::abs(ya);
		goOn = resid > check;
	}

	return std::make_pair(a, (iter < maxIt));
}

