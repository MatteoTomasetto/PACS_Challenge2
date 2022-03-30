#ifndef HH__ZERO_FUN__HH
#define HH__ZERO_FUN__HH

#include <string>
#include <functional>
#include <memory>
#include <tuple>


// Traits with the types used in the child classes
struct SolverTraits
{
	using InputType = double;
	using OutputType = double;
	using FunType = std::function<OutputType (const InputType &)>;
	using Real = double;
	using Uint = unsigned int;
	using Interval = std::pair<InputType,InputType>;
	using SolverOutput = std::pair<InputType, bool>;
};


// Abstract base class for methods that find the zero of a function
class SolverBase: public SolverTraits
{
	public: SolverBase(FunType f_, const Real & tol_, const Uint & maxIt_) 
			: f(f_), tol(tol_), maxIt(maxIt_) {};
	
			SolverBase(FunType f_) : f(f_), tol(1.e-5), maxIt(200) {};
			
			virtual SolverOutput solve() = 0;
						
			inline void set_f(FunType f_){ f = f_; };
			
			inline FunType get_f() const { return f; };
			
			virtual ~SolverBase() = default;
			
	protected:  FunType f;
				const Real tol;
				const Uint maxIt;
};


// Abstract base class for methods that need a bracket interval
class SolverBaseInterval: public SolverBase
{
	public: SolverBaseInterval(FunType f_, const Real & tol_, const Uint & maxIt_, Interval interval_, 
			const Uint & maxIt_interval_, InputType h_interval_) 
			: SolverBase(f_, tol_, maxIt_), interval(interval_), maxIt_interval(maxIt_interval_), h_interval(h_interval_) {};
	
			SolverBaseInterval(FunType f_, Interval interval_) : SolverBase(f_), interval(interval_), 
																 maxIt_interval(200), h_interval(0.1) {};
			
			virtual SolverOutput solve() = 0;
			
			inline void set_interval(Interval interval_){ interval = interval_; };
			inline void set_h_interval(InputType h_interval_){ h_interval = h_interval_; };
			
			inline Interval get_interval() const { return interval; };
			inline InputType get_h_interval() const { return h_interval; };
			
			virtual ~SolverBaseInterval() = default;
			
	protected: Interval interval;
			   const Uint maxIt_interval;
			   InputType h_interval;
			   	
			   std::pair<Interval, bool> bracketInterval(InputType x1);
	
			   std::pair<SolverTraits::Interval, bool> CheckInterval();
};

// Bisection method
class Bisection final: public SolverBaseInterval
{

	public: Bisection(FunType f_, const Real & tol_, const Uint & maxIt_, Interval interval_, 
			const Uint maxIt_interval_, InputType h_interval_) 
			: SolverBaseInterval(f_, tol_, maxIt_, interval_, maxIt_interval_, h_interval_) {};
			
			Bisection(FunType f_, Interval & interval_) 
			: SolverBaseInterval(f_, interval_) {};
	
			SolverOutput solve() override;

};


// RegulaFalsi method
class RegulaFalsi final: public SolverBaseInterval
{

	public: RegulaFalsi(FunType f_, const Real & tol_, const Uint & maxIt_, const Real tola_, Interval interval_,
			const Uint maxIt_interval_, InputType h_interval_) 
			: SolverBaseInterval(f_, tol_, maxIt_, interval_, maxIt_interval_, h_interval_), tola(tola_) {};
			
			RegulaFalsi(FunType f_, Interval & interval_) 
			: SolverBaseInterval(f_, interval_), tola(1.e-10) {};
	
			SolverOutput solve() override;
				
	private: const Real tola;
};


// Brent method
class Brent final: public SolverBaseInterval
{

	public: Brent(FunType f_, const Real & tol_, const Uint & maxIt_, Interval interval_,
			const Uint maxIt_interval_, InputType h_interval_) 
			: SolverBaseInterval(f_, tol_, maxIt_, interval_, maxIt_interval_, h_interval_) {};
			
			Brent(FunType f_, Interval & interval_) 
			: SolverBaseInterval(f_, interval_) {};
	
			SolverOutput solve() override;

};


// Secant method
class Secant final: public SolverBase
{

	public: Secant(FunType f_, const Real & tol_, const Uint & maxIt_, const Real & tola_, const Interval & interval_) 
			: SolverBase(f_, tol_, maxIt_), tola(tola_), interval(interval_) {};
			 
	
			Secant(FunType f_, const Interval & interval_)
			: SolverBase(f_), tola(1.e-10), interval(interval_) {};
			
			SolverOutput solve() override;
			
			inline void set_interval(Interval interval_){ interval = interval_; };
			
			inline Interval get_interval() const { return interval; };
	

	private: const Real tola;
			 Interval interval;

};


// Newton method
class Newton : public SolverBase
{
	public: Newton(FunType f_, const Real & tol_, const Uint & maxIt_, const Real & tola_, const InputType & x_, 
			FunType df_) 
			:	SolverBase(f_, tol_, maxIt_), tola(tola_), x(x_), df(df_) {}; 
			
			Newton(FunType f_, const InputType & x_, FunType df_)
			:	SolverBase(f_), tola(1e-10), x(x_), df(df_) {};
			
			SolverOutput solve() override;
			
			inline void set_x(InputType x_){ x = x_; };
			inline void set_df(FunType df_){ df = df_; };
			
			inline InputType get_x() const { return x; };
			inline FunType get_df() const { return df; };
			
	protected: const Real tola;
			   InputType x;
		 	   FunType df;
		 	   
		 	   Newton(FunType f_, const Real & tol_, const Uint & maxIt_, const Real & tola_, const InputType & x_)
		 	   :	SolverBase(f_, tol_, maxIt_), tola(tola_), x(x_) {};
		 	   
		 	   Newton(FunType f_, const InputType & x_)
		 	   :	SolverBase(f_), tola(1e-10), x(x_) {};
};


// QuasiNewton method
class QuasiNewton final: public Newton
{
	public: QuasiNewton(FunType f_, const Real & tol_, const Uint & maxIt_, const Real & tola_, const InputType & x_, 
			const InputType & h_) 
			:	Newton(f_, tol_, maxIt_, tola_, x_), h(h_) 
			{
				df = [*this](const double & x){ return (f(x + h) - f(x - h)) / (2. * h); };
			};
				
			QuasiNewton(FunType f_, const InputType & x_, const InputType & h_)
			:	Newton(f_, x_), h(h_) 
			{
				df = [*this](const double & x){ return (f(x + h) - f(x - h)) / (2. * h); };
			};
						
			QuasiNewton(FunType f_, const InputType & x_)
			:	Newton(f_, x_), h(1.e-2)
			{
				df = [*this](const double & x){ return (f(x + h) - f(x - h)) / (2. * h); };
			};
			
			SolverOutput solve() override;
			
			inline void set_h(InputType h_){ h = h_; };
					
			inline InputType get_h() const { return h; };
						
	private: InputType h;
};



// Parameteres structure
struct Parameters
{	
	SolverTraits::FunType	f;
	SolverTraits::Real		tol = 1.e-5;
	SolverTraits::Uint		maxIt = 200;
	SolverTraits::Real		tola = 1.e-10;
	SolverTraits::Interval	interval;
	SolverTraits::InputType	x;
	SolverTraits::InputType	h = 1.e-2;
	SolverTraits::Uint		maxIt_interval = 200;
	SolverTraits::InputType	h_interval = 0.1;
	SolverTraits::FunType	df;
};


// SolverFactory to retrieve a pointer to a class method object
class SolverFactory
{

	public: SolverFactory(Parameters & param_) : param(param_) {};
	
			std::unique_ptr<SolverBase>	operator()(std::string solver_name) const;
			
			void set_param(Parameters & param_);
			
	private: Parameters param;
};

#endif
