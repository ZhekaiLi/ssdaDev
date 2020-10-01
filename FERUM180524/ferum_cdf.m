function P = ferum_cdf(type,x,varargin)

% Cumulative Density Function
%
%   P = ferum_cdf(type,x,varargin)
%
%   Evaluates the cumulative distribution function and returns the probability.
%
%   Output: - P           = cumulative density value 
%   Input:  - type        = probability distribution type (1: normal, 2: lognormal, ...)
%           - x           = 'Abscissa' value(s)
%           - varargin{1} = parameter #1 of the random variable
%           - varargin{2} = parameter #2 of the random variable (optional)
%           - varargin{3} = parameter #3 of the random variable (optional)
%           - varargin{4} = parameter #4 of the random variable (optional)

switch type
   
   case 1  % Normal marginal distribution
      
      mean = varargin{1};
      stdv = varargin{2};
      
      P = 0.5+erf(((x-mean)/stdv)/sqrt(2))/2;
      
   case 2  % Lognormal marginal distribution
      
      lambda = varargin{1};
      zeta   = varargin{2};
      
      z = ( log(x) - lambda ) / zeta;
      P = 0.5+erf(z/sqrt(2))/2;
      
   case 3  % Gamma distribution
      
      lambda = varargin{1};
      k      = varargin{2};
      
      P = gammainc(lambda*x,k);
      
   case 4  % Shifted exponential distribution
      
      lambda = varargin{1};
      x_zero = varargin{2};
      
      P = 1 - exp( -lambda*(x-x_zero) );
      
   case 5  % Shifted Rayleigh distribution
      
      a      = varargin{1};
      x_zero = varargin{2};
      
      P = 1 - exp( -0.5*( (x-x_zero)/a ).^2 );
      
   case 6  % Uniform distribution
      
      a = varargin{1};
      b = varargin{2};
      
      P = (x-a) / (b-a);
      
   case 7  % Beta distribution
      
      q = varargin{1};
      r = varargin{2};
      a = varargin{3};
      b = varargin{4};
      
      x01 = (x-a) / (b-a);
      P = betainc(x01,q,r);
      
   case 8  % Chi-square distribution
      
      nu     = varargin{1};
      lambda = 0.5;
      k      = nu/2;
      
      P = gammainc(lambda*x,k);
      
   case 11  % Type I largest value distribution ( same as Gumbel distribution )
      
      u_n = varargin{1};
      a_n = varargin{2};
      
      P = exp( -exp( -a_n*(x-u_n) ) );
      
   case 12  % Type I smallest value distribution
      
      u_1 = varargin{1};
      a_1 = varargin{2};
      
      P = 1 - exp( -exp( a_1*(x-u_1) ) );
      
   case 13  % Type II largest value distribution
      
      u_n = varargin{1};
      k   = varargin{2};
      
      P = exp( -(u_n./x).^k );
      
   case 14  % Type III smallest value distribution
      
      u_1     = varargin{1};
      k       = varargin{2};
      epsilon = varargin{3};
      
      P = 1 - exp( -( (x-epsilon) / (u_1-epsilon) ).^k );
      
   case 15  % Gumbel distribution ( same as type I largest value distribution )
      
      u_n = varargin{1};
      a_n = varargin{2};
      
      P = exp( -exp( -a_n*(x-u_n) ) );
      
   case 16  % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
      
      u_1 = varargin{1};
      k   = varargin{2};
      
      P = 1 - exp( -(x/u_1).^k );
      
   case 18  % (Reserved for Laplace marginal distribution)
      
   case 19  % (Reserved for Pareto marginal distribution)
      
   case 51  % Truncated normal marginal distribution
      
      mean = varargin{1};
      stdv = varargin{2};
      xmin = varargin{3};
      xmax = varargin{4};
      
      P = ( erf(((x-mean)/stdv)/sqrt(2))/2 - erf(((xmin-mean)/stdv)/sqrt(2))/2 ) / ...
          ( erf(((xmax-mean)/stdv)/sqrt(2))/2 - erf(((xmin-mean)/stdv)/sqrt(2))/2 );
      Imin = find(x<xmin); P(Imin) = 0;
      Imax = find(x>xmax); P(Imax) = 1;
      
   otherwise
      
end
