function  x = u_to_x(u,probdata)
%U_TO_X    Transformation between u and x space
%
%   x = u_to_x(u,marg,varargin)
%
%   Function 'u_to_x' perform the transformation between standard
%   normal space and original space. 
%
%   Output: - x           = values of r.v.'s in original space
%   Input:  - u           = values of r.v.'s in standard normal space
%           - probdata

transf_type = probdata.transf_type;
marg = probdata.marg;
nx = size(u,2);
nrv = size(marg,1);

x = zeros(nrv,nx);

switch transf_type
   
   case 1
      
   case 2
      
      %for i = 1 : nrv
      %   switch marg(i,1)
      %      case 1  % Normal distribution
      %         x(i,:) = u(i,:) * marg(i,3) + marg(i,2);
      %      case 2  % Lognormal distribution
      %         zeta = ( log( ( 1 + (marg(i,3)/marg(i,2))^2 ) ) )^0.5;
      %         lambda = log( marg(i,2) ) - 0.5 * zeta^2;
      %         x(i,:) =  exp(  u(i,:) * zeta + lambda  ) ;
      %      case 4  % Shifted exponential marginal distribution
      %         b = marg(i,3);
      %         a = marg(i,2) - b;
      %         x(i,:) = a - b * log(1-ferum_cdf(1,u(i,:),0,1));
      %      case 5  % Shifted Rayleigh marginal distribution
      %         b = marg(i,3) / (1-pi/4)^0.5;
      %         a = marg(i,2) - b*pi^0.5/2;
      %         x(i,:) = a + b * ( -log( 1-ferum_cdf(1,u(i,:),0,1) ) ).^0.5;
      %      case 6  % Uniform marginal distribution
      %         a = marg(i,2) - sqrt(3)*marg(i,3);
      %         b = marg(i,2) + sqrt(3)*marg(i,3);
      %         x(i,:) = a + (b-a) * ferum_cdf(1,u(i,:),0,1);
      %      case 11 % Type I Largest Value or Gumbel marginal distribution
      %         b = 6^0.5/pi*marg(i,3);
      %         a = marg(i,2) - 0.5772156649*b;
      %         x(i,:) = a - b * log( -log( ferum_cdf(1,u(i,:),0,1) ) );
      %      otherwise
      %   end
      %end
      
   case 3
      
      Lo = probdata.Lo;
      z = Lo * u;
      for i = 1 : nrv
         switch marg(i,1)
            case 1  % Normal distribution
               x(i,:) = z(i,:) * marg(i,3) + marg(i,2);
            case 2  % Lognormal distribution
               lambda = marg(i,5);
               zeta = marg(i,6);
               x(i,:) = exp( z(i,:) * zeta + lambda  );
            case 3  % Gamma distribution
               lambda = marg(i,5);
               k = marg(i,6);
               mean = marg(i,2);
               for j = 1 : nx
                  normal_val = ferum_cdf(1,z(i,j),0,1);
                  %x(i,j) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val); % Doesn't work
                  x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
               end
            case 4  % Shifted exponential distribution
               lambda = marg(i,5);
               x_zero = marg(i,6);
               x(i,:) = x_zero + 1/lambda * log( 1 ./ ( 1 - ferum_cdf(1,z(i,:),0,1) ) );
            case 5  % Shifted Rayleigh distribution
               a = marg(i,5);
               x_zero = marg(i,6);
               x(i,:)= x_zero + a * ( 2*log( 1 ./ (1-ferum_cdf(1,z(i,:),0,1)) ) ) .^0.5;
            case 6  % Uniform distribution
               a = marg(i,5);
               b = marg(i,6);
               x(i,:) = a + (b-a) * ferum_cdf(1,z(i,:),0,1);
            case 7  % Beta distribution
               q = marg(i,5);
               r = marg(i,6);
               a = marg(i,7);
               b = marg(i,8);
               mean = marg(i,2);
               for j = 1 : nx
                  normal_val = ferum_cdf(1,z(i,j),0,1);
                  %x01 = fzero('zero_beta',(mean-a)/(b-a),optimset('fzero'),q,r,normal_val); % Doesn't work
                  x01 = fminbnd('zero_beta',0,1,optimset('fminbnd'),q,r,normal_val);
                  % Transform x01 from [0,1] to [a,b] interval
                  x(i,j) = a + x01 * ( b - a );
               end
            case 8 % Chi-square distribution
               lambda = 0.5;
               nu = marg(i,5);
               k = nu/2 ;
               mean = marg(i,2);
               for j = 1 : nx
                  normal_val = ferum_cdf(1,z(i,j),0,1);
                  %x(i,j) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val); % Doesn't work
                  x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
               end
            case 11 % Type I largest value distribution ( same as Gumbel distribution )
               u_n = marg(i,5);
               a_n = marg(i,6);
               x(i,:) = u_n - (1/a_n) * log( log( 1 ./ ferum_cdf(1,z(i,:),0,1) ) ); 
            case 12 % Type I smallest value distribution
               u_1 = marg(i,5);
               a_1 = marg(i,6);
               x(i,:) = u_1 + (1/a_1) * log( log( 1 ./ ( 1 - ferum_cdf(1,z(i,:),0,1) ) ) );
            case 13 % Type II largest value distribution
               u_n = marg(i,5);
               k = marg(i,6);
               x(i,:) = u_n * log( 1 ./ ferum_cdf(1,z(i,:),0,1) ) .^ (-1/k);
            case 14 % Type III smallest value distribution
               u_1 = marg(i,5);
               k = marg(i,6);
               epsilon = marg(i,7);
               x(i,:) = epsilon + ( u_1 - epsilon ) * log( 1 ./ ( 1 - ferum_cdf(1,z(i,:),0,1) ) ) .^(1/k);
            case 15 % Gumbel distribution ( same as type I largest value distribution )
               u_n = marg(i,5);
               a_n = marg(i,6);
               x(i,:) = u_n - (1/a_n) * log( log( 1 ./ ferum_cdf(1,z(i,:),0,1) ) ); 
            case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
               u_1 = marg(i,5);
               k = marg(i,6);
               x(i,:) = u_1 * log( 1 ./ ( 1 - ferum_cdf(1,z(i,:),0,1) ) ) .^(1/k);
            case 18 % (Reserved for Laplace distribution)
            case 19 % (Reserved for Pareto distribution)
            case 51 % Truncated normal marginal distribution
               mean = marg(i,5);
               stdv = marg(i,6);
               xmin = marg(i,7);
               xmax = marg(i,8);
               x(i,:) = mean + stdv * inv_norm_cdf( ...
                           ferum_cdf(1,(xmin-mean)/stdv,0,1) + ...
                           (ferum_cdf(1,(xmax-mean)/stdv,0,1)-ferum_cdf(1,(xmin-mean)/stdv,0,1)) * ferum_cdf(1,z(i,:),0,1) ...
                                                  );
               Imin = find(x(i,:)<xmin); x(i,Imin) = xmin;
               Imax = find(x(i,:)>xmax); x(i,Imax) = xmax;
            otherwise
         end
      end
      
   case 4
      
      %opt = optimset([]);
      %
      %str = '';
      %for i = 1:nrv
      %   eval(['x(' num2str(i) ') = fzero(''zero_u_to_x'',marg(' num2str(i) ',2),opt,u,marg,' num2str(i) str ');']);
      %   str = [str ',x(' num2str(i) ')'];
      %end
      %
      %% x(1) = fzero('zero_u_to_x',marg(1,2),opt,u,marg,1);
      %% x(2) = fzero('zero_u_to_x',marg(2,2),opt,u,marg,2,x(1));
      %% ...      opt = optimset([]);

   otherwise
      
end
