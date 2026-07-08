function [fval, fval_lower, p, error] = computeFunctionValueAll(z, v, Sigma, error_level, m)
%computeFunctionValueAll Computation of function values (considering error) of v'z+log(F(z)) at the point z

plane_point = z*v';

%lower integration limit column vector
a = -inf*ones(1,length(v))'; 
%uses a randomized quasi-random rule with m points
m=getParameter('NumberOfEstimationPointsFUNCTION');
%function to estimate MVN probability, created by Alan Genz
[ p error ] = mvnmethod( m, Sigma, a, z );

% if error_level<10^-4
%     options=statset('MaxFunEvals',1e14,'TolFun',error_level);
%     [ p error ] = mvncdf(z,zeros(1,length(z)),Sigma, options);
% else
%     [ p error ] = mvncdf(z,zeros(1,length(z)),Sigma);
% end

p_lower=p-error;
p_upper=p+error;

%returns the p multipied with -1 because we need the maximum but instead we 
%use minimization function (e.g. fminsearch)
phi       = log(p);
phi_lower = log(p_lower);

fval       = -1*(plane_point + phi);
fval_lower = -1*(plane_point + phi_lower);


