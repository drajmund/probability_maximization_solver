function [fval, fval_lower, grad, grad_lower, grad_upper, p, gradF, error] = computeFunctionGradAll(z, v, Sigma, error_level, m)
%computeFunctionGradAll Computation of the gradients (considering error) of the function v'z+log(F(z))=v'z-fi  at the point z, fi=-log(F(z))
%
%[fval, fval_lower, grad, grad_lower, grad_upper, p, gradF, error] = computeFunctionGradAll(z, v, Sigma, error_level)
%PARAMETERS:
%IN
%   z           - point z
%   v           - linear component vector
%   Sigma       - covariance matrix
%   error_level - error leve for the computation of the MVN CDF
%OUT
%   fval        - function value
%   fval_lower  - lower estimation of fval
%   grad        - gradient vector
%   grad_lower  - lower estimation of grad
%   grad_upper  - higher estimation of grad
%   p           - 
%   gradF       - 
%   error       - error level

plane_point = z*v';

%lower integration limit column vector
a = -inf*ones(1,length(v))'; 
%uses a randomized quasi-random rule with m points
m=getParameter('NumberOfEstimationPointsFUNCTION');
%function to estimate MVN probability, created by Alan Genz
[ p error ] = mvnmethod( m, Sigma, a, z ); %disp([ p e ]);
% options=statset('MaxFunEvals',1e14,'TolFun', min(error_level));
% [ p error ] = mvncdf(z,zeros(1,length(z)),Sigma,options);

p_lower=p-error;
p_upper=p+error;

%returns the p multipied with -1 because we need the maximum but instead we 
%use minimization function (e.g. fminsearch)
phi       = log(p);
phi_lower = log(p_lower);
phi_upper = log(p_upper);

fval       = -1*(plane_point + phi);
fval_lower = -1*(plane_point + phi_lower);
fval_upper = -1*(plane_point + phi_upper);

m=getParameter('NumberOfEstimationPointsGRADIENT');
%m=m*10;
[ gradF error_gradF ] = gradientMVN(Sigma, z, m, error_level);
gradF_lower = gradF - error_gradF;
gradF_upper = gradF + error_gradF;

parc_phi       = (gradF/(p));
parc_phi_lower = (gradF_lower/(p_upper));
parc_phi_upper = (gradF_upper/(p_lower));

grad =       -1*(v + parc_phi);
grad_upper = -1*(v + parc_phi_lower);   %!!
grad_lower = -1*(v + parc_phi_upper);   %!!

%check
for i=1:length(v)
    if isnan(grad(i)) || isinf(grad(i))
        grad(i)=-1*v(i);
    end
    if isnan(grad_lower(i)) || isinf(grad_lower(i))
        grad_lower(i)=-1*v(i);
    end
    if isnan(grad_upper(i)) || isinf(grad_upper(i))
        grad_upper(i)=-1*v(i);
    end
end