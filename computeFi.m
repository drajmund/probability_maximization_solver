function [fi_ret] = computeFi(z, Sigma, type)
%computeFi Computing -log(p+e)
%
%[fi_ret] = computeFi(z, Sigma)
format long

a = -inf*ones(1,length(z))'; 
%uses a randomized quasi-random rule with m points
if strcmp(type,'zero')
    m=getParameter('NumberOfEstimationPointsFUNCTION0');
elseif strcmp(type,'not_zero')
    m=getParameter('NumberOfEstimationPointsFUNCTION');
elseif strcmp(type,'model')
    m=getParameter('NumberOfEstimationPointsMODEL');
end

%function to estimate MVN probability
[ p e ] = mvnmethod( m, Sigma, a, z );
fi_ret=-log(p);
