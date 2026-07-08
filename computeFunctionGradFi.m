%computation of v'z+log(F(z)) at the point z
function [fval, grad, p, gradFi, e,error] = computeFunctionGradFi(z, v, Sigma,m)
%function [fval, grad, p, gradFi, e,error] = computeFunctionGradFi(z, v, Sigma, testNumber,m)
%computeFunctionGrad Computation of the gradient of the function v'z+log(F(z))=v'z-fi  at the point z, fi=-log(F(z))
%
%[fval, grad, p, gradF, e] = computeFunctionGrad(z, v, Sigma, testNumber)
%
%Kivettem. Edit
%[Sigma, EX, std]=getSigmaAndExpectedValue(testNumber);

plane_point = z*v';

%lower integration limit column vector
a = -inf*ones(1,length(z))'; 
%uses a randomized quasi-random rule with m points
m=getParameter('NumberOfEstimationPointsPRECISE');
%function to estimate MVN probability
[ p e ] = mvnmethod( m, Sigma, a, z );

%returns the p multipied with -1 because we need the maximum but instead we 
%use minimization function (e.g. fminsearch)
fval = -1*(plane_point+log(p));

[ gradF error] = gradientMVN(Sigma, z, m );
gradFi=-1*(gradF/(p));
grad=-1*(v + (gradF/(p)));
for i=1:length(v)
    if isnan(grad(i)) || isinf(grad(i))
        grad(i)=-1*v(i);
    end
end