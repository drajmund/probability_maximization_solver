function [ grad e ] = gradientMVN(r, point, numberOfPoints, error_level )
%gradientMVN Finding an element of the gradient of a multidimensional probability distribution function
%
%Implemented according to Prekopa's formula from his book of "Stochastic Programming", 
%pages 203-204., section 6.6.4 Gradients
%
%[ grad error ] = gradientMVN(r, point, numberOfPoints, error_level )
%PARAMETERS:
%IN
%   r               - covariance matrix, 
%   point           - point, where we request derivative,
%   numberOfPoints  - (not used, for old Genz code),
%   DIM             - dimension where we calculate derivative, 
%   abs_err         - mvncdf error level
%OUT
%   grad            - gradient vector
%   error           - error level returned by mvncdf
%
%   See also derivativeMVN.

% rand_vec = normrnd(0,1,1,length(r));
% grad = 0.045*rand_vec./sqrt(rand_vec*rand_vec');
% e = 0;
% return

if size(r,1)~=size(r,2)
    error('Input parameter *r* is not a square matrix')
end
if size(r,1)~=length(point)
    error('Input parameter *point* does not have a proper size')
end

n=size(r,1); %dimenstion
%n - z1, z2, ..., zn
%Fi(z1,z2,...,zn|zi)
parfor i=1:n
    grad_(i)=derivativeMVN(r, point, numberOfPoints, i, 0); 
end
grad = grad_;
e = 0;      
end

% fprintf('[');
% fprintf('%g, ', grad(1:end-1));
% fprintf('%g]\n', grad(end));


