function [ grad e ] = gradientMVN_simple(r, point, numberOfPoints, error_level )
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
for i=1:n
    h=1;
    for m=1:n
      if m==i 
          continue
      end
        point_new(h)=( point(m)-r(m,i)*point(i) ) / sqrt( 1 - (r(m,i))^2 );
        h=h+1;
    end

    %  Ri - (n-1)x(n-1) correlation matrix with sjk components
      jj=1;
      for j=1:n
          if j==i 
              continue
          end

          kk=1;
          for k=1:n
              if k==i 
                  continue
              end

              r_new(jj,kk)=( r(j,k)-r(j,i)*r(k,i) ) / ( sqrt( 1-( r(j,i) )^2 ) * sqrt( 1-( r(k,i) )^2 ));
              kk=kk+1;
          end
          jj=jj+1;
      end

      % ith element of the gradient is the a partial derivative: n-1 dimension cdf multiplied with pdf(point(i)) 
      a = -inf*ones(1,length(point))'; 
      b = point_new';
      
      [ p e(i) ] = qsimvnv( numberOfPoints, r_new, a, b ); %n-1 dimension cdf
%       if error_level(1)<10^-4
%         options=statset('MaxFunEvals',1e12,'TolFun',error_level(n));
%         [p e(i)] = mvncdf(b,zeros(1,length(b))',r_new,options);
%       else
%         [p e(i)] = mvncdf(b,zeros(1,length(b))',r_new);
%       end

      grad(i)= p*normpdf(point(i));
      
end

% fprintf('[');
% fprintf('%g, ', grad(1:end-1));
% fprintf('%g]\n', grad(end));


