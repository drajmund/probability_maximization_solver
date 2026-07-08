function [ derivative, error ] = derivativeMVN(r, point, numberOfPoints, DIM, abs_err)
%derivativeMVN Finding an element of the gradient of a multidimensional probability distribution function
%F(z1,...,zn)
%Implemented according to Prekopa's formula from his book of "Stochastic Programming", 
%pages 203-204., section 6.6.4 Gradients
%
%[ derivative, error ] = derivativeMVN(r, point, numberOfPoints, DIM, abs_err)
%PARAMETERS:
%IN
%   r               - covariance matrix, 
%   point           - point, where we request derivative,
%   numberOfPoints  - (not used, for old Genz code),
%   DIM             - dimension where we calculate derivative, 
%   abs_err         - mvncdf error level
%OUT
%   derivative      - derivative value
%   error           - error level returned by mvncdf
%
%   See also gradientMVN.

if size(r,1)~=size(r,2)
    error('Input parameter *r* is not a square matrix')
end
if size(r,1)~=length(point)
    error('Input parameter *point* does not have a proper size')
end

n=size(r,1); %dimenstion
%n - z1, z2, ..., zn
%Fi(z1,z2,...,zn|zi)
for i=DIM
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
      
      [ p error ] = qsimvnv( numberOfPoints, r_new, a, b ); %n-1 dimension cdf
      %options=statset('MaxFunEvals',1e12,'TolFun',abs_err);
      %[p error] = mvncdf(b,zeros(1,length(b))',r_new,options);

      derivative= p*normpdf(point(i));
end

