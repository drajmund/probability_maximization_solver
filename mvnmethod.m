function [ p e ] = mvnmethod( m, Sigma, a, z)

method=getParameter('MethodOfCalculationofMVNCDF');

if method==1
    [ p e ] = qsimvnv( m, Sigma, a, z );
    %'qsimvnv'
else
    seed=floor(rand()*10000);    
    [ p ] = calllib('sztnorm','sztnorcall',m,length(Sigma),z',getSigmaSymmetricPacked(Sigma),seed);
    e=0;
    %'sztnorm'
end