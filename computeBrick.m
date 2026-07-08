function [z_min, z_max] = computeBrick(Sigma, p_goal)
%computeBrick Computes the size of the brick with probability p_goal
%
%[z_min z_max] = computeBrick(Sigma, p_goal)
%
format long

lambda0=0.1;
j=1;
p=0;
%uses a randomized quasi-random rule with m points
m=getParameter('NumberOfEstimationPointsBRICK');

while p<p_goal
    z_min=-1*(1+j*lambda0)*ones(1,length(Sigma));
    z_max=1*(1+j*lambda0)*ones(1,length(Sigma));
    %function to estimate MVN probability
    [ p error ] = mvnmethod( m, Sigma, z_min, z_max); 
    j=j+1;
end

lambda=0.005;
i=(j-2)*lambda0/lambda;
p=0;
while p<p_goal
    z_min=-1*(1+i*lambda)*ones(1,length(Sigma));
    z_max=1*(1+i*lambda)*ones(1,length(Sigma));
    %function to estimate MVN probability
    [ p error ] = mvnmethod( m, Sigma, z_min, z_max); 
    i=i+1;
end
