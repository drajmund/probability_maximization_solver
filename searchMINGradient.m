function [z_actual, f_z_actual, e, Number_of_Genz_calls, Number_of_inner_iterations, grad_z_starter] = searchMINGradient(v_vector, z_min, z_max, Sigma, tol_stop, fi0, z0, theta, iteration_k,m) 
%searchMINGradient Calculates minimum of the function v'z+log(F(z))=v'z-fi, fi=-log(F(z))
%
%[z_actual, fval, e, fineSearch] = searchMINGradient(v_vector, z_min, z_max, Sigma, tol_stop, fi0, z0, theta, iteration_k)
%PARAMETERS:
%IN
%   v_vector    - linear component vector
%   z_min       - "smallest" corner point of the brick
%   z_max       - "largest" corner point of the brick
%   Sigma       - covariance matrix
%   tol_stop    - stopping tilerance of the LP problem solver 
%   fi0         - starting value of Phi
%   z0          - starting value of z (start point of the search)
%   theta       - theta value of dual solution
%   iteration_k - iteration number of the LP problem solver 
%
%OUT
%   f_z_actual  - the value of the objective function: max v'z+log(F(z))
%               in the solution point z
%   z_actual    - the solution z vector 
%   e           - error level of the solution

dimension=size(z0,2);
z_actual=z0; %starting point
error_levels=ones(1,dimension)*10^-4;
error_level=10^-4;

%%Calculate grad_z_actual
[f_z_actual, f_z_actual_lower, grad_z_actual, grad_z_actual_lower, grad_z_actual_upper, p_z_actual, gradF_z_actual, e] = ...
    computeFunctionGradAll(z_actual, v_vector, Sigma, error_levels,m);
%-------------------------------------------------------------------------
%GRADIENT SEARCH
%-------------------------------------------------------------------------

%multiplier of Brick
BrickSizeMultiplier = 1;
   
z_counter_previous = z_actual;
f_counter_previous = f_z_actual;

%Gradient search with the modified brick
[z_actual_new, better, Number_of_Genz_calls, Number_of_inner_iterations, grad_z_starter] = ...
    searchMINGoldenSection(v_vector, Sigma, z_actual, BrickSizeMultiplier*z_min, BrickSizeMultiplier*z_max, iteration_k,m);

if better %f_z_actual_new<z_actual
    z_actual = z_actual_new; 
    %we need to recalculate p_z_actual
    [f_z_actual, f_z_actual_lower, grad_z_actual, grad_z_actual_lower, grad_z_actual_upper, p_z_actual, gradF_z_actual, e] = ...
       computeFunctionGradAll(z_actual, v_vector, Sigma, error_levels,m);
end                  
         
