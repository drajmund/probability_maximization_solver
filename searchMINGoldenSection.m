function [z_actual, better, Number_of_Genz_calls, Number_of_inner_iterations, grad_z_starter, z_counter, zzz, fff] = searchMINGoldenSection(v, Sigma, z0, z_min, z_max, iteration_k, m) 
%searchMINGoldenSection Calculates maximum of the function v'z+log(F(z))=v'z-fi, fi=-log(F(z)) 
%Gradient search only in the steepest descent gradient direction with Golden
%section seteps
%
%[z_actual, better] = searchMINGoldenSection(v, Sigma, z0, z_min, z_max, testNumber) 

Number_of_Genz_calls=0;
Number_of_inner_iterations=0;

number_of_LEFT_side_skip=0;
number_of_RIGHT_side_skip=0;
LEFT_side_skip=0;
RIGHT_side_skip=0;

allowed_number_of_LEFT_side_skip=getParameter('AllowedNumberOfLeftSideSkip');

dimension=size(z0,2);
error_levels=ones(1,dimension)*10^-4;
error_level=10^-4;
z_actual=z0; %starting point
gr = (sqrt(5) - 1) / 2;
better=1;
%calculate function value (f_z_actual), gradient (grad_z_actual) for z_actual
[f_z_actual, f_z_actual_lower, grad_z_actual, grad_z_actual_lower, grad_z_actual_upper, p_z_actual, gradF_z_actual, e] = ...
    computeFunctionGradAll(z_actual, v, Sigma, error_levels,m);
%p_z_actual
grad_z_starter = grad_z_actual;
Number_of_Genz_calls=Number_of_Genz_calls+1+dimension;
max_DIM_value=max(-1*(grad_z_actual));
DIM = find(-1*(grad_z_actual)==max_DIM_value,1,'first');

%calculate counter point
[z_counter, diagonal_z_min_z_max] = ...
    calculateZCounterSimple(dimension, z_min, z_max, z_actual, grad_z_actual);
%tol_stop=diagonal_z_min_z_max/10^9;%10^-1;

if max(abs(z_counter-z_actual))>10^-6    
    % a - c - d - b
    %calculate inner points with Golden section
    [a, b, c, d] = calculateInnerPoints(dimension, z_actual, z_counter);
    
    %calculate function values
    [f_c, f_c_lower, p_c, error_c] = computeFunctionValueAll(c, v, Sigma, error_level,m);
    Number_of_Genz_calls=Number_of_Genz_calls+1;
    [f_d, f_d_lower, p_d, error_d] = computeFunctionValueAll(d, v, Sigma, error_level,m);
    Number_of_Genz_calls=Number_of_Genz_calls+1;
    
    while (1)
        Number_of_inner_iterations=Number_of_inner_iterations+1;
        zzz(Number_of_inner_iterations, :)=[a, b, c, d];
        fff(Number_of_inner_iterations, :)=[f_c, f_d];
        %calculate funcion values        
        if LEFT_side_skip==0
            [f_c, f_c_lower, p_c, error_c] = computeFunctionValueAll(c, v, Sigma, error_level,m);
            Number_of_Genz_calls=Number_of_Genz_calls+1;
        end
        if RIGHT_side_skip==0
            [f_d, f_d_lower, p_d, error_d] = computeFunctionValueAll(d, v, Sigma, error_level,m);
            Number_of_Genz_calls=Number_of_Genz_calls+1;
        end
     
            if f_c < f_d
                    %'skip right side of the region'
                    number_of_RIGHT_side_skip=number_of_RIGHT_side_skip+1;                  
                    RIGHT_side_skip=1;
                    LEFT_side_skip=0;
                    
                    b = d;
                    f_b = f_d;
                    d = c;
                    f_d = f_c;
                    c = b - gr * (b - a);                  
            else
                    %'skip left side of the region'
                    number_of_LEFT_side_skip=number_of_LEFT_side_skip+1;
                    LEFT_side_skip=1;
                    RIGHT_side_skip=0;
                    if number_of_LEFT_side_skip > allowed_number_of_LEFT_side_skip
                       %'too many number_of_LEFT_side_skip'
                        break
                    end
                    
                    a = c;
                    f_a = f_c;
                    c = d;
                    f_c = f_d;
                    d = a + gr * (b - a);                 
            end         
          
    end

                z_actual=c;
                f_c;
                if f_z_actual<f_c
                    f_z_actual;
                    better=0;
                end
   
end
end

%calculatin inner points according to golden section 
% a - c - d - b
function [a, b, c, d] = calculateInnerPoints(dimension, z_actual, z_counter)
gr = (sqrt(5) - 1) / 2;
    for i=1:dimension
        a(i) = z_actual(i);
        b(i) = z_counter(i);
        c(i) = b(i) - gr * (b(i) - a(i));
        d(i) = a(i) + gr * (b(i) - a(i));    
    end
end