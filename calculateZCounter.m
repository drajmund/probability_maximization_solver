function [z_counter, DIM] = calculateZCounter(dimension, z_min, z_max, z_actual, grad_z_actual)
%calculateZCounter Claculating the counter point on the wall for the gradient search
%
%[z_counter, DIM] = calculateZCounter(dimension, z_min, z_max, z_actual, grad_z_actual)
%

%tolerance to determine whether wall is too close or not
%tol=10^-4;

z_counter=zeros(1,dimension);
grad_z_actual_vector=-grad_z_actual;

for i=1:dimension
    if grad_z_actual(i) > 0
%         if((z_actual(i)-z_min(i))<tol)
%             %wall is too close
%             grad_z_actual_vector(i)=0;
%             step_number(i)=inf;
%         else
            step_number(i)=(z_actual(i)-z_min(i))/grad_z_actual_vector(i);
%        end            
    else      
%         if((z_max(i)-z_actual(i))<tol)
%             %wall is too close
%             grad_z_actual_vector(i)=0;
%             step_number(i)=inf;
%         else
            step_number(i)=(z_max(i)-z_actual(i))/grad_z_actual_vector(i);
%        end
    end                
end     
    
step_min=min(abs(step_number));  
DIM = find(abs(step_number)==step_min,1,'first');
z_counter=z_actual+grad_z_actual_vector*step_min;

if step_min==inf
    z_counter=z_actual;
end

end