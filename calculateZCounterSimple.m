function [z_counter, diagonal_z_min_z_max] = calculateZCounterSimple(dimension, z_min, z_max, z_actual, grad_z_actual)
%calculateZCounter Claculating the counter point on the wall for the gradient search
%
%[z_counter, DIM] = calculateZCounter(dimension, z_min, z_max, z_actual, grad_z_actual)
%

%tolerance to determine whether wall is too close or not
%tol=10^-4;
    diagonal_z_min_z_max=2*z_max(1)*sqrt(length(z_min));
    %grad_z_actual
    z_counter=z_actual-grad_z_actual*diagonal_z_min_z_max;
    
    %kiválasztjuk a maximum eltérést a vektorban és
    %visszaskálázzuk a vektort ennek segítségével, hogy ez
    %a maximum érték a falon legyen
    skalazo_faktor = 1.0;
    max_z_parameter = max(abs(z_counter));
    for i=1:length(z_min)
        if((z_counter(i) < z_min(i)) && (abs(z_counter(i)) == max_z_parameter) ) 
            skalazo_faktor = z_counter(i)/z_min(i);
            
        end
        if((z_counter(i) > z_max(i)) && (z_counter(i) == max_z_parameter)) 
            skalazo_faktor = z_counter(i)/z_max(i);
        end
    end
    %disp('skálázás után:')
    z_counter = z_counter./skalazo_faktor;
%     for i=1:length(z_min)
%         if(z_counter(i) < z_min(i) ) 
%             z_counter(i) = z_min(i);
%         end
%         if(z_counter(i) > z_max(i) ) 
%             z_counter(i) = z_max(i);
%         end
%     end
end