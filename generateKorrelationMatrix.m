function [R] = generateKorrelationMatrix(dimension_, R_type, R_value)

R_value = 1.0/dimension_;

if strcmp(R_type, 'decr')
    R = zeros(dimension_);
    start_column = 2;
    for row=1:dimension_
        value = 1.0-R_value;
        for column=start_column:dimension_
            R(row,column) = value;
            value = value - R_value;
        end
        start_column = start_column + 1;
    end
elseif strcmp(R_type, 'big')
    R = zeros(dimension_);
    start_column = 2;
    for row=1:dimension_
        value = 0.95;
        for column=start_column:dimension_
            R(row,column) = value;
            value = value - 0.0;
        end
        start_column = start_column + 1;
    end
elseif strcmp(R_type, 'fix')
    R = zeros(dimension_);
    start_column = 2;
    for row=1:dimension_
        value = R_value;
        for column=start_column:dimension_
            R(row,column) = value;
            value = value - 0.0;
        end
        start_column = start_column + 1;
    end
elseif strcmp(R_type,'no_corr')
    R = zeros(dimension_)
end

%make symmetric and 1-s on main 
R = R + R' + eye(dimension_);