function [ retval ] = getParameter(parameterName)
retval = 0;
if strcmp(parameterName,'NumberOfEstimationPointsLOW') 
    retval = 100;
elseif strcmp(parameterName,'NumberOfEstimationPointsHIGH') 
    retval = 1000;
elseif strcmp(parameterName,'NumberOfEstimationPointsPRECISE') 
    retval = 100;
elseif strcmp(parameterName,'NumberOfEstimationPointsBRICK') 
    retval = 1000;
elseif strcmp(parameterName,'NumberOfEstimationPointsFUNCTION0') 
    retval = 100;
elseif strcmp(parameterName,'NumberOfEstimationPointsFUNCTION') 
    retval = 200;
elseif strcmp(parameterName,'NumberOfEstimationPointsMODEL') 
    retval = 1000;
elseif strcmp(parameterName,'NumberOfEstimationPointsGRADIENT') 
    retval = 400;
elseif strcmp(parameterName,'AllowedNumberOfLeftSideSkip') 
    retval = 2;
elseif strcmp(parameterName,'MethodOfCalculationofMVNCDF') 
    retval = 1;  %qsimvnv
    %retval = 2; %sztnorm
end