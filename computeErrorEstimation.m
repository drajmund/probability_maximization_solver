function [minStop]=computeErrorEstimation(zmin,zmax,grad,z)
    

f=grad';
Aineq=[];
bineq=[];
Aeq=[];
beq=[];
lb=zmin;
ub=zmax;
x0=[];
%keyboard;
options = cplexoptimset('TolFun', 1e-8);

[x,fval,exitflag,output,lambda]= cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub, x0);
%keyboard;
minStop=fval-(grad)*z';






end