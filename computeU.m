% nT=T sorainak száma
function [u,theta,x,p,fval]=computeU(f,Aineq,bineq,Aeq,beq,lb,ub,nT,x0)

    [x,fval,exitflag,output,lambda]= cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub, x0);
    u=lambda.ineqlin(1:nT,1)*(-1);
    theta=lambda.eqlin*(-1);
    p=exp(-fval);
    
end


