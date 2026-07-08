function [znull]=computeZnull(A,b,T,bT,xUp,xDn, zmin,zmax)
[nT,mT]=size(T);
[nA,mA]=size(A);

Aineq=[eye(nT) -1*T zeros(nT,1); zeros(nA,nT) A zeros(nA,1)];
%keyboard;
bineq=[bT; b];
%keyboard;

f=[zeros(1,nT) zeros(1,mA) 1];
Aeq=[eye(nT) zeros(nT,mA) zmax'-zmin'];
beq=[zmax'];
lb=[ones(1,nT)*-1*inf   xDn      0];
ub=[ones(1,nT)*inf      xUp   1];

xnull=cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);

lambda=xnull(nT+mA+1,1);
znull=lambda*(zmin-zmax)'+zmax';
%keyboard;
end