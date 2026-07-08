%% CSE változtatások (2021.03.20.) 
% computeFunctionGradFi.m

%% DR változások (2017.08.23.)
%új függvény: getSigmaSymmetricPacked, mvnmethod
%módosított függvények:
%getParameter - MethodOfCalculationofMVNCDF
%computeFi - mvnmethod hívása qsimvnv helyett
%computeBrick - mvnmethod hívása qsimvnv helyett
%computeFunctionGradAll - mvnmethod hívása qsimvnv helyett
%computeFunctionValueAll - mvnmethod hívása qsimvnv helyett

%új paraméter a num, ami a futtatás sorszámát jelzi,
%a callCoffeProblem-bõl ciklusban hívom a CoffeeProblem-et és
%a futtatás sorszáma a képzõdõ xls fájlnév elejére kerül 1_...

%a gradiens pontosság beállításához már nem a getParameter-bõl szedi a
%NumberOfEstimationPointsLOW és NumberOfEstimationPointsHIGH értékeket,
%hanem az m=... értéket itt lehet beállítani a CoffeProblem-ben, amit
%tovább is adunk a computeMaximum-nak

%kezdeti irány visszaadása és kiírása Excelbe
%grad_z_starter=grad_z_actual
%grad_z_starters(:,i)=grad_z_starter;

%error_levels
%új paraméter x0_all, DR 20250607
function [x,Number_of_Genz_calls, Number_of_inner_iterations,pvec,grad_norm_objective,diff_sum_fi,timeCMax,x0_all]=DemandSideProblem(probLevel,num, m_m0,ptol,fval,A,b,T,bT,xUp,xDn, Sigma, zmin, zmax, R_num)
%function [x,Number_of_Genz_calls, Number_of_inner_iterations,pvec,grad_norm_objective,diff_sum_fi]=DemandSideProblem(probLevel,num, m_m0, ptol, Sigma, testNumber, zmin, zmax)
 global R_type;
 global std_type;

 global maxItNumber; 
 global diagZ;
 global timeGenz;
 global WrightDim;
 timeGenz=zeros(1,maxItNumber);
 global kintVan;
 kintVan=false;
 global use_batteries;
 
    tic
    dimension=length(Sigma);
    error_levels=ones(1,dimension)*10^-4;
    
    %%SZTNORM DLL
    if getParameter('MethodOfCalculationofMVNCDF') == 2 && not(libisloaded('sztnorm'))
        loadlibrary('sztnorm')
    end
    %libfunctions('sztnorm')

    %% Állítható paraméterek
    %ptol=0.99
    %testNumber=10; %Wim simple %testNumber=11; %Wim cash_matching %testNumber=12; %Wim Isr48
    %Probability level of the optimal solution
    %probLevel=0.99;

    tolStop=0.01;  

    %%  A mátrixok felépítése
    [nT,mT]=size(T);
    [nA,mA]=size(A);

    w=zeros(nT,nT+2);
    w(:,1)=zmax';
    w(:,end)=zeros(nT,1);
    fiw=zeros(1,nT+2);
    for i=1:nT
        for j=2:(nT+1)
            if i==(j-1)
                w(i,j)=0;
            else
                w(i,j)=zmax(j-1);
            end
        end
    end
    for i=1:(nT+2)
        [fiw(i)]=computeFi(w(:,i)',Sigma, 'not_zero');
    end
    tcompZnull_start=tic;
    [znull]=computeZnull(A,b,T,bT,xUp,xDn, zmin,zmax);
    tcompZnull=toc(tcompZnull_start);
    [finull]=computeFi(znull',Sigma, 'zero');
    %keyboard
    %finull = 0.6
    f=[finull];
    % zVectors=[(znull-E)./std];
    zVectors=[znull];
    for i=(nT+2):-1:1
        f=[f fiw(i)];
        % zVectors=[zVectors (w(i)-E)./std];
        zVectors=[zVectors w(:,i)];
    end

    f=[f fval];

    Aineq=[zVectors (-1)*T;
        zeros(nA,nT+3) A];
    bineq=[bT; b];
    Aeq=[ones(1,nT+3) zeros(1,mT)];

    beq=1;

    % x pozitív-e,lambda pozitív-e
    [nf,mf]=size(f);
    % y-okra lehet negatív is!!!!
    if use_batteries
        lb=[zeros(1,mf-2*nT) -100000*ones(1,nT) zeros(1,nT)];
    else
        lb=[zeros(1,mf)];
    end
    ub=[];
    u_z_Fi_p_fval=[];

    %Used in ExcelExport
    f0=f; Aineq0=Aineq; bineq0=bineq; Aeq0=Aeq;
    beq0=beq; lb0=lb; ub0=lb0;

    %% Számolás
    itNumber=0; %iterációk száma
   
    diff_sum_fi=zeros(1,maxItNumber);
    p_error=zeros(1,maxItNumber);
    grad_error=zeros(dimension,maxItNumber);
    zteto=zeros(dimension,maxItNumber);
    time_iter=zeros(1,maxItNumber);
    m_m=zeros(1,maxItNumber);
    errorEstimation= zeros(1,maxItNumber);
    Bvector= zeros(1,maxItNumber);
   % m_m(i)=m_m0;
    m_actual=m_m0;
    %DR
    x0 = [];
    x0_all = [];
    %x0 = zeros(1,53);
    for i=1:maxItNumber
        itStart=tic;
        [u,theta,x,p,fval]=computeU(f,Aineq,bineq,Aeq,beq,lb,ub,nT,x0);
        
        %DR??? x0=x;
        [nzVectors,mzVectors]=size(zVectors);
        lambda=x(1:(mzVectors));
        z0=(zVectors*lambda);%zteto
        zteto(:,i)=z0';
         
        z0T=z0';
        if  z0T<zmin 
            kintVan=true;
        end
         if  z0T>zmax
            kintVan=true;
        end
        fi0=computeFi(z0',Sigma,'not_zero');%fi_zteto
        fizi=f(1:(mzVectors));
        
        %diff_sum_fi(i)=fval-fi0;
        diff_sum_fi(i)=fizi*lambda-fi0;
          
        [fval2, grad, p2, gradFi, e,error] = computeFunctionGradFi(z0', u', Sigma, m_actual);
   
        grad_norm_objective(i)=norm(u'-gradFi);
        grad_norm(i)=norm(gradFi);
        u_norm(i)=norm(u);
        p_error(i)=e;
        
         grad_error(:,i)=error';
            if(m_m(i)<m_actual)
                m_m(i)=m_actual;
            else
                m_actual=m_m(i);
            end
        
        startCMax=tic;
        itNumber
        [z, f_z_actual, e Number_of_Genz_calls(i), Number_of_inner_iterations(i), grad_z_starter]=...
            searchMINGradient(u', zmin, zmax, Sigma,tolStop,fi0,z0',theta,itNumber,m_m(i));
        Fi=computeFi(z, Sigma, 'model');

        timeCMax(i)=toc(startCMax);
       
        if  z<zmin
            kintVan=true;
        end
        if  z>zmax
            kintVan=true;
        end
        
        grad_z_starters(:,i)=grad_z_starter;
        zVectors=[z' zVectors];

        s=[u;z';Fi;0;p;fval];
        p_prev=p;
        pvec(i)=p;
        u_z_Fi_p_fval=[u_z_Fi_p_fval s];

        f=[Fi f];
        itNumber=itNumber+1;
        s=[z'; zeros(nA,1)];
        %baj lehet a sebességgel, a mátrixok mérete mindig változik
        Aineq=[s Aineq];
        Aeq=[1 Aeq];

        % x pozitív-e,lambda pozitív-e
        %feltételek
        lb=[0 lb];
        [minStop]=computeErrorEstimation(zmin,zmax,grad_z_starter,z0');
        errorEstimation(i)=minStop;   
        Bvector(i)=(fval-fi0)+minStop+diagZ*sqrt(sum(error.*error));%(57) képlet
        time_iter(i)=toc(itStart);

        %DR 20250607 x0 gyûjtése, hogy megnézzük, mi a gond a 96 dim
        [nx,mx]=size(x);
        x0=x((nx-mA+1):nx,1);
        x0_all(i,:)=x0;
    end
    %keyboard;
    %% Adatok és eredmények excelbe
    
    time=toc;
    time_tcompZnull_diagZ=[time tcompZnull diagZ];
    TimeIter_TimeGenz=[time_iter;timeGenz];
    %ExcelExport(x,u_z_Fi_p_fval, f0, Aineq0, bineq0, Aeq0, beq0, lb0, ub0, itNumber,probLevel,ptol,Number_of_Genz_calls, Number_of_inner_iterations, time_tcompZnull_diagZ, num, grad_z_starters, m_m, grad_norm_objective, grad_norm, u_norm,diff_sum_fi,p_error, grad_error,time_iter,zteto,errorEstimation,TimeIter_TimeGenz,Bvector,R_num);
    [nx,mx]=size(x);
    x
    x0=x((nx-mA+1):nx,1)
    T*x0
    itNumber

    time_tcompZnull_diagZ
    kintVan
    timeCMax
    
    if use_batteries
        [fval,A,b,T,bT,xUp,xDn, Sigma, EX_, std_]=...
            DataloadDemandSideProblem3D_y_all(dimension, R_type, R_num, std_type, 1, x0, x0_all, pvec);
    else
        [fval,A,b,T,bT,xUp,xDn, Sigma, EX_, std_]=...
        DataloadDemandSideProblem3D_cont(dimension, R_type, R_num, std_type, 1, x0, x0_all, pvec);
    end
end