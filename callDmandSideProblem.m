%%DemandSideProble - computeMaximum delted, it calls directly searchMINGradient 
%computeFi0 and computeFiModel deleted - computeFi last parameter true/false determines

%clear all
%rng(1);
% dimension = 24;
% global R_type;
% R_type = 'no_corr';
% global std_type;
% std_type = 'constant';
% 
% global save_directory;
% save_directory = '.\CEJOR_exp_20260625\';
% global experiment_name;
% experiment_name = 'exp1_basic';
% global maxItNumber;
% global diagZ;
% maxItNumber=50;
% global use_batteries;
% use_batteries = false;

for index_run = 1:1    
    R_num=0.03
    runNumber=1;
    probLevel=0.99;
    ptol=0.999;
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community1262\cplex\matlab\x64_win64');
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community1262\cplex\examples\src\matlab');
    savepath
    
    %Load parameters
    if use_batteries
        [fval,A,b,T,bT,xUp,xDn, Sigma, EX_, std_]=...
            DataloadDemandSideProblem3D_y_all(dimension,R_type,R_num,std_type, 0, [], [], []);
    else    
        [fval,A,b,T,bT,xUp,xDn, Sigma, EX_, std_]=...
            DataloadDemandSideProblem3D_cont(dimension,R_type,R_num,std_type, 0, [], [], []);
    end
    [zmin, zmax]=computeBrick(Sigma,ptol);

    diagZ=sqrt(sum((zmin-zmax).*(zmin-zmax)));

    m_vec=[250];
    for j=1:length(m_vec)
        t=zeros(runNumber,1);

        for i=1:runNumber
           tic;
           [x,Number_of_Genz_calls, Number_of_inner_iterations, pvec,grad_norm_objective,diff_sum_fi,timeCMax,x0_all]=DemandSideProblem(probLevel,i,m_vec(j),ptol,fval,A,b,T,bT,xUp,xDn, Sigma, zmin, zmax, R_num);

           %p(i,:)=pvec;
           %time(i)=toc;
           %timeCMaxs(i,:)=timeCMax;
        end
    end
end
