function [fval,A0,b0,T,bT,xUp,xDn,R,  EX_, std_] = DataloadDemandSideProblem3D_cont(dimension_, R_type, R_num, std_type, plotting, x0, x0_all, p)

global save_directory;
global experiment_name;
global maxItNumber;
global cons_op;

%power of a single controllable applience	
h = 1000 %W
%number of time slots for consecutive operation/charging			
k = cons_op
%parameters of the cost fn.			
alpha = 1.2 %(slope)		
beta = 2 % (b)	

R = generateKorrelationMatrix(dimension_, R_type, R_num);

%Expected value, upper limit, standard deviations of each variable
[Mu, e, std] = generateMuEStd(dimension_, std_type);

%number of controllable devices	
M = (dimension_/24)*200/k

%upper limit of cost	
d = (dimension_/24)*1000 %24, 12->2000

%Expected value shifted
Mu_ = Mu - e;
EX_ = Mu_;

D = diag(std); %dim x dim mátrix

% kovariancia mátrix
Sigma = D*R*D;
std_=sqrt(diag(Sigma))';
%-------------------------------------

%idõszakok száma		
N = size(std_,2); %idõszak

%helyettesítõ		
b = d - N*beta;

%-------------------------------------
fval = repmat(0,N,1)';

%C=diag(ones(N-2,1),-2)+diag(ones(N-1,1),-1)+diag(ones(N,1));
%C=diag(ones(N,1));
C = diag(ones(N,1));
for i=2:k
    C=C+diag(ones(N-(i-1),1),-(i-1));
end

%A -> a vektor
A = (C*repmat(alpha,N,1))';

%T = T_eredeti/std
T = C*diag(repmat(-h,1,N)./std_);


%bT = 0 - Mu_/std
bT = repmat(0,N,1) - (Mu_./std_)';

%felsõ korlát
xUp = repmat(M,1,N);
%alsó korlát
xDn = repmat(0,1,N); 

%1*x = M*k
ones_ = repmat(1,1,N);
A0=[A;
	ones_*C; 
    -ones_*C;
    C;
    -1*C]

b0=[b; 
    M*k;
    -M*k; 
    xUp';  
    -xDn']

if plotting
    x = (C*x0)'
    
    %Bar Chart
    %https://www.mathworks.com/help/matlab/creating_plots/bar-chart-with-error-bars.html
    periods = 1:N;
    time_unit = 1 %0.5; %30 minutes
    errhigh = time_unit*std_;
    errlow  = time_unit*std_;

    gcf = figure;
    bar(periods,time_unit*e,'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625])
    hold on
    bar(periods,time_unit*(Mu+x*h),'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911])
    hold on
    bar(periods,time_unit*Mu, 'FaceColor',[0.39215686917305 0.474509805440903 0.635294139385223])
    hold on
    er = errorbar(periods,time_unit*(Mu+x*h),errlow,errhigh, 'LineWidth',1, 'Color',[0.152941182255745 0.227450981736183 0.372549027204514]);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    xlabel('time windows')
    ylabel('consumption [W]')
    legend1 = legend('upper limit', 'controllable appliances','random consumption', 'std of random consumption')
    set(legend1,...
        'Position',[0.158630952380949 0.127380952380952 0.319642857142857 0.188888888888889]);
    grid on
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_consumption'];
    savefig([save_str,'.fig']);
    saveas(gcf, [save_str,'.eps']);
    hold off

    %X Plot
    gcf = figure;
    bar(C*x0*h) %bar(x0*h/2)
    xlabel('time windows')
    ylabel('consumption [W]')
    grid on
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_x'];
    savefig([save_str,'.fig']);
    saveas(gcf, [save_str,'.eps']);
    
    %Probablity Plot
    gcf = figure;
    plot(p','b-o','LineWidth',3)
    xlabel('iterations')
    ylabel('probability')
    grid on
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_p'];
    savefig([save_str,'.fig']);
    saveas(gcf, [save_str,'.eps']);


    gcf = figure;
    surf(x0_all)
    grid on
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_3d'];
    savefig([save_str,'.fig']);
    saveas(gcf, [save_str,'.eps']);


end

e-Mu-std
sum(e-Mu-std_)
sum(e-Mu+std_)
M*k*h
C