function [fval,A0,b0,T,bT,xUp,xDn,R,  EX_, std_] = DataloadDemandSideProblem3D_y_all(dimension_, R_type, R_num, std_type, plotting, x0, x0_all, p)

global save_directory;
global experiment_name;
global maxItNumber;
global cons_op;

%power of a single controllable applience
h = 1000 %W
%power of chargeable/dischargeable batteries
s = 500 %W
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
%number of batteries
Y = 20
%initial number of batteries
y0 = 20
%upper limit of cost	
d = (dimension_/24)*1000

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
fval = repmat(0,2*N,1)';

%A -> a vektor
A = repmat(alpha,2*N,1)';

%T = T_eredeti/std
T = [diag(repmat(-s,1,N)./std_) diag(repmat(-h,1,N)./std_)];

%bT = 0 - Mu_/std
bT = repmat(0,N,1) - (Mu_./std_)';

%korlát
Y_vec = repmat(Y,1,N)
y0_vec = repmat(y0,1,N);
%korlát
Y_y0_vec = repmat(Y-y0,1,N);

%
ones_ = repmat(1,1,N);
zeros_= repmat(0,1,N);

%felsõ korlát
%            y              x
xUp = [repmat(Y,1,N) repmat(M,1,N)]
%            y              x
%alsó korlát
xDn = [repmat(-Y,1,N) repmat(0,1,N)] 

%y-ok tölthetõ/kisüthetõ akksik
L = tril(ones(N));

A0=[A;
    [L 0*L]; 
    [-L 0*L];
    [eye(N) 0*eye(N)];
    [-1*eye(N) 0*eye(N)];
    [zeros_ ones_]; 
    [zeros_ -ones_];
    ]

b0=[b; 
    Y_y0_vec';
    y0_vec';
    Y_vec';  
    Y_vec';
    M*k;
    -M*k; 
    ]
b
if plotting

    %X Plot
    gcf = figure;
    y=x0(1:dimension_)'
    x=x0(dimension_+1:end)'
    bar(x*h,'g','BarWidth',0.8)
    hold on
    bar(y*s,'y','BarWidth',0.6)
    grid on
    xlabel('time windows')
    ylabel('discharge/charge/consumption [W]')
    legend2 = legend('batteries','controllable appliances')
    set(legend2,...
    'Position',[0.269615186458728 0.671140797121318 0.319259259259259 0.184429819635678]);
    periods = 1:N;
    time_unit = 1; %0.5; %30 minutes
    errhigh = time_unit*std_;
    errlow  = time_unit*std_;
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_x'];
    savefig([save_str,'.fig']);
    saveas(gcf, [save_str,'.eps']);

    %Bar Chart
    %https://www.mathworks.com/help/matlab/creating_plots/bar-chart-with-error-bars.html
    gcf = figure;
    time_unit = 1;
    errhigh = time_unit*std_;
    errlow  = time_unit*std_;

    bar(periods,time_unit*e,'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625])
    hold on
    bar(periods,time_unit*(Mu+x*h+y*s),'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911])
    hold on
    bar(periods,time_unit*(Mu+y*s), 'FaceColor',[0.39215686917305 0.474509805440903 0.635294139385223])
    hold on
    er = errorbar(periods,time_unit*(Mu+x*h+y*s),errlow,errhigh, 'LineWidth',2, 'Color',[0.152941182255745 0.227450981736183 0.372549027204514]);
    %er = errorbar(periods,Mu,errlow,errhigh, 'LineWidth',2);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    xlabel('time windows')
    ylabel('consumption [W]')
    legend1 = legend('upper limit', 'controllable devices','random consumption', 'std')
    %title('iterations')
    set(legend1,...
    'Position',[0.183072574762946 0.155704473775105 0.298333333333333 0.187993680884676]);
    grid on
    save_str = [save_directory,experiment_name,'_dim=', num2str(dimension_),'_M=',num2str(M),'_it=',num2str(maxItNumber),'_',R_type,'_',std_type,'_consumption'];
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
end

e-Mu-std
sum(e-Mu-std_)
sum(e-Mu+std_)
M*k*h