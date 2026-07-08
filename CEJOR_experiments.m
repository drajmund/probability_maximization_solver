dimension = 24;
global R_type;
global std_type;
global save_directory;
save_directory = '.\CEJOR_exp_20260625\';
global experiment_name;
global maxItNumber;
global diagZ;
maxItNumber=50;
global use_batteries;
global cons_op;
cons_op = 1;

%Model I Basic - Exp1
R_type = 'decr';
std_type = 'constant';
experiment_name = 'exp1_basic';
use_batteries = false;

tic, callDmandSideProblem, toc

%Model I Basic - Exp2
% R_type = 'no_corr';
% std_type = 'increasing';
% experiment_name = 'exp2_basic';
% use_batteries = false;
% 
% tic, callDmandSideProblem, toc

%Model II Battery - Exp1
% R_type = 'no_corr';
% std_type = 'constant';
% experiment_name = 'exp1_battery';
% use_batteries = true;
% 
% tic, callDmandSideProblem, toc

%Model II Battery - Exp2
% R_type = 'no_corr';
% std_type = 'increasing';
% experiment_name = 'exp2_battery';
% use_batteries = true;
% 
% tic, callDmandSideProblem, toc


%Model III Consecutive Operation - Exp1
% maxItNumber=100;
% R_type = 'no_corr';
% std_type = 'constant';
% experiment_name = 'exp1_cons';
% use_batteries = false;
% cons_op = 10;
% 
% tic, callDmandSideProblem, toc

%Model III Consecutive Operation - Exp2
% maxItNumber=100;
% R_type = 'no_corr';
% std_type = 'increasing';
% experiment_name = 'exp2_cons';
% use_batteries = false;
% cons_op = 10;
% tic, callDmandSideProblem, toc

%48 dim - Exp1
% dimension = 48;
% maxItNumber=100;
% R_type = 'decr';
% std_type = 'increasing';
% experiment_name = 'exp1_48';
% use_batteries = true;
% cons_op = 1;
% tic, callDmandSideProblem, toc
