function [Mu, e, std] = generateMuEStd(dimension_, type)

sequence = [1:dimension_];
s1 = -abs((sequence - floor(length(sequence)/3)).^1.3);
Mu = ones(1,dimension_)*20000 + (s1 - min(s1))*1000;
e = ones(1,dimension_)*max(Mu).*1.5;
s1 = -abs(sequence - floor(length(sequence)/3));

if strcmp(type,'peak')
    std = ones(1,dimension_)*2000 + (s1 - min(s1))*1000;
elseif strcmp(type,'constant')
    std = ones(1,dimension_)*sum(sequence*1000)/dimension_;
elseif strcmp(type,'increasing')
    std_const = ones(1,dimension_)*sum(sequence*1000)/dimension_;
    std = sequence*300+(sum(std_const)-sequence*300)/dimension_;
elseif strcmp(type,'saw_constant')
    std = ones(1,dimension_)*sum(sequence*1000)/dimension_;
    Mu_x = 2*sum(Mu)/(dimension_*1.5);
    Mu(1:2:end) = Mu_x;
    Mu(2:2:end) = Mu_x/2;
elseif strcmp(type,'saw_peak')
    std = ones(1,dimension_)*2000 + (s1 - min(s1))*1000;
    Mu_x = 2*sum(Mu)/(dimension_*1.5);
    Mu(1:2:end) = Mu_x;
    Mu(2:2:end) = Mu_x/2;
end

 