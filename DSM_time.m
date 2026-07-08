clear all
%global maxItNumber;
%it = [5 10 20 50 100 200]
time_=[1 2 3 4 5];
for n=1:5
    %maxItNumber = it(n);
    tic;
    callDmandSideProblem;
    time_(n)=toc;
end