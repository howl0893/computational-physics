function [ output_args ] = PHYS_HW2( input_args )
% Matthew Howlett
%   A mass is released from rest with an initial position of 10 mm.  
%   This program is to make a single plot of the position as a function
%   of time for all of the trials, from t = 0 sec to t = 10 sec.


t = 0: .1 : 10;

gammaA = 0;
gammaB = 0.01;
gammaC = 0.30;
gammaD = 1;
gammaE = 2*pi;

trialA = 0.01*exp((-gammaA/2)*t)+((gammaA/2)*0.01)*exp((-gammaA/2)*t);
trialB = 0.01*exp((-gammaB/2)*t)+((gammaB/2)*0.01)*exp((-gammaB/2)*t);
trialC = 0.01*exp((-gammaC/2)*t)+((gammaC/2)*0.01)*exp((-gammaC/2)*t);
trialD = 0.01*exp((-gammaD/2)*t)+((gammaD/2)*0.01)*exp((-gammaD/2)*t);
trialE = 0.01*exp((-gammaE/2)*t)+((gammaE/2)*0.01)*exp((-gammaE/2)*t);

plot(t, trialA, '-r','LineWidth',2); hold on;
plot(t, trialB, '-b','LineWidth',2); 
plot(t, trialC, '-m','LineWidth',2);
plot(t, trialD, '-g','LineWidth',2);
plot(t, trialE, '-k','LineWidth',2);
hold off;

end

