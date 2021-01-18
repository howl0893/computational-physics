%% ForwardEulerMethod.m
%
% Algorithm designed to implement the forward euler method to predict the exponential 
% cell growth of a bacteria
%% Declaration of Parameters

t = 0:0.1:4.8;     %2 days
dt = 0.1;         %time step
N = ones(1,49);   %create vector for cells
N(1) = 5;
a = N(1);
b = 2;             %positive growth factor constant
T = 1;       %time interval for x to increase by a factor of b
tau = T/log(2);    %time constant
%% Calculations

%exact solution
x = a.*b.^(t/T);

%Euler forward method
for i = (1: 1: 48)
    t(i+1) = t(i) + dt;           %increase time by step size each iteration
    N(i+1) = N(i) * (1 + dt/tau); %calculate new number of cells each iteration
end

%% Plotting

figure(1);
hold on;

plot(t,N);                      %plot the theoretical solution
plot(t,x);                      %plot the exact solution

xlabel('Time [h/20]');             %label x-axis
ylabel('Number of cells');      %label y-axis
title('Cell growth over time'); %label title

xlim([0 2.4]);

legend('Theoretical','Exact');

hold off;