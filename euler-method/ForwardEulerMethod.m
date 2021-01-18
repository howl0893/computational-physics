%% ForwardEulerMethod.m
% Algorithm designed to implement the forward euler method to predict the exponential 
% cell growth of a bacteria

%% Declaration of Parameters

dt = 0.05;        % step time
length = (48/dt);
N(1) = 5;       % number of cells
a = N;
t = linspace(0, dt, 961);
b = 2;
T = (40/60);       % time [h] for cells to double
tau = T/log(2);   % time constant

size(t)
size(N)
%exact solution
x =a.*b.^(t/T);

%% Calculations

%Euler forward method
for i = 1:length
    %t(i+1) = t(i) + dt;          %increase time by .1 each iteration
    N(i+1) = (N(i)*(1 + dt/tau)); %calculate new number of cells each iteration
end

%% Plotting

figure(1);
hold on                         %hold to include informaiton on graph

plot(t,N);                      %plot the theoretical solution
plot(t,x);                      %plot the exact solution

xlabel('Time [days]');          %label x-axis
ylabel('Number of cells');      %label y-axis
title('Cell growth over time'); %label title

xlim([0 2]);

legend('Theoretical Solution', 'Exact Solution'); % create legend

hold off