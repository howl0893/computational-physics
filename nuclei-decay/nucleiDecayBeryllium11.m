%% Beryllium-11 Nuclei Decay
%
% This algorithm is design to model the decay of beryllium-11 nuclei with a
% starting amount of 5000 total nuclei.

%% Parameters
N0 = 5000;                  % initial number of nuclei
N = zeros(1,N0);            % create vector to store decaying nuclei
halfLife = 13.81;           % half life of Beryllium-11 [seconds]
multiplier = 10;            % Number of half lifes being evaluated
timestep = 0.1;             % days
lambda = log(2)/halfLife;   % decay constant
nucleiN = 1:N0;             % used to plot nuclei
p = lambda*timestep;        % probability of decay
time = 0:timestep:(halfLife*multiplier);  % time elapsed 
index = 0;                  % used for indexing number of undecayed matrix
pIndex = 0;                 % used for indexing percent error matrix

%% Calculations
% Theoretical Solution
exactSol=N0.*exp(-lambda.*time);

% Monte Carlo Approximation
for t = 0:timestep:(halfLife*multiplier)
    for i = 1:N0
        r = rand;  
        if r < p && N(i)~= 1
            N(i) = 1;
%             plot(nucleiN,N,'bo');
%             pause
        end
    end
    index=index+1;
    Nundecayed(index)=N0-sum(N);
     
    % calculate percent error at each point for first 5 half lives
    if t <= (halfLife*5)
        percentErr(index) = abs((Nundecayed(index)-exactSol(index))/exactSol(index)).*100;
        pIndex = pIndex+1;
    end
end

% Average Percent Error
avgPercentErr = abs((sum(Nundecayed)-sum(exactSol))./sum(exactSol)).*100;
disp('Average Percent Error = ')
disp(avgPercentErr)

%% Plotting
figure(1);
hold on;
plot(time, Nundecayed);
plot(time,exactSol);
xlim([0 (halfLife*multiplier)]);
titleStr = ['Non-decayed Beryllium Nuclei Over ', num2str(multiplier),' Half Lives'];
title(titleStr);
xlabel('Time elapsed [s]');
ylabel('Number of Undecayed Nuclei');
legend('Monte Carlo solution','Analytical solution');

figure(2);
hold on;
plot(1:pIndex,percentErr);
title('Percent Error at Each Point for First 5 Half Lives');
xlabel('Point');
ylabel('Percent Error [%]');
