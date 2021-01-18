% Unit 5 Homework 3
% Joshua Wylie

clear

x0 = [0]; %setting an initial zero value

nWalkers = 2; % sets the number of walkers
stepMax = 5; % initial number of steps

for stepMultiple = 1:1:3 % will increase the number of steps
    stepSize = 10; % in nm
    tmax = 50; % in s
    dt = 1; % in s
    time = 0; % in s
    x = x0; % in nm
    
    for j = 1:nWalkers
        
        for t = 1:1:stepMax
            
            random = rand();
            if rem((t/6),1) == 0 % find if the step number is a multiple of 6
                random = -random; %moves in the opposite direction at a multiple of 6
            end
            magnitude = random*stepSize; % determines how far the Myosin moves
            
            x(t+1) = x(t) + magnitude; % moves the Myosin
            
            time(t+1) = time(t) +dt; % increment time/ step number
        end
        displacement(j) = x(t) - x(1); % finds the displacement for each walker
    end
    
    displacementAvg(stepMultiple) = sum(displacement)/nWalkers; % finds the average displacement
    numberOfSteps(stepMultiple) = stepMax; %saves the maximum step number
    
%     plots position with respect to time/step number
    figure; plot(time,x)
    title(['Myosin After ' num2str(stepMax) ' Steps'])
    xlabel('Number of Steps')
    ylabel('X-Position of Myosin')
    legend('Myosin')
    
    stepMax = stepMax + 5*stepMultiple;
end

% plots the average displacement with respect to the maximum steps
figure; plot(numberOfSteps,displacementAvg)
title(['Average Displacement vs Number of Steps'])
xlabel('Number of Steps')
ylabel('Displacement')
legend('Myosin Displacement')
