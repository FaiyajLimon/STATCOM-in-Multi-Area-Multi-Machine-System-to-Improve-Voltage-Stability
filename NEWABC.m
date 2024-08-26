clear all
close all
clc

% Define the details of the table design problem
nVar = 6;                 % number of variables  
ub = [3 3 3 3 3 3];       % upper Bound
lb = [-3 -3 -3 -3 -3 -3]; % lower bound

fobj = @tunning;          % Objective function Name

% Define the ABC parameters 
noEmployedBees = 5;        % number of employed bees
maxIter = 50;             % maximum iterations
limit = 5;                % limit for abandoning solution
maxTrials = 50;          % maximum trials for scouting
nOnlookerBees = noEmployedBees; % number of onlooker bees (same as employed bees)

% Initialize the population (employed bees)
for i = 1 : noEmployedBees
    EmployedBees(i).X = (ub - lb) .* rand(1, nVar) + lb; 
    EmployedBees(i).O = fobj(EmployedBees(i).X);
    EmployedBees(i).trial = 0;
end

% Initialize the best solution (global best)
GBEST.X = EmployedBees(1).X;
GBEST.O = EmployedBees(1).O;

% ABC main loop
for iter = 1 : maxIter
    % Employed bees phase
    for i = 1 : noEmployedBees
        % Choose a random parameter to be modified
        paramToModify = randi(nVar);
        
        % Choose another solution different from i
        k = randi(noEmployedBees);
        while k == i
            k = randi(noEmployedBees);
        end
         
        % Generate a new solution based on neighbor information
        newSolution = EmployedBees(i).X;
        newSolution(paramToModify) = EmployedBees(i).X(paramToModify) + (rand - 0.5) * 2 * (EmployedBees(i).X(paramToModify) - EmployedBees(k).X(paramToModify));
        
        % Check if the new solution is within the boundaries
        newSolution = max(newSolution, lb);
        newSolution = min(newSolution, ub);
        
        % Evaluate the objective function
        newObj = fobj(newSolution);
        
        % Greedy selection
        if newObj < EmployedBees(i).O
            EmployedBees(i).X = newSolution;
            EmployedBees(i).O = newObj;
            EmployedBees(i).trial = 0;
        else
            EmployedBees(i).trial = EmployedBees(i).trial + 1;
        end
    end
    
    % Onlooker bees phase
    onlookerProbabilities = computeProbabilities(EmployedBees);
    onlookerIndex = 1;
    t = 0;
    while t < noEmployedBees
        if rand < onlookerProbabilities(onlookerIndex)
            t = t + 1;
            
            i = onlookerIndex;
            
            % Choose a random parameter to be modified
            paramToModify = randi(nVar);
            
            % Choose another solution different from i
            k = randi(noEmployedBees);
            while k == i
                k = randi(noEmployedBees);
            end
            
            % Generate a new solution based on neighbor information
            newSolution = EmployedBees(i).X;
            newSolution(paramToModify) = EmployedBees(i).X(paramToModify) + (rand - 0.5) * 2 * (EmployedBees(i).X(paramToModify) - EmployedBees(k).X(paramToModify));
            
            % Check if the new solution is within the boundaries
            newSolution = max(newSolution, lb);
            newSolution = min(newSolution, ub);
            
            % Evaluate the objective function
            newObj = fobj(newSolution);
            
            % Greedy selection
            if newObj < EmployedBees(i).O
                EmployedBees(i).X = newSolution;
                EmployedBees(i).O = newObj;
                EmployedBees(i).trial = 0;
            else
                EmployedBees(i).trial = EmployedBees(i).trial + 1;
            end
        end
        
        onlookerIndex = onlookerIndex + 1;
        if onlookerIndex > noEmployedBees
            onlookerIndex = 1;
        end
    end
    
    % Scout bees phase
    for i = 1 : noEmployedBees
        if EmployedBees(i).trial >= limit
            EmployedBees(i).X = (ub - lb) .* rand(1, nVar) + lb;
            EmployedBees(i).O = fobj(EmployedBees(i).X);
            EmployedBees(i).trial = 0;
        end
    end
    
    % Update the global best solution
    for i = 1 : noEmployedBees
        if EmployedBees(i).O < GBEST.O
            GBEST.X = EmployedBees(i).X;
            GBEST.O = EmployedBees(i).O;
        end
    end
    
    % Display iteration information
    disp(['Iteration #', num2str(iter), ', GBEST.O = ', num2str(GBEST.O)]);
    
    % Store GBEST.O for plotting
    cgCurve(iter) = GBEST.O;
end

% Plot convergence curve
semilogy(cgCurve);
xlabel('Iteration#');
ylabel('Weight');

% % Objective function (replace with your own function)
% function objective = tunning(variables)
%     % Objective function code here
%     objective = sum(variables);
% end
