clear all
close all
clc

nVar = 6;                                              % number of variables (Kp, Ki, Kd)
ub = [10 2000 50 50 50 50]; %upper Bound
lb = [0 0 0 0 0 0]; % lower bound

fobj = @tunning;            % objective function


% Define the HABC's parameters
noFood = 5;              % number of food sources
limit = nVar;               % limit of trials
maxIter =5;              % maximum iterations
% Define the PSO's parameters
noP = 5;                  % number of particles for initialization
maxIterPSO = 1;           % maximum iterations
wMax = 1;
wMin = 0.1;
c1 = 2;
c2 = 2;
vMax = (ub - lb) .* 0.2;
vMin = -vMax;
% Define the number of runs
numRuns =1;
% Main loop for multiple runs
for run = 1:numRuns
% Initialize the food sources
% Initialize the food sources
for i = 1 : noFood
    Foods(i).X = (ub - lb) .* rand(1, nVar) + lb;
    Foods(i).O = fobj(Foods(i).X);
    Foods(i).trial = 0; % Initialize trial value for each food source
   
    % Calculate fitness value for the food source
    Foods(i).fitness = 1 / (1 + Foods(i).O); % Modify this line based on your fitness function
end

% Initialize the best food source found so far
GBestFood.X = zeros(1, nVar);
GBestFood.O = inf;

% Initialize the counter for scout bee phase
scoutBeeCounter = zeros(noFood, 1);

% Initialize matrix to store trial values
trialValues = zeros(maxIter, noFood);
bestObjValuesPSO = zeros(1, noFood);  % Initialize the array to store best objective values from PSO

% Store the results of each run
meanValues = zeros(1, numRuns);
stdValues = zeros(1, numRuns);
bestObjValues = zeros(1, numRuns);% Initialize the new
% Initialize an array to store the objective function values for each iteration in all runs
objValuesAllRuns = zeros(maxIter, numRuns);


% Main loop
for t = 1 : maxIter
    % Employ the employed bees phase
    for i = 1 : noFood
            if Foods(i).trial <= limit
        % The parameter to be changed is determined randomly
        Param2Change = fix(rand * nVar) + 1;
       
        % A randomly chosen solution is used in producing a mutant solution of the solution i
        neighbour = fix(rand * noFood) + 1;
       
        % Randomly selected solution must be different from the solution i
        while neighbour == i
            neighbour = fix(rand * noFood) + 1;
        end
       
        sol = Foods(i).X;
       
        % v_{ij} = x_{ij} + \phi_{ij} * (x_{kj} - x_{ij})
    phi = -1 + 2 * rand; % Generate a random value between -1 and 1
     sol(Param2Change) = sol(Param2Change) + phi * (Foods(i).X(Param2Change) - Foods(neighbour).X(Param2Change));
        % If generated parameter value is out of boundaries, it is shifted onto the boundaries
        ind = find(sol < lb);
        sol(ind) = lb(ind);
        ind = find(sol > ub);
        sol(ind) = ub(ind);
       
        % Evaluate new solution
        ObjValSol = fobj(sol);
        FitnessSol = 1 / (1 + ObjValSol);
       
        % A greedy selection is applied between the current solution i and its mutant
        if FitnessSol > Foods(i).fitness % If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i
            Foods(i).X = sol;
            Foods(i).O = ObjValSol;
            Foods(i).fitness = FitnessSol;
            Foods(i).trial = 0;
        else
            Foods(i).trial = Foods(i).trial + 1; % If the solution i can not be improved, increase its trial counter
        end
       
        trialValues(t, i) = Foods(i).trial; % Store trial value
        
    end    
    
    end
    % Update the best food source found so far
    for i = 1 : noFood
        if Foods(i).O < GBestFood.O
            GBestFood = Foods(i);
        end
    end
  prob = ([Foods.fitness] ./ max([Foods.fitness]));


    % Employ the onlooker bees phase
    for i = 1 : noFood
 
        % Randomly select a food source to produce a mutant solution
        j = randi([1 noFood], 1);
        while j == i
            j = randi([1 noFood], 1);
        end
       
        if rand < prob(i)
            sol = Foods(i).X;
           
            % The parameter to be changed is determined randomly
            Param2Change = fix(rand * nVar) + 1;
           
            % A randomly chosen solution is used in producing a mutant solution of the solution i
            neighbour = fix(rand * noFood) + 1;
           
            % Randomly selected solution must be different from the solution i
            while neighbour == i
                neighbour = fix(rand * noFood) + 1;
            end
           
            % v_{ij} = x_{ij} + \phi_{ij} * (x_{kj} - x_{ij})
            phi = -1 + 2 * rand; % Generate a random value between -1 and 1
            sol(Param2Change) = sol(Param2Change) + phi * (Foods(i).X(Param2Change) - Foods(neighbour).X(Param2Change));
           
            % If generated parameter value is out of boundaries, it is shifted onto the boundaries
            ind = find(sol < lb);
            sol(ind) = lb(ind);
            ind = find(sol > ub);
            sol(ind) = ub(ind);
           
            % Evaluate new solution
            ObjValSol = fobj(sol);
            FitnessSol = 1 / (1 + ObjValSol);
           
            % A greedy selection is applied between the current solution i and its mutant
            if FitnessSol > Foods(i).fitness % If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i
                Foods(i).X = sol;
                Foods(i).O = ObjValSol;
                Foods(i).fitness = FitnessSol;
                Foods(i).trial = 0;
            else % The old food source is better
            end
           
        end
    end
   
    % Update the best food source found so far
    for i = 1 : noFood
        if Foods(i).O < GBestFood.O
            GBestFood = Foods(i);
        end
    end
   
    % Employ the scout bees phase
    for i = 1 : noFood
        if Foods(i).trial > limit
            scoutBeeCounter(i) = scoutBeeCounter(i) + 1; % Increment counter for the food source
            Foods(i).X = (ub - lb) .* rand(1, nVar) + lb;
            Foods(i).O = fobj(Foods(i).X);
            Foods(i).trial = 0; % Reset trial value after a scout bee phase
           
            % Calculate fitness value for the food source
            Foods(i).fitness = 1 / (1 + Foods(i).O); % Modify this line based on your fitness function
        end
    end


% Initialize the PSO particles using onlooker bee sources from ABC
Particles = struct('X', {}, 'V', {}, 'P', {}, 'O', {}, 'pBest', {}, 'OBest', {});
for i = 1 : noFood
    Particle.X = Foods(i).X;
    Particle.V = (ub - lb) .* rand(1, nVar) - (ub - lb) .* rand(1, nVar);
    Particle.P = Particle.X;
    Particle.O = Foods(i).O;
    Particle.pBest = Particle.X;
    Particle.OBest = Particle.O;
    Particles(i) = Particle;
end

% Employ the PSO optimization
for i = 1 : noFood
    Particle = Particles(i);
    X = Particle.X;
    V = Particle.V;
    P = Particle.P;
    O = Particle.O;
    pBest = Particle.pBest;
    OBest = Particle.OBest;
    for j = 1 : maxIterPSO
        w = wMax - (wMax - wMin) * j / maxIterPSO;
        V = w * V + c1 * rand(1, nVar) .* (pBest - X) + c2 * rand(1, nVar) .* (GBestFood.X - X);
        V = max(V, vMin);
        V = min(V, vMax);
        X = X + V;
        X = max(X, lb);
        X = min(X, ub);
        O = fobj(X);
        if O < OBest
            pBest = X;
            OBest = O;
        end
    end
    Particle.X = pBest;
    Particle.O = OBest;
    Particle.pBest = pBest;
    Particle.OBest = OBest;
    Particles(i) = Particle;
    if OBest < GBestFood.O
        GBestFood.X = pBest;
        GBestFood.O = OBest;
    end
end

    disp(['IterationPSO ' num2str(t) ': Best Objective Value = ' num2str(GBestFood.O)]);
   bestObjValuesPSO(t) = GBestFood.O;
   
  % Store the objective function value for the current iteration and run
        objValuesAllRuns(t, run) = GBestFood.O;
        
       
   
end
            bestObjValues(run) = min(bestObjValuesPSO);

disp(['Iteration ' num2str(maxIter) ': Best Objective Value = ' num2str(GBestFood.O)]);
%     disp('Scout Bee Phase Counter:');
%     disp(scoutBeeCounter);
% 
%     % Display trial values
%     disp('Trial Values:');
%     disp(trialValues);

    disp('Best Objective Values (PSO):');
    disp(bestObjValuesPSO);
    meanValue = mean(bestObjValuesPSO);
    stdValue = std(bestObjValuesPSO);

    disp(['Mean of Best Objective Values (PSO): ' num2str(meanValue)]);
    disp(['Standard Deviation of Best Objective Values (PSO): ' num2str(stdValue)]);
    
    % Store mean and standard deviation values
    meanValues(run) = meanValue;
    stdValues(run) = stdValue;
    
    disp('===================================');
       
end

% Calculate mean and standard deviation of all runs
meanAllRuns = mean(meanValues);
stdAllRuns = std(meanValues);

disp(['Mean of All Runs: ' num2str(meanAllRuns)]);
disp(['Standard Deviation of All Runs: ' num2str(stdAllRuns)]);

meanValue = mean(bestObjValues);
stdValue = std(bestObjValues);

disp(['Mean of Best Objective Values (PSO): ' num2str(meanValue)]);
disp(['Standard Deviation of Best Objective Values (PSO): ' num2str(stdValue)]);
