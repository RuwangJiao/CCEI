function minF = CCEI(problemIndex)
    %
    % INPUT:
    % - problemIndex: the problem to be solved
    %
    % OUTPUT
    % - minF   : the optimized objective function value
    %
    
    % *** General settings here ****
    maxFEs = 100; % Maximum number of function evaluations
    nTest  = 30;  % Population size of DE to optimize the acquisition function
    maxGen = 200; % Maximum number of generations to optimize the acquisition function
    % ******************************
    
    format long
    format compact
    addpath(genpath('scripts'));

    [nO, nC, nD, lu] = problem(problemIndex);
        
    % Latin hypercube design
    rand('seed', sum(100*clock));
    iniSize = 11*nD - 1;
    x_train = lhsdesign(iniSize, nD, 'criterion','maximin', 'iteration',100);
    P       = repmat(lu(1,:), iniSize, 1) + x_train.*repmat((lu(2,:) - lu(1,:)), iniSize, 1);
    [objF, conV] = fitness(P, problemIndex);
    y_train = [objF, conV];
        
    % Normalizaion (data pre-processing)
    maxY    = max(abs(y_train)); 
    trainY  = y_train./maxY;

    FEs     = 1;
    nTask   = nO + nC;
    minF = inf; minNorG = inf; minNorF = inf;
    boolFeasible = max(max(0, conV), [], 2) == 0;
    if size(objF(boolFeasible), 1) == 0   
        % Infeasible case
        feasibleFlag = 0;
        G            = trainY(:, 2:end);
        minNorG      = min(max(max(0, G), [], 2), [], 1);
        [model, data] = Training(x_train, G(:), nTask-1);
    else
        % Feasible case
        feasibleFlag = 1;
        minF         = min(objF(boolFeasible), [], 1);
        F            = trainY(:,1);
        minNorF      = min(F(boolFeasible), [], 1);
        [model,data] = Training(x_train, trainY(:), nTask);
    end
        
    while FEs <= maxFEs   
        x_test = rand(nTest, nD); 
        EI_parent = calEI(nTask, x_test, feasibleFlag, nC, minNorF, minNorG, model, data);
        for g = 1:maxGen
            x_child  = DEgenerator(x_test, [zeros(1,nD); ones(1,nD)]);
            EI_child = calEI(nTask, x_child, feasibleFlag, nC, minNorF, minNorG, model,data);
            EI_parent(EI_child>=EI_parent) = EI_child(EI_child>=EI_parent);
            x_test(EI_child>=EI_parent,:)  = x_child(EI_child>=EI_parent,:);
        end 
        [~, max_index] = max(EI_parent);
        ind           = repmat(lu(1,:), size(x_test(max_index,:), 1), 1) + x_test(max_index,:).*repmat((lu(2,:) - lu(1,:)), size(x_test(max_index,:), 1), 1);
        [indf, indv]  = fitness(ind, problemIndex);
        x_train       = [x_train; x_test(max_index,:)];
        objF          = [objF; indf];
        conV          = [conV; indv];
        FEs           = FEs + 1;
        y_train       = [objF, conV];
        % Normalizaion
        maxY    = max(abs(y_train)); 
        trainY  = y_train./maxY;
            
        % Judge whether there are feasible solutions in the database
        boolFeasible = max(max(0, conV), [], 2) == 0;
        if size(objF(boolFeasible), 1) == 0  
            % Infeasible case
            feasibleFlag = 0;
            G = trainY(:, 2:end);
            minNorG = min(max(max(0, G), [], 2), [], 1);
            [model,data]  =  Training(x_train, G(:), nTask-1);
        else
            % Feasible case
            feasibleFlag = 1;
            minF = min(objF(boolFeasible), [], 1);
            F = trainY(:,1);
            minNorF = min(F(boolFeasible), [], 1);
            [model,data]  =  Training(x_train, trainY(:), nTask);
        end
    end
end