clc
clear all
format long
format compact
tic

addpath(genpath('scripts'));
warning('off');

problemSet = [1:44];
maxFEs   = 100;
for problemIndex = [44]
    prob      = problemSet(problemIndex);
    time      = 1;
    totalTime = 30;  
    while time <= totalTime 
        [nO, nC, nD, lu] = problem(prob);
        iniSize = 11*nD - 1;
        rand('seed', sum(100*clock));
        x_train = lhsdesign(iniSize, nD, 'criterion','maximin', 'iteration',100);
        P       = repmat(lu(1,:), iniSize, 1) + x_train.*repmat((lu(2,:) - lu(1,:)), iniSize, 1);
        [objF, conV] = fitness(P, prob);
        y_train = [objF, conV];
        
        % Normalizaion (data pre-processing)
        maxY    = max(abs(y_train)); 
        trainY  = y_train./maxY;

        FEs     = 1;
        nTask   = nO + nC;
        minG = inf; minF = inf; minNorG = inf; minNorF = inf;
        boolFeasible = max(max(0, conV), [], 2) == 0;
        if size(objF(boolFeasible), 1) == 0   % Infeasible case
            feasibleFlag = 0;
            minG         = min(max(max(0, conV), [], 2), [], 1);
            G            = trainY(:, 2:end);
            minNorG      = min(max(max(0, G), [], 2), [], 1);
            disp(['Infeasible! Minimum CV: ', num2str(minG)]);
            [model, data] = Training(x_train, G(:), nTask-1);
        else                                  % Feasible case
            feasibleFlag = 1;
            minF         = min(objF(boolFeasible), [], 1);
            minG         = 0;
            F            = trainY(:,1);
            minNorF      = min(F(boolFeasible), [], 1);
            disp(['Feasible! Best obj: ', num2str(minF)]);
            [model,data] = Training(x_train, trainY(:), nTask);
        end
        %Save data
        savePath = strcat('./Data/', 'g', num2str(problemIndex), '-', num2str(time), '.txt');
        fid      = fopen(savePath, 'wt');
        fprintf(fid, '%g\t', FEs);
        fprintf(fid, '%g\t', minF);  
        fprintf(fid, '%g\n', minG);
        while FEs <= maxFEs   
            disp(['Number of function evaluations: ', num2str(FEs)]);
            nTest  = 30;
            x_test = rand(nTest, nD);
            
            P1       = repmat(lu(1,:),nTest,1) + x_test.*repmat((lu(2,:) - lu(1,:)), nTest, 1);
            [objF1, conV1] = fitness(P1, prob);
            
            EI_parent = calEI(nTask, x_test, feasibleFlag, nC, minNorF, minNorG, model, data);
            for g = 1:200
                x_child  = DEgenerator(x_test, [zeros(1,nD); ones(1,nD)]);
                EI_child = calEI(nTask, x_child, feasibleFlag, nC, minNorF, minNorG, model,data);
                EI_parent(EI_child>=EI_parent) = EI_child(EI_child>=EI_parent);
                x_test(EI_child>=EI_parent,:)  = x_child(EI_child>=EI_parent,:);
            end 
            [~, max_index] = max(EI_parent);
            %disp(max(EI_parent));
            ind           = repmat(lu(1,:), size(x_test(max_index,:), 1), 1) + x_test(max_index,:).*repmat((lu(2,:) - lu(1,:)), size(x_test(max_index,:), 1), 1);
            [indf, indv]  = fitness(ind, prob);
            x_train       = [x_train; x_test(max_index,:)];
            objF          = [objF; indf];
            conV          = [conV; indv];
            FEs           = FEs + 1;
            y_train       = [objF, conV];
            disp(['Predict obj: ', num2str(indf), ';  Predict CV: ', num2str(max(max(0, indv), [], 2))]);
            % Normalizaion
            maxY    = max(abs(y_train)); 
            trainY  = y_train./maxY;
            % Judge whether there are feasible solutions in the database
            boolFeasible = max(max(0, conV), [], 2) == 0;
            if size(objF(boolFeasible), 1) == 0  % Infeasible case
                feasibleFlag = 0;
                minG = min(max(max(0, conV), [], 2), [], 1);
                G = trainY(:, 2:end);
                minNorG = min(max(max(0, G), [], 2), [], 1);
                disp(['Infeasible! Minimum CV: ', num2str(minG)]);
                [model,data]  =  Training(x_train, G(:), nTask-1);
            else                                 % Feasible case
                feasibleFlag = 1;
                minF = min(objF(boolFeasible), [], 1);
                minG = min(max(max(0, conV), [], 2), [], 1);
                F = trainY(:,1);
                minNorF = min(F(boolFeasible), [], 1);
                disp(['Feasible! Best obj: ', num2str(minF)]);
                [model,data]  =  Training(x_train, trainY(:), nTask);
            end
            fprintf(fid, '%g\t', FEs);
            fprintf(fid, '%g\t', minF);  
            fprintf(fid, '%g\n', minG);
        end
        time = time + 1;
        fclose(fid);
    end
end