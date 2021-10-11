function [nO, nC, n, lu] = problem(problem_index)
    % load the propery of the test problems
    %
    % INPUT:
    % - problem_index: the index of the test problem
    %
    % OUTPUT
    % - nO   : number of objective
    % - nC   : number of constraint
    % - n    : number of decision variable
    % - lu   : the lower and upper bounds of the decision variable
    % ******************************
    
    switch problem_index
        case 1   
            % G03 problem
            n  = 2;
            lu = [zeros(1, n); ones(1, n)]; %\\\-1\ 
            nO = 1;
            nC = 2;
        case 2  
            % G06 problem
            lu = [13 0; 100 100];
            n  = 2;  %\\-6961.81388\  
            nO = 1;
            nC = 2;
        case 3
            % G09 problem
            lu = [-10 -10 -10 -10 -10 -10 -10; 10 10 10 10 10 10 10];
            n  = 7; %\680.6300573\   
            nO = 1;
            nC = 4;
        case 4 
            % G11 problem
            lu = [-1 -1; 1 1];
            n  = 2;  %\\0.75\
            nO = 1;
            nC = 2;
        case 5 
            % G12 problem
            lu = [0 0 0; 10 10 10];
            n  = 3; %\\-1\
            nO = 1;
            nC = 1;
        case 6
            % G24 problem
            lu = [0 0;  3 4];
            n  = 2; %\-5.5080132716\
            nO = 1;
            nC = 2;
        case 7
            % Pressure vessel design (PVD) problem
            lu = [0.0625 0.0625 10 10; 99*0.0625 99*0.0625 200 200];
            n  = 4; 
            nO = 1;
            nC = 4;
        case 8
            % Tension/compression spring design (TSD) problem
            lu = [0.05 0.25 2; 2 1.3 15];
            n  = 3; 
            nO = 1;
            nC = 4;
        case 9
            % Three-bar Truss Design (TBTD) Problem
            lu = [0 0; 1 1];
            n  = 2; 
            nO = 1;
            nC = 3;
    end
end