function P = DEgenerator(P, lu)
    % Use DE operator to optimiza the acquisition function
    %
    % INPUT:
    % - P   : parent population
    % - lu  : lower and upper bounds of variables
    %
    % OUTPUT
    % - P   : offspring population
    %
    
    % *** General settings here ****
    F  = 0.5;
    CR = 0.9;
    % ******************************

    [popsize, D] = size(P);
    for i = 1:popsize
        % Randomly select 3 solutions
        indexset     = 1:popsize;
        r1           = floor(rand*(popsize-1)) + 1;
        xr1          = indexset(r1);
        indexset(r1) = [];
        r2           = floor(rand*(popsize-2)) + 1;
        xr2          = indexset(r2);
        indexset(r2) = [];
        r3           = floor(rand*(popsize-3)) + 1;
        xr3          = indexset(r3);
        
        % Mutation
        v            = P(xr1,:) + F*( P(xr2,:) - P(xr3,:) );
        
        % Handle the elements of the mutant vector which violate the boundary
        vioLow       = find(v < lu(1, : ));
        v(1, vioLow) = 2.*lu(1, vioLow) - v(1, vioLow);
        vioLowUpper  = find(v(1, vioLow) > lu(2, vioLow));
        v(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
        vioUpper     = find(v > lu(2, : ));
        v(1, vioUpper) = 2 .* lu(2, vioUpper) - v(1, vioUpper);
        vioUpperLow  = find(v(1, vioUpper) < lu(1, vioUpper));
        v(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));

        % Binomial crossover
        jRand       = floor(rand * D) + 1;
        t           = rand(1, D) < CR;
        t(1, jRand) = 1;
        t_          = 1 - t;
        u           = t .* v + t_ .* P(i, : );
        P(i, : )    = u;
    end
end