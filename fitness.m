function [objF, conV] = fitness(P, problem)
    popsize = size(P, 1);
    g = [];
    % g denotes the constraint
    % f denotes the objective function
    switch problem
        case 1
            % G03 problem
            g(:, 1) = sum(P.^2, 2) - 1 - 0.0001;
            g(:, 2) = -(sum(P.^2, 2) - 1) - 0.0001;
            f       = -(2.^0.5)^2 * prod(P, 2);

        case 2
            % G06 problem
            g(:, 1) = -(P(:, 1) - 5).^2 - (P(:, 2) - 5).^2 + 100;
            g(:, 2) = (P(:, 1) - 6).^2 + (P(:, 2) - 5).^2 - 82.81;
            f       = (P(:, 1) - 10).^3 + (P(:, 2) - 20).^3;

        case 3  
            % G09 problem
            g(:, 1) = -127 + 2 * P(:, 1).^2 + 3 * P(:, 2).^4 + P(:, 3) + 4 * P(:, 4).^2 + 5 * P(:, 5);
            g(:, 2) = -282 + 7 * P(:, 1) + 3 * P(:, 2) + 10 * P(:, 3).^2 + P(:, 4) - P(:, 5);
            g(:, 3) = -196 + 23 * P(:, 1) + P(:, 2).^2 + 6 * P(:, 6).^2 - 8 * P(:, 7);
            g(:, 4) = 4 * P(:, 1).^2 + P(:, 2).^2 - 3 * P(:, 1).* P(:, 2) + 2 * P(:, 3).^2 + 5 * P(:, 6) - 11 * P(:, 7);
            f       = (P(:, 1) - 10).^2 + 5 * (P(:, 2) - 12).^2 + P(:, 3).^4 + 3 * (P(:, 4) - 11).^2 + 10 * P(:, 5).^6 + ...
                7 * P(:, 6).^2 + P(:, 7).^4 - 4 * P(:, 6).* P(:, 7) - 10 * P(:, 6) - 8 * P(:, 7);

        case 4
            % G11 problem
            g(:, 1) = (P(:, 2) - P(:, 1).^2) - 0.0001;
            g(:, 2) = -(P(:, 2) - P(:, 1).^2) - 0.0001;
            f       = P(:, 1).^2 + (P(:, 2) - 1).^2;

        case 5
            % G12 problem
            l = 1;
            for i = 1:9
                for j = 1:9
                    for k = 1:9
                        aaa(l, :) = [i j k];
                        l = l+1;
                    end
                end
            end 
            f = -(100 - (P(:, 1) - 5).^2 - (P(:, 2) - 5).^2 - (P(:, 3) - 5).^2)/100;
            for j = 1:popsize
                g(j, 1) = min(sum((repmat(P(j, :), 9 * 9 * 9, 1) - aaa).^2, 2)) - 0.0625;
            end

        case 6
            % G24 problem
            g(:, 1) = -2 * P(:, 1).^4 + 8 * P(:, 1).^3 - 8 * P(:, 1).^2 + P(:, 2) - 2;
            g(:, 2) = -4 * P(:, 1).^4 + 32 * P(:, 1).^3 - 88 * P(:, 1).^2 + 96 * P(:, 1) + P(:, 2) - 36;
            f       = -P(:, 1) - P(:, 2);
            
        case 7   
            % Pressure vessel design (PVD) problem
            x1     = P(:,1); x2 = P(:,2); x3 = P(:,3); x4 = P(:,4);
            f      = 0.6224.*x1.*x3.*x4 + 1.7781.*x2.*x3.^2 + 3.1661.*x1.^2.*x4 + 19.84.*x1.^2.*x3;
            g(:,1) = -x1 + 0.0193.*x3;
            g(:,2) = -x2 + 0.00954.*x3;
            g(:,3) = -pi.*x3.^2.*x4 - 4/3.*pi.*x3.^3 + 1296000;
            g(:,4) = x4 - 240;
       
        case 8  
            % Tension/compression spring design (TSD) problem
            x1     = P(:,1); x2 = P(:,2); x3 = P(:,3); 
            f      = x1.^2.*x2.*(x3+2);
            g(:,1) = 1-x2.^3.*x3./(71785*x1.^4);
            g(:,2) = (4*x2.^2-x1.*x2)./(12566.*x1.^3.*(x2-x1)) + 1./(5108.*x1.^2) - 1;
            g(:,3) = 1-140.45.*x1./(x3.*x2.^2);
            g(:,4) = (x1+x2)/1.5 - 1;
            
        case 9 
            % Three-bar Truss Design (TBTD) Problem
            x1     = P(:,1); x2 = P(:,2); l=100; p = 2; sigma = 2; 
            f      = l*(x2+2*sqrt(2).*x1);
            g(:,1) = x2.*p./(2*x2.*x1+sqrt(2).*x1.^2)-sigma;
            g(:,2) = (x2+sqrt(2).*x1)*p./(2*x2.*x1+sqrt(2).*x1.^2)-sigma;
            g(:,3) = p./(x1+sqrt(2)*x2)-sigma;
    end

    % Obtain the fitness
    objF = f;
    conV = g;
end