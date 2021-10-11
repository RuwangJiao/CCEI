function EIvalue = calEI(nTask, x_test, feasibleFlag, nC, minNorF, minNorG, Model, data)
    % Calculcate constrained expected improvement value (acquisition function)
    %
    % INPUT:
    % - nTask       : number of tasks
    % - x_test      : test points
    % - feasibleFlag: whether the current samples include feasible points 
    % - nC          : number of constraints
    % - minNorF     : the minimum normilized objective value
    % - minNorG     : the minimum normilized constraint violation value
    % - Model       : hyperparameters of MTGP model
    % - data        : cell data for learning and prediction
    %
    % OUTPUT
    % - EIvalue     : constrained expected improvement value
    %
    % Author        : Ruwang Jiao

    format long;
    EIvalue = zeros(size(x_test, 1), 1);
    n = 10000;             % The sample number of Monte Carlo 
    %syms f g1 g2 g3 z;
    if feasibleFlag == 1   % feasible situation
        [ymu, ys2] =  predict_mtgp_all_tasks(Model, data, x_test);
        for i = 1:size(ys2,1)
            mu    = ymu(i,:);
            sigma = reshape(ys2(i,:), nTask, nTask);
            sigma = tril(sigma) + tril(sigma)' - diag(diag(sigma));
            while rank(sigma) < nTask
                sigma = sigma + diag(diag(ones(nTask,nTask)*1e-8));
            end
            
            f11 = rand(1,n)*(minNorF - (mu(:,1)-3*sqrt(abs(sigma(1,1))))) + mu(:,1) - 3*sqrt(abs(sigma(1,1)));
            X   = [f11];
            tmp = 1;
            for t=2:nC+1
                gtt = rand(1,n)*(- (mu(:,t)-3*sqrt(abs(sigma(t,t))))) + mu(:,t) - 3*sqrt(abs(sigma(t,t)));
                X   = [X; gtt];
                tmp = tmp*(-(mu(:,t)-3*sqrt(abs(sigma(t,t)))));
            end
            EIvalue(i,:) = (minNorF-(mu(:,1)-3*sqrt(abs(sigma(1,1)))))*tmp*sum((minNorF-f11).*MonteCarloFun(X(:,:),mu,sigma))/n;
        end
    else                   % infeasible situation
        [ymu, ys2] =  predict_mtgp_all_tasks(Model, data, x_test);
        for i = 1:size(ys2, 1)
            mu    = ymu(i, 1:end);
            sigma = reshape(ys2(i,:), nTask-1, nTask-1);
            sigma = tril(sigma) + tril(sigma)' - diag(diag(sigma));
            while rank(sigma) < nTask - 1
                sigma = sigma + diag(diag(ones(nTask-1,nTask-1)*1e-8));
            end
            if nC == 1     % one constraint  
                EIvalue(i,:) = integral(@(z) normcdf((z-mu(:,1))/sigma(1,1)), 0, minNorG) - minNorG.*normcdf((0-mu(:,1))/sigma(1,1));
            else           % correlated EI (MonteCarlo integration)
                zSample = rand(1, n)*minNorG;
                Z       = [];
                ind     = 1;
                TMP1    = minNorG;
                TMP2    = minNorG;
                for t=1:nC
                    gSample = rand(1,n)*(minNorG - (mu(:,t)-3*sigma(t,t))) + mu(:,t) - 3*sigma(t,t);
                    Z       = [Z; gSample];
                    ind     = ind&(gSample>=mu(:,t)-3*sigma(t,t))&(gSample<=zSample);
                    TMP1    = TMP1*(minNorG - (mu(:,t)-3*sigma(t,t)));
                    TMP2    = TMP2*(-(mu(:,t)-3*sigma(t,t)));
                end
                TMP1    = TMP1*sum(MonteCarloFun(Z(:,ind),mu,sigma))/n;
                TMP2    = TMP2*sum(MonteCarloFun(zeros(nC,1),mu,sigma))/n;
                EIvalue(i,:) = TMP1 - TMP2;  
            end
        end
    end
end


function f = MonteCarloFun(x, mu, sigma)
    f = zeros(1, size(x, 2));
    for i = 1:size(x, 2)
        f(:,i) = mvnpdf(x(:,i), mu', sigma);
    end
end