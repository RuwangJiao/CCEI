function [ Ypred, Vpred ] = predict_mtgp_all_tasks(logtheta_all, data, xtest)
    %PREDICT_MTGP_ALL_TASKS Makes predictions at all points xtest for all tasks
    %
    % INPUT:
    % - logtheta_all : all hyperparameters
    % - data         : cell data in the order 
    %                  [covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train, ind_kx_train]
    % - xtest        : Test points 
    %
    % OUTPUT
    % - YPred        : (Ntest x M) Matrix of Mean MTGP Predictions 
    % - Vpred        : (Ntest x M) Matrix of MTGP Variances
    %                   Where M is the number of tasks and Ntest: number of test points            
    % Author         : Ruwang Jiao

    [covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train, ind_kx_train] = deal(data{:});
    Ntest = size(xtest, 1);
    [alpha, Kf, L, K, Kxstar, Kss] = alpha_mtgp(logtheta_all, covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train,ind_kx_train, xtest);
    all_Kxstar = Kxstar(ind_kx_train,:);

    Ypred = zeros(Ntest,M);
    Vpred = zeros(Ntest,M*M);
    for k = 1:Ntest
       Kstar      = kron(Kf, all_Kxstar(1:size(xtrain,1),k));
       Ypred(k,:) = (Kstar'*alpha)';
       tmp        = Kf - Kstar'*solve_chol(L',Kstar);
       tmp        = tril(tmp) + tril(tmp)' - diag(diag(tmp));
       tmp        = tmp + diag(diag(ones(M,M)*1e-8));
       while rank(tmp) < M || prod(diag(tmp))<=0  || det(tmp)<=0
            tmp   = tmp + diag(diag(ones(M,M)*1e-8));            % Regularization to avoid program crashes
       end
       Vpred(k,:) = tmp(:)';
    end
end