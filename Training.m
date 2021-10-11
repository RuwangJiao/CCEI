function [logtheta_all, data] = Training(xtrain, ytrain, nTask)
    % Train the MTGP model
    %
    % INPUT:
    % - xtrain      : training input data
    % - ytrain      : training output data
    % - nTask       : number of tasks
    %
    % OUTPUT
    % - logtheta_all: learned hyper-parameters
    % - data        : cell data for learning and prediction
    %
    % Author        : Ruwang Jiao
    
    addpath(genpath('scripts'));
    addpath(genpath('gpml'));

    covfunc_x    = {'covSEard'};
    irank        = nTask;                      % rank for Kf (1, ... M). irank=M -> Full rank
    ntrain       = size(xtrain, 1);
    nx           = ones(ntrain*nTask, 1);      % observations on each task-input point
    ind_kx_t     = linspace(1, size(xtrain,1), size(xtrain,1))';
    ind_kx_train = [];
    for i=1:nTask
        ind_kx_train = [ind_kx_train; ind_kx_t];
    end
    ind_kf_train = ones(size(ytrain,1), 1);
    for i=1:nTask
        ind_kf_train((i-1)*size(xtrain,1)+1:i*size(xtrain,1),:) = ind_kf_train((i-1)*size(xtrain,1)+1:i*size(xtrain,1),:)*i;
    end

    %% Assigns cell data for learning and prediction
    data = {covfunc_x, xtrain, ytrain, nTask, irank, nx, ind_kf_train, ind_kx_train};

    %% Hyper-parameter learning
    [logtheta_all, deriv_range] = init_mtgp_default(xtrain, covfunc_x, nTask, irank);
    logtheta_all                = learn_mtgp(logtheta_all, deriv_range, data);
end