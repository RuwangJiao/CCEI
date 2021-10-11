function [ logtheta_all, deriv_range] = init_mtgp_default(xtrain, covfunc_x, M, irank  )
    % Initializes parameters of mtgp by default
    % You should pay careful attention to modifying this as it may be a 
    % bad initialization for your problem!
    % 
    % INPUT:
    % - xtrain: Input training data
    % - covfunc_x: Input covariance function
    % - M: Number of tasks
    % - irank: Rank required for Kf
    %
    % OUTPUT:
    % - logtheta_all: Vector of all hyper-paramters
    % - deriv_range: Indices of hyper-parameters to optimize for
    %
    % Edwin V. BOnilla

    nlf           = irank*(2*M - irank +1)/2;    % Number of parameters for Lf
    theta_lf0     = init_Kf(M, irank);
    theta_kx0     = init_kx(xtrain, covfunc_x);
    theta_sigma0  = init_sigma(M);
    logtheta_all  = [theta_lf0; theta_kx0; theta_sigma0];

    % assumes that we don't want to optimize the signal variance of kx
    deriv_range   = 1 : length(logtheta_all);
    deriv_range(nlf+length(theta_kx0)) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_kx0 = init_kx(xtrain, covfunc_x)
    % Init the length scales {li}1 ≤ i ≤ D and the signal variance sigma_f
    D         = size(xtrain,2);
    L         = eval(feval(covfunc_x{:}));
    theta_kx0 = log(ones(L,1)); 
    %theta_kx0 = log(ones(L,1)./2); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_lf0 = init_Kf(M,irank)
    % Init to diagonal matrix (No task correlations)
    Kf0 = eye(M);                
    Lf0 = chol(Kf0)';
    theta_lf0 = lowtri2vec_inchol(Lf0,M,irank);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_sigma0 = init_sigma(M) 
    % Init the noise variances sigma
    theta_sigma0 =  (1e-7)*rand(M,1);  
    %theta_sigma0 =  (0.0001)*rand(M,1); 
end



