% This function is used to calculate the Spearman‘s rank correlation
% coefficient among the objective and constraint functions

problemSet = [1:9];
for problemIndex = 9
    prob         = problemSet(problemIndex)
    [nO, nC, nD, lu] = problem(prob);
    sampleSize   = 1000000;
    rand('seed', sum(100*clock));
    x_train      = lhsdesign(sampleSize, nD, 'criterion','maximin', 'iteration',100);
    P            = repmat(lu(1,:),iniSize,1) + x_train.*repmat((lu(2,:) - lu(1,:)), iniSize, 1);
    [objF, conV] = fitness(P, prob);
    r01          = corr(objF, conV(:,1), 'type', 'Spearman')
%     r02           = corr(objF, conV(:,2), 'type', 'Spearman')
%     r03           = corr(objF, conV(:,3), 'type', 'Spearman')
%     r04           = corr(objF, conV(:,4), 'type', 'Spearman')
%     r12           = corr(conV(:,1), conV(:,2), 'type', 'Spearman')
%     r13           = corr(conV(:,1), conV(:,3), 'type', 'Spearman')
%     r14           = corr(conV(:,1), conV(:,4), 'type', 'Spearman')
%     r23           = corr(conV(:,2), conV(:,3), 'type', 'Spearman')
%     r24           = corr(conV(:,2), conV(:,4), 'type', 'Spearman')
%     r34           = corr(conV(:,3), conV(:,4), 'type', 'Spearman')
end