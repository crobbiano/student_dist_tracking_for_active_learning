function [kdiv] = kldiv(A, B)
    % Compute the Kullback-Leibler divergence between two probability distributions.
    %
    % Input
    % -----
    % A, B : array-like
    % Probability distributions of equal length that sum to 1
    
    kdiv = nansum(A.*log2(A./B));
end