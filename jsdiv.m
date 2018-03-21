function [div] = jsdiv(P, Q)
    % Compute the Jensen-Shannon divergence between two probability distributions.
    %
    % Input
    % -----
    % P, Q : array-like
    % Probability distributions of equal length that sum to 1
    
    M = 0.5 * (P + Q);
    div = 0.5 * (kldiv(P, M) + kldiv(Q, M));
end