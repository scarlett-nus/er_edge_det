function [Est_triv, Est_DD] = DD_Algo(Adj_Mat, t, p)
%DD_Algo Estimates an unknown graph using the trivial and DD algorithm
%
%   Learning Erdös-Rényi Random Graphs via Edge Detecting Queries
% 
%   Zihan Li 1, Matthias Fresacher 2, and Jonathan Scarlett 1
% 
%   	1 - National University of Singapore
%   	2 - University of Adelaide
% 
%   Corresponding authors:
%   Matthias Fresacher - matthias.fresacher@adelaide.edu.au
%   Jonathan Scarlett - scarlett@comp.nus.edu.sg
% 
%   Published at the 33rd Conference on Neural Information Processing
%   Systems (NeurIPS 2019), Vancouver, Canada.
% 
%   December 2019
% 
%   For a digital copy of the paper visit: arxiv.org/abs/1905.03410
%%
%   Estimates an unknown graph using the trivial and DD algorithm using t 
%   tests each containing l nodes per sample
%
%   Inputs:
%       Adj_Mat - adjacency matrix or unknown graph
%       t - integer, number of tests
%       p - probability of a given node being included in a test
%
%   Outputs:
%       Est_triv - adjacency matrix of estimated trivial graph
%       Est_DD   - adjacency matrix of estimated DD graph

% Number of nodes
n = length(Adj_Mat);
 
% Generate estimated graph with all edges present
Est_triv = triu(ones(length(Adj_Mat)),1);
Est_DD = zeros(length(Adj_Mat));

% Pre allocates memory for samples
sample = zeros(t,n);
 
% For each test
for test = 1:t
    
    % Select random sample
    sample(test,:) = binornd(1,p,[1,n]);
    indices = find(sample(test,:) == 1);
    
    % Calculates number of edges in sample
    edge_num = sum(sum(Adj_Mat(indices,indices)));
    
    % If no edges are present
    if edge_num == 0
        
        % Remove edges from estimate
        Est_triv(indices, indices) = 0;
    end
    
end

% For each test
for test = 1:t
    
    % Selects sample
    indices = find(sample(test,:) == 1);
    
    % Calculates number of edges in sample in trivial estimate
    edge_num = sum(sum(Est_triv(indices, indices)));
    
    % If exactly one edge is present
    if edge_num == 1
        
        % Insert edges to estimate
        Est_DD(indices, indices) = Est_triv(indices, indices);
    end
    
end
 
end