function Adj_Mat = ER_Graph_Gen(n, q)
%ER_Graph_Gen Generates an Erdos-Renyi Random Graph
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
%   Generates an Erdos-Renyi random graph with n nodes where the
%   probability an edge is present is q.
%
%   Inputs:
%       n - integer, number of nodes
%       q - probability, 0 < q <1, probability of an edge being present
%
%   Outputs:
%       Adj_Mat - n x n adjacency matrix of graph
 
% Number of edges
ed = nchoosek(n, 2);
 
% Random numbers
ran = rand(1, ed);
 
% Generate zeros
z = zeros(1, ed);
 
% Generate ones
o = ones(1, ed);
 
% Pre-allocate binary random numbers
ran_bin = zeros(1, ed);
 
% Convert to binary numbers
ran_bin(ran >= q) = z(ran >= q);
ran_bin(ran < q) = o(ran < q);
 
% Generate adjacency matrix
Adj_Mat = triu(ones(n),1);
 
% Populate upper triangle of matrix with edges
Adj_Mat(Adj_Mat == 1) = ran_bin;
 
end