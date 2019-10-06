function [Est_LP, Est_ILP] = LP_Algo(Adj_Mat, t, p, ilp)
%LP_Algo Estimates an unknown graph using the LP algorithm
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
%   Estimates an unknown graph using the LP algorithm using t tests
%   each containing l nodes per sample
%
%   Inputs:
%       Adj_Mat - adjacency matrix of unknown graph
%       t   - integer, number of tests
%       p   - probability of a given node being included in a test
%       ilp - logical, 0 - no ILP output, 1 - include ILP output
%
%   Outputs:
%       Est_LP  - adjacency matrix of estimated graph using LP
%       Est_ILP - adjacency matrix of estimated graph using ILP
 
%% Sets parameters
 
% Number of nodes
n = length(Adj_Mat);
 
% Calculates number of possible edges
edges = nchoosek(n,2);
 
% Creates sample matrix
sample = zeros(t, n);
 
% Creates test outcome vector
y = zeros(t, 1);
 
% Creates constraints type vector
const_type = char(t, 1);
 
% Pre allocates memory for vector of edges included in test
edge_matr = zeros(t, edges);
 
%% Determines samples and test outcomes
 
% For each test
for test = 1:t
    
    % Select random sample
    sample(test,:) = binornd(1,p,[1,n]);
    indices = find(sample(test,:) == 1);
    
    % Calculates number of edges in sample
    edge_num = sum(sum(Adj_Mat(indices,indices)));
    
    % If no edges are present
    if edge_num == 0
        
        % Store 0 test outcome
        y(test) = 0;
        
        % Sets constraints type to equal
        const_type(test) = '=';
        
    % If edges are present
    else
        
        % Store 1 test outcome
        y(test) = 1;
        
        % Sets constraints type to more than or equal
        const_type(test) = '>';
        
    end
    
    % Create vector for edges in test
    edge_vec = zeros(1, edges);
    
    % Sets index for edge vector
    edge_vec_index = 1;
    
    % Create vector of edges from node i
    node_edge_vec = zeros(1,n);
            
    % Set entry 1 one for edges between nodes in sample
    node_edge_vec(indices) = ones(1, length(indices));
    
    % For each node except last one
    for node_i = 1:n-1
        
        % Check if node is  in sample
        if isempty(find(indices == node_i, 1)) == 0

            % Add node to edge vector
            edge_vec(edge_vec_index:edge_vec_index + n - node_i - 1) = ...
                node_edge_vec(node_i + 1:end);
 
        end
        
        % Update edge vector index
        edge_vec_index = edge_vec_index + n - node_i;
            
        
    end
    
    % Add edges vector to matrix
    edge_matr(test, :) = edge_vec;
    
end
 
%% Solves LP
 
% Sets constraint matrix
model.A = sparse(edge_matr);
 
% Sets constraints
model.rhs = y;
 
% Sets constraints types
model.sense = const_type;
 
% Sets the objective function
model.obj = ones(1, edges);
 
% Sets type of optimisation
model.modelsense = 'Min';
 
% Sets lower bound
model.lb = zeros(1, edges);
 
% Sets upper bound
model.ub = ones(1, edges);
 
% Sets variables as continuous
model.vtype = 'C';

% Stops working output
params.OutputFlag = 0;

% Solves LP
result_LP = gurobi(model, params);
 
% Saves estimates of edges
edges_est_LP = result_LP.x;

% If ILP output is desired
if ilp == 1
    % Sets variables as binary
    model.vtype = 'B';

    % Solves ILP
    result_ILP = gurobi(model, params);

    % Saves estimates of edges
    edges_est_ILP = result_ILP.x;
    
else
    
    % Create dummy output
    edges_est_ILP = zeros(size(edges_est_LP));
    
end
 
%% Formats results
 
% Generate estimated adjacency matrices
Est_LP = triu(ones(n), 1);
Est_LP = Est_LP';
Est_ILP = Est_LP;
 
% Populate upper triangle of matrix with edges for LP
Est_LP(Est_LP == 1) = edges_est_LP;
Est_LP = Est_LP';
 
% Populate upper triangle of matrix with edges for ILP
Est_ILP(Est_ILP == 1) = edges_est_ILP;
Est_ILP = Est_ILP';
 
end