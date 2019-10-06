% Learning Erdös-Rényi Random Graphs via Edge Detecting Queries
% 
% Zihan Li 1, Matthias Fresacher 2, and Jonathan Scarlett 1
% 
%   1 - National University of Singapore
%   2 - University of Adelaide
% 
% Corresponding authors:
% Matthias Fresacher - matthias.fresacher@adelaide.edu.au
% Jonathan Scarlett - scarlett@comp.nus.edu.sg
% 
% Published at the 33rd Conference on Neural Information Processing Systems
% (NeurIPS 2019), Vancouver, Canada.
% 
% December 2019
% 
% For a digital copy of the paper visit: arxiv.org/abs/1905.03410
%
%%
% Reproduces Figure 2 Left
 
close all;
 
%% Set parameters

% Number of trials
trials = 2000;

% Number of nodes
n = 50;

% Mean number of edges
k = 5;

% Sets parameter v for p
v = 1;

% Sets number of tests
t_vec = 25:25:450;

% Calculates number of possible edges
edges = nchoosek(n,2);

% Sets probability of edge
q = k/edges;

% Sets sample probability
p = sqrt(v/k);
 
% Adds required functions to working directory
addpath('Functions');
 
% Loads seed for random values to ensure repeatable results
load('seed');
rng(seed);
 
% Pre allocates success probabilities of algorithms
ps_COMP = zeros(size(t_vec));
ps_DD = zeros(size(t_vec));
ps_LP = zeros(size(t_vec));
ps_SSS = zeros(size(t_vec));

% Starts timer
tic;
 
% For each number of tests
for t_index = 1:length(t_vec)
    
    %% Runs simulations
    
    % Number of tests
    t = t_vec(t_index);
    % Resets number of successes for trivial algorithm
    succ_triv = 0;
 
    % Resets number of successes for LP rounded down algorithm
    succ_LP_down = 0;
 
    % Resets number of successes for LP rounded up algorithm
    succ_LP_up = 0;
 
    % Resets number of successes for LP rounded algorithm
    succ_LP = 0;
 
    % Resets number of successes for ILP algorithm
    succ_ILP = 0;
 
    % Resets number of successes for DD algorithm
    succ_DD = 0;
 
    % For each trial
    for tr = 1:trials
 
        if(mod(tr,5) == 0)
            disp(['Trial: ', num2str(tr), ' of ', num2str(trials), ...
                ' (with t = ', num2str(t), ')         [Total progress: ', ...
                num2str((t_index - 1)*(trials) + tr), '/', ...
                num2str(length(t_vec)*trials), '  ', ...
                num2str(round((t_index - 1)*(trials) + tr)*100/...
                (length(t_vec)*trials), 4), '%]']);
            
            % Display elapsed time
             disp(['Elapsed time is ', datestr(toc/(24*60*60), 'HH:MM:SS'), ...
                ' (hh:mm:ss).']);
        end
 
        % Generate random graph
        G_unknown = ER_Graph_Gen(n, q);
 
        % Produces LP estimate of graph
        [G_est_LP, G_est_ILP] = LP_Algo(G_unknown, t, p, 1);
 
        % Produces trivial estimate of graph
        [G_est_triv, G_est_DD] = DD_Algo(G_unknown, t, p);
 
        % Produced rounded down estimate
        G_est_down = floor(G_est_LP);
 
        % Produced rounded up estimate
        G_est_up = ceil(G_est_LP);
 
        % Produced rounded estimate
        G_est_round = round(G_est_LP);
 
        % Tests if trivial estimate is correct
        if G_unknown == G_est_triv
 
            % If correct
            succ_triv = succ_triv + 1;
        end
 
        % Tests if LP rounded down estimate is correct
        if G_unknown == G_est_down
 
            % If correct
            succ_LP_down = succ_LP_down + 1;
        end
 
        % Tests if LP rounded up estimate is correct
        if G_unknown == G_est_up
 
            % If correct
            succ_LP_up = succ_LP_up + 1;
        end
 
        % Tests if LP rounded estimate is correct
        if G_unknown == G_est_round
 
            % If correct
            succ_LP = succ_LP + 1;
        end
 
        % Tests if ILP estimate is correct
        if G_unknown == G_est_ILP
 
            % If correct
            succ_ILP = succ_ILP + 1;
        end
 
        % Tests if DD estimate is correct
        if G_unknown == G_est_DD
 
            % If correct
            succ_DD = succ_DD + 1;
        end
    end
    
    % Calculates success probability
    ps_COMP(t_index) = succ_triv / trials;
    ps_DD(t_index) = succ_DD / trials;
    ps_LP(t_index) = succ_LP / trials;
    ps_SSS(t_index) = succ_ILP / trials;
end
 
%% Display Summary
 
disp(' ');
disp(['Number of nodes:        ', num2str(n)]);
disp(['Mean number of edges:   ', num2str(k)]);
disp(['v in setting of p:      ', num2str(v)]);
disp(' ');
disp(['Total elapsed time is ', datestr(toc/(24*60*60), 'HH:MM:SS'), ...
        ' (hh:mm:ss).']);
 
%% Saves Results
 
save('Data_Fig_2_Left'); 
 
%% Creates Figure
 
figure(2);
hold on;
grid on;
xlim([25, 450]);
xticks(50:50:450);
yticks(0:0.2:1);
plot(t_vec, ps_SSS, 'ko-', 'LineWidth' , 1);
plot(t_vec, ps_LP, 'r', 'LineWidth' , 2);
plot(t_vec, ps_DD, 'b-.', 'LineWidth' , 2);
plot(t_vec, ps_COMP, 'g--', 'LineWidth' , 2);
xlabel('Number of tests');
ylabel('Success probability');
legend('SSS', 'LP', 'DD', 'COMP', 'Location', 'SouthEast');