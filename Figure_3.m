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
% Reproduces Figure 3
 
close all;
 
%% Set parameters

% Number of trials
trials = 1000;

% Number of nodes
n_vec = 80:20:140;

% Sets parameter v for p
v = 1;

% Sets number of tests
t_vec = 25:25:800;

% Adds required functions to working directory
addpath('Functions');

% Loads seed for random values to ensure repeatable results
load('seed');
rng(seed);

% Pre allocates success probabilities of algorithms
    % Each row corresponds to a given value of n
    % Each colomn corresponds to a given value of t
ps_COMP = zeros(length(n_vec), length(t_vec));
ps_DD = zeros(length(n_vec), length(t_vec));
ps_LP = zeros(length(n_vec), length(t_vec));

% Pre allocates normalised tests
    % Each row corresponds to a given value of n
    % Each colomn corresponds to a given value of t
t_norm = zeros(length(n_vec), length(t_vec));

% Starts timer
tic;

% For each number of nodes
for n_index = 1:length(n_vec)
    
    % Number of nodes
    n = n_vec(n_index);

    % Mean number of edges
    k = n/10;

    % Calculates number of possible edges
    edges = nchoosek(n,2);

    % Sets probability of edge
    q = k/edges;

    % Sets sample probability
    p = sqrt(v/k);

    % For each number of tests
    for t_index = 1:length(t_vec)

        %% Runs simulations

        % Number of tests
        t = t_vec(t_index);
        
        % Number of normalised tests
        t_norm(n_index, t_index) = t/(k*log(1/q));
        
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
                disp(['Nodes: ', num2str(n), ', Trial: ', num2str(tr), ...
                    ' of ', num2str(trials), ' (with t = ', num2str(t), ...
                    ')         [Total progress: ', ...
                    num2str(length(t_vec)*trials*(n_index - 1) + ...
                    (t_index - 1)*(trials) + tr), '/', ...
                    num2str(length(t_vec)*trials*length(n_vec)), '  ', ...
                    num2str(round(length(t_vec)*trials*(n_index - 1) + ...
                    (t_index - 1)*trials + tr)*100/...
                    (length(t_vec)*trials*length(n_vec)), 4), '%]']);

                % Display elapsed time
                 disp(['Elapsed time is ', datestr(toc/(24*60*60), ...
                     'HH:MM:SS'), ' (hh:mm:ss).']);
            end

            % Generate random graph
            G_unknown = ER_Graph_Gen(n, q);

            % Produces LP estimate of graph
            [G_est_LP, ~] = LP_Algo(G_unknown, t, p, 0);

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

            % Tests if DD estimate is correct
            if G_unknown == G_est_DD

                % If correct
                succ_DD = succ_DD + 1;
            end
        end

        % Calculates success probability
        ps_COMP(n_index, t_index) = succ_triv / trials;
        ps_DD(n_index,t_index) = succ_DD / trials;
        ps_LP(n_index,t_index) = succ_LP / trials;
    end

    %% Display Summary

    disp(' ');
    disp(['Number of nodes:        ', num2str(n)]);
    disp(['Mean number of edges:   ', num2str(k)]);
    disp(['v in setting of p:      ', num2str(v)]);
    disp(' ');
    disp(['Total elapsed time is ', datestr(toc/(24*60*60), 'HH:MM:SS'), ...
            ' (hh:mm:ss).']);
        
end
 
%% Saves Results
 
save('Data_Fig_3'); 
 
%% Creates Figures

% Based on example code by Cris LaPierre in order to have two legends.
% Source code available at:
% https://au.mathworks.com/matlabcentral/answers/430791-how-to-add-a-second-legend-box-to-a-figure-without-new-plots
 
f3 = figure(3);
firstax3 = axes(f3); 
hold on;
L1_3 = plot(t_vec, ps_LP(1,:), 'r-', 'LineWidth' , 2, 'Parent', firstax3);
L2_3 = plot(t_vec, ps_DD(1,:), 'b-', 'LineWidth' , 2, 'Parent', firstax3);
L3_3 = plot(t_vec, ps_COMP(1,:), 'g-', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_LP(2,:), 'r--', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_LP(3,:), 'r-.', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_LP(4,:), 'r:', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_DD(2,:), 'b--', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_DD(3,:), 'b-.', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_DD(4,:), 'b:', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_COMP(2,:), 'g--', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_COMP(3,:), 'g-.', 'LineWidth' , 2, 'Parent', firstax3);
plot(t_vec, ps_COMP(4,:), 'g:', 'LineWidth' , 2, 'Parent', firstax3);
set(firstax3, 'Box', 'off','XLim',[25 800], 'XTick', 50:50:800, ...
    'YTick', 0:0.2:1);
leg3 = legend(firstax3, {'LP', 'DD', 'COMP'}, 'Location', 'east');
pos3 = [0, -0.22, 0, 0];
set(leg3, 'Position', leg3.Position + pos3);
secondax3 = copyobj(firstax3, gcf);
delete(get(secondax3, 'Children'))

% Hidden dummy plots to create second legend
H1_3 = plot([-1 0], [0 0], 'k-', 'LineWidth' , 2, 'Parent', secondax3);
H2_3 = plot([-1 0], [0 0], 'k--', 'LineWidth' , 2, 'Parent', secondax3);
H3_3 = plot([-1 0], [0 0], 'k-.', 'LineWidth' , 2, 'Parent', secondax3);
H4_3 = plot([-1 0], [0 0], 'k:', 'LineWidth' , 2, 'Parent', secondax3);
set(secondax3, 'Color', 'none', 'XTick', []) 
legend ([H1_3 H2_3 H3_3 H4_3], {'$n = 80$', '$n = 100$', '$n = 120$', ...
    '$n = 140$'}, 'Location', 'southeast', 'Color', 'white', ...
    'Interpreter', 'latex');
grid on;
hold off;
xlabel('Number of tests');
ylabel('Success probability');

f4 = figure(4);
firstax4 = axes(f4); 
hold on;
L1_4 = plot(t_norm(1,:), ps_LP(1,:), 'r-', 'LineWidth' , 2, ...
    'Parent', firstax4);
L2_4 = plot(t_norm(1,:), ps_DD(1,:), 'b-', 'LineWidth' , 2, ...
    'Parent', firstax4);
L3_4 = plot(t_norm(1,:), ps_COMP(1,:), 'g-', 'LineWidth' , 2, ...
    'Parent', firstax4);
plot(t_norm(2,:), ps_LP(2,:), 'r--', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(3,:), ps_LP(3,:), 'r-.', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(4,:), ps_LP(4,:), 'r:', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(2,:), ps_DD(2,:), 'b--', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(3,:), ps_DD(3,:), 'b-.', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(4,:), ps_DD(4,:), 'b:', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(2,:), ps_COMP(2,:), 'g--', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(3,:), ps_COMP(3,:), 'g-.', 'LineWidth' , 2, 'Parent', firstax4);
plot(t_norm(4,:), ps_COMP(4,:), 'g:', 'LineWidth' , 2, 'Parent', firstax4);
set(firstax4, 'Box', 'off','XLim',[0 16], 'XTick', 0:2:16, ...
    'YTick', 0:0.2:1);
leg4 = legend(firstax4, {'LP', 'DD', 'COMP'}, 'Location', 'east');
pos4 = [0, -0.22, 0, 0];
set(leg4, 'Position', leg4.Position + pos4);
secondax4 = copyobj(firstax4, gcf);
delete(get(secondax4, 'Children'))

% Hidden dummy plots to create second legend
H1_4 = plot([-1 0], [0 0], 'k-', 'LineWidth' , 2, 'Parent', secondax4);
H2_4 = plot([-1 0], [0 0], 'k--', 'LineWidth' , 2, 'Parent', secondax4);
H3_4 = plot([-1 0], [0 0], 'k-.', 'LineWidth' , 2, 'Parent', secondax4);
H4_4 = plot([-1 0], [0 0], 'k:', 'LineWidth' , 2, 'Parent', secondax4);
set(secondax4, 'Color', 'none', 'XTick', []) 
legend ([H1_4 H2_4 H3_4 H4_4], {'$n = 80$', '$n = 100$', '$n = 120$', ...
    '$n = 140$'}, 'Location', 'southeast', 'Color', 'white', ...
    'Interpreter', 'latex');
grid on;
hold off;
xlabel('Normalized number of tests');
ylabel('Success probability');