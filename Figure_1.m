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
% Reproduces Figure 1
 
close all;
 
% Sets values for theta
thetas = linspace(0, 1, 1000);
 
% Pre allocates variables
Rs_COMP = zeros(size(thetas));
Rs_DD = zeros(size(thetas));
Rs_CONV = ones(size(thetas));
Rs_SSS = zeros(size(thetas));
 
% For each theta
for i=1:length(thetas)
    
    % Sets value of theta
    theta = thetas(i);
    
    % Calculates values of parameters
    const = 2*(1-theta) / log(2);
    Rs_COMP(i) = const / (2 * exp(1));
    
    if(theta > 3/4)
        
        val = 2 - theta;
        
    elseif(theta > 1/2)
        
        val = 3*theta - 1;
        
    else
        
        val = theta;
        
    end
    
    Rs_DD(i) = const / max( max(2-2*theta,val)*exp(1), 2*theta*exp(1) );
    Rs_SSS(i) = min(1,const / (2*theta*exp(1)));
end
 
% Creates figure
figure(1)
hold on;
grid on;
xticks(0:0.2:1);
ylim([0, 1.06]);
yticks(0:0.2:1.1);
plot(thetas,Rs_CONV,'k', 'LineWidth' , 2);
plot(thetas,Rs_COMP,'g--', 'LineWidth' , 2);
plot(thetas,Rs_DD,'b-.', 'LineWidth' , 2);
plot(thetas,Rs_SSS,'r:', 'LineWidth' , 2);
xlabel('Sparsity Parameter \theta');
ylabel('Relative Inverse #Tests');
text(0.77, 0.96, 'Converse');
text(0.46, 0.7, 'SSS');
text(0.29, 0.56, 'DD');
text(0.275, 0.31, 'COMP');
legend('off');