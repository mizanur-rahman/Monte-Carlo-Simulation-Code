

x = randn(50,1)*1+50;                      % Create Data
% SEM = std(x)/sqrt(length(x));               % Standard Error
% ts = tinv([0.95],length(x)-1);      % T-Score
% CI = mean(x) + ts*SEM; 
ci = 0.95;
alpha = 1 - ci;
ybar=mean(x);

n = length(x); %number of elements in the data vector
T_multiplier = tinv(.975, n-1);
 
ci95 = T_multiplier*std(x)/sqrt(n);

Lowest_value_of_95percent_confidence_interval= ybar - ci95
Highest_value_of_the_95percent_cofidence_interval= ybar + ci95