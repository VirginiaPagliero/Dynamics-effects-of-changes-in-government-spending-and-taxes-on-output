%% Intro
clear 
clc
close all
clear global 

%% Global variables
global p             % # of VAR lags
global T             % # of observations  
global Sigma_u       % error covariance matrix
global ParamNumberA  % structural parameters of A to estimate 
global ParamNumberB  % structural parameters of B to estimate 
global ParamNumber

%% Importing dataset 
importdata CK_RESTUD_DATASET.xlsx ;   

VariablesSelected = [3 4 2];  % selects TAX, G and GDP
StartSample = 1; 
EndSample = 119;

data = ans.data.FINAL(StartSample:EndSample,VariablesSelected);
dates1 = ans.data.FINAL(StartSample:EndSample, 1); %date

%% Define variables
TAX = data(:,1);
G = data(:,2);
GDP = data(:,3);

p=4;
T=size(data,1)-p;
M=size(data,2);
DuplicationMatrix = DuplicationMatrixFunction(M);
mDD=(DuplicationMatrix'*DuplicationMatrix)^(-1)*DuplicationMatrix'; 
mNN=DuplicationMatrix*mDD;

%% ANALYSIS VARIABLES IN LEVELS 
exp_tax = exp(TAX);
exp_gdp = exp(GDP);
exp_g = exp(G);

%% Plot
figure
plot(exp_tax,'LineWidth',2);
hold on
plot(exp_g,'LineWidth',2);
hold on
exp_series = plot(exp_gdp,'LineWidth',2);
legend('TAX','G','GDP')
xlabel('Quarters')
title '{\bf Levels of Tax, Spending and GDP }';
axis tight
grid on
saveas(exp_series,'exp_series','png');

%% Ratio
ratio_tax = exp_tax./exp_gdp;
ratio_g = exp_g./exp_gdp;

%% Plot ratio
figure
plot(ratio_tax,'LineWidth',2);
hold on
ratio_levels = plot(ratio_g,'LineWidth',2);
legend('TAX/GDP','G/GDP')
xlabel('Quarters')
title '{\bf Levels of Tax/GDP and G/GDP}';
axis tight
grid on
saveas(ratio_levels,'ratio_levels','png');

%% ANALYSIS VARIABLES IN LOGS
plot(TAX,'LineWidth',2);
hold on
plot(G,'LineWidth',2);
hold on
series = plot(GDP,'LineWidth',2);
legend('TAX','G','GDP')
xlabel('Quarters')
title '{\bf Tax, Spending and GDP}';
axis tight
grid on
saveas(series,'series','png'); 

%% Tax on GDP
figure
tax_gdp = plot(TAX-GDP,'LineWidth',2);
xlabel 'Quarters';
title '{\bf Ratio of tax on GDP}';
axis tight
grid on
saveas(tax_gdp,'tax_gdp','png'); 

%% Spending on GDP
figure
g_gdp = plot(G-GDP,'LineWidth',2);
xlabel 'Quarters';
title '{\bf Ratio of spending on GDP}';
axis tight
grid on
saveas(g_gdp,'g_gdp','png');
