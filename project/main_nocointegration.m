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

%% Ratio
ratio_tax = exp_tax./exp_gdp;
ratio_g = exp_g./exp_gdp;

%% VAR - NO COINTEGRATION
VAR_Const = [NaN; NaN; NaN];
VAR_trend = [NaN; NaN; NaN];

Pi1 =[NaN NaN NaN;     % Pi matrices collecting parameters associated with lags
       NaN NaN NaN; 
       NaN NaN NaN];

Pi2 =[NaN NaN NaN; 
       NaN NaN NaN; 
       NaN NaN NaN];      

Pi3 =[NaN NaN NaN; 
       NaN NaN NaN; 
       NaN NaN NaN];

Pi4 =[NaN NaN NaN; 
       NaN NaN NaN; 
       NaN NaN NaN];

VAR_Pi = {Pi1 Pi2 Pi3 Pi4};

VAR= varm('Constant',VAR_Const,'AR',VAR_Pi, 'Trend', VAR_trend); % Estimate VAR
[EstVAR,EstSE,logLikVAR,Residuals] = estimate(VAR,data);

const = EstVAR.Constant;
trend = EstVAR.Trend;
mP1 = EstVAR.AR{1,1};
mP2 = EstVAR.AR{1,2};
mP3 = EstVAR.AR{1,3};
mP4 = EstVAR.AR{1,4};

Sigma_u = (Residuals'*Residuals)/T;  % estimated covariance matrix
Sigma_u_Sample = Sigma_u;

%% Parameters
Sigma_u_sample = Sigma_u;

ParamNumberA =[3 ...     
              6 ...     
              ]';

ParamNumberB =[1 ...
               4 5 ...
               9]';

ParamNumber = [ParamNumberA; ParamNumberB];

Matrix_SelectedA = zeros(size(Sigma_u_sample,1),size(Sigma_u_sample,1));
Matrix_SelectedB = zeros(size(Sigma_u_sample,1),size(Sigma_u_sample,1));
 
for c_parA = 1 : size(ParamNumberA,1)
    Matrix_SelectedA(ParamNumberA(c_parA,1)) = 1;    % this loop puts parameters in matrix form
end

for c_parB = 1 : size(ParamNumberB,1)
    Matrix_SelectedB(ParamNumberB(c_parB,1)) = 1;    % this loop puts parameters in matrix form
end

StructuralParamA = size(ParamNumberA,1);        
StructuralParamB = size(ParamNumberB,1);
StructuralParam = StructuralParamA + StructuralParamB; % dimension of vector of structural parameters

%% Maximum Likelihood
InitialValue = (randn(StructuralParam,1)/10);

options = optimset('MaxFunEvals', 200000,'TolFun',1e-500,'MaxIter', 200000,'TolX',1e-50);   
[StructuralParam_Estimation_MATRIX,Likelihood_SVAR,exitflag,output,grad,Hessian_MATRIX] = fminunc('Likelihood_SVAR2',InitialValue',options);
SE_Hessian_MATRIX = diag(inv(Hessian_MATRIX)).^0.5;  % computes Hessian-based standard errors

%% Estimated parameters
A = eye(size(Sigma_u,1),size(Sigma_u,1));
A(1,3) = -2.08;
B = zeros(size(Sigma_u,1),size(Sigma_u,1));

SE_A = zeros(size(Sigma_u,1),size(Sigma_u,1));
SE_B = zeros(size(Sigma_u,1),size(Sigma_u,1));

HSelectionA = zeros(M*M,StructuralParamA);       % matrix S_A
HSelectionB = zeros(M*M,StructuralParamB);       % matrix S_B

for a_par = 1 : size(ParamNumberA,1);
    A(ParamNumberA(a_par,1)) = StructuralParam_Estimation_MATRIX(a_par); % puts the estimated elements of A in the right place       
    SE_A(ParamNumberA(a_par,1))= SE_Hessian_MATRIX(a_par);
    HSelectionA(ParamNumberA(a_par,1),a_par) = 1;                  % puts "1" in the correct place of selection matrix S_A 
end

i = StructuralParam_Estimation_MATRIX(size(ParamNumberA,1)+1:StructuralParam)
for b_par = 1:size(ParamNumberB,1)
    B(ParamNumberB(b_par,1)) = i(b_par);
    SE_B(ParamNumberB(b_par,1))= SE_Hessian_MATRIX(b_par);
    HSelectionB(ParamNumberB(b_par,1),b_par) = 1;
end        

%% Sign normalization
if B(1,1)<0
   B(:,1)=-B(:,1);
end

if B(2,2)<0
   B(:,2)=-B(:,2); 
end

if B(3,3)<0
   B(:,3)=-B(:,3);
end    

%% IRFs
HorizonIRF = 20;

C_IRF = A^(-1)*B;                                   % instantaneous impact at h=0
J=[eye(M) zeros(M,M*(p-1))];                        % selection matrix J used in IRF computation 
CompanionMatrix = [mP1 mP2 mP3 mP4;                 % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA(:,:,h+1)=J*CompanionMatrix^h*J'*C_IRF;
    end

for i = 0 : HorizonIRF       
    tax_tax(1,i+1) = TETA(1,1,i+1);
    spend_tax(1,i+1) = TETA(1,2,i+1);
    gdp_tax(1,i+1) = TETA(1,3,i+1);  
    tax_spend(1,i+1) = TETA(2,1,i+1);
    spend_spend(1,i+1) = TETA(2,2,i+1);
    gdp_spend(1,i+1) = TETA(2,3,i+1);
    tax_gdp(1,i+1) = TETA(3,1,i+1);
    spend_gdp(1,i+1) = TETA(3,2,i+1);
    gdp_gdp(1,i+1) = TETA(3,3,i+1);
end

%% Dynamic Multipliers
ratioGDPTAX = exp_gdp./exp_tax;
mratioGDPTAX = mean(ratioGDPTAX);

ratioGDPG = exp_gdp./exp_g;
mratioGDPG = mean(ratioGDPG);

multTAX_noCoint = (tax_gdp./tax_tax(1,1))*mratioGDPTAX; 
multG_noCoint = (spend_gdp./spend_spend(1,1))*mratioGDPG; 

%% Bootstrap
mean_residuals = mean(Residuals);

Recentered_residuals = Residuals - mean_residuals;

Z = 999;

resampledResidualsArray= cell(Z,1);
bootstrappedDataArray = cell(Z,1);

% Matrices to store different estimates in the 999 BS samples
Trend_starArray = cell(Z,1);
Const_starArray = cell (Z,1);
mP1_starArray = cell (Z,1);
mP2_starArray = cell (Z,1);
mP3_starArray = cell (Z,1);
mP4_starArray = cell (Z,1);
Sigma_u_starArray = cell (Z,1);
C_starArray = cell(Z,1);
SE_C_starArray = cell(Z,1);
necessary_and_sufficent_condition_starArray = cell(Z,1);
missing_condition_starArray = cell(Z,1);
TETA_starArray = cell(Z,1);
taxongdp_starArray = cell (Z,1);
spendingongdp_starArray = cell (Z,1);
DYNAMICMULTIPLIERTAX_starArray = cell (Z,1);
DYNAMICMULTIPLIERSPENDING_starArray = cell(Z,1);

%%
rng(126)

InitialValue = (randn(StructuralParam,1)/10);
% Loop through each BS sample
for j=1:Z
   % IID resampled residuals matrix
   % Initialize the matrix with the first value of the recentered residuals vector
   resampledResiduals = Recentered_residuals(1, :);
    % Loop through each time point
   for t = 2:T
      % Draw from recentered residuals with probability 1/T
      n = floor((T-1)*rand(1,1))+1;

      % Assign the residuals at time n to the current time point in the resampled residuals matrix
      resampledResiduals(t, :) = Recentered_residuals(n, :);
   end

   % Store the resampled residuals matrix in the array
   resampledResidualsArray{j} = resampledResiduals;

   bootstrappedData = zeros(T, M);

%Fix initial conditions for BS dataset
bootstrappedData(1:p, :) = data(1:p, :);

% Loop through remaining time points
for t = p+1:T
   % Generate the bootstrapped data using the VAR model
   bootstrappedData(t, :) = const + trend*t + mP1*(bootstrappedData(t-1, :)') + mP2*(bootstrappedData(t-2, :)') + mP3*(bootstrappedData(t-3, :)') + mP4*(bootstrappedData(t-4, :)') + resampledResiduals(t, :)';
end
bootstrappedDataArray{j} = bootstrappedData;


[EstVAR,EstSE,logLikVAR,Residuals] = estimate(VAR,bootstrappedData);

mPi1c = EstVAR.AR{1,1};
mPi2c = EstVAR.AR{1,2};
mPi3c = EstVAR.AR{1,3};
mPi4c = EstVAR.AR{1,4};
mconst = EstVAR.Constant;
mtrend = EstVAR.Trend;


Trend_starArray_{j}= mtrend;
Const_starArray{j} = mconst;
mP1_starArray{j}= mPi1c;
mP2_starArray{j} = mPi2c;
mP3_starArray{j} = mPi3c;
mP4_starArray{j} = mPi4c;

% Define covariance matrix in the BS sample

Sigma_u = EstVAR.Covariance;
Sigma_u_starArray{j} = Sigma_u;

options = optimset('MaxFunEvals', 200000,'TolFun',1e-500,'MaxIter', 200000,'TolX',1e-50);   
[StructuralParam_Estimation_MATRIX,Likelihood_SVAR,exitflag,output,grad,Hessian_MATRIX] = fminunc('Likelihood_SVAR2',InitialValue',options);
SE_Hessian_MATRIX = diag(inv(Hessian_MATRIX)).^0.5;  % computes Hessian-based standard errors


% Estimated parameters
A = eye(size(Sigma_u,1),size(Sigma_u,1));
A(1,3) = -2.08;
B = zeros(size(Sigma_u,1),size(Sigma_u,1));

SE_A = zeros(size(Sigma_u,1),size(Sigma_u,1));
SE_B = zeros(size(Sigma_u,1),size(Sigma_u,1));

HSelectionA = zeros(M*M,StructuralParamA);       % matrix S_A
HSelectionB = zeros(M*M,StructuralParamB);       % matrix S_B

for a_par = 1 : size(ParamNumberA,1)
    A(ParamNumberA(a_par,1)) = StructuralParam_Estimation_MATRIX(a_par); % puts the estimated elements of A in the right place       
    SE_A(ParamNumberA(a_par,1))= SE_Hessian_MATRIX(a_par);
    HSelectionA(ParamNumberA(a_par,1),a_par) = 1;                  % puts "1" in the correct place of selection matrix S_A 
end

l = StructuralParam_Estimation_MATRIX(size(ParamNumberA,1)+1:StructuralParam);
for b_par = 1:size(ParamNumberB,1)
    B(ParamNumberB(b_par,1)) = l(b_par);
    SE_B(ParamNumberB(b_par,1))= SE_Hessian_MATRIX(b_par);
    HSelectionB(ParamNumberB(b_par,1),b_par) = 1;
end        

% Sign normalization
if B(1,1)<0
   B(:,1)=-B(:,1);
end

if B(2,2)<0
   B(:,2)=-B(:,2); 
end

if B(3,3)<0
   B(:,3)=-B(:,3);
end    

HorizonIRF = 20;

C_IRF = A^(-1)*B;                                   % instantaneous impact at h=0
J=[eye(M) zeros(M,M*(p-1))];                        % selection matrix J used in IRF computation 
CompanionMatrix = [mPi1c mPi2c mPi3c mPi4c;                 % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
     

    for h = 0 : HorizonIRF
    TETA(:,:,h+1)=J*CompanionMatrix^h*J'*C_IRF;
    end


for k = 0 : HorizonIRF       
    tax_tax(1,k+1) = TETA(1,1,k+1);
    spend_tax(1,k+1) = TETA(1,2,k+1);
    gdp_tax(1,k+1) = TETA(1,3,k+1);  
    tax_spend(1,k+1) = TETA(2,1,k+1);
    spend_spend(1,k+1) = TETA(2,2,k+1);
    gdp_spend(1,k+1) = TETA(2,3,k+1);
    tax_gdp(1,k+1) = TETA(3,1,k+1);
    spend_gdp(1,k+1) = TETA(3,2,k+1);
    gdp_gdp(1,k+1) = TETA(3,3,k+1);
end

TETA_starArray{j} = TETA;
taxongdp_starArray{j} = tax_gdp;
spendingongdp_starArray{j} = spend_gdp;

MultTAX = (tax_gdp./tax_tax(1,1))*mratioGDPTAX;
MultG = (spend_gdp./spend_spend(1,1))*mratioGDPG;

DYNAMICMULTIPLIERTAX_starArray{j} = MultTAX;
DYNAMICMULTIPLIERSPENDING_starArray{j} = MultG;
end

%%
DYNAMICMULTIPLIERTAX_starMatrix = cell2mat(DYNAMICMULTIPLIERTAX_starArray);
DYNAMICMULTIPLIERSPENDING_starMatrix = cell2mat(DYNAMICMULTIPLIERSPENDING_starArray);
alpha = 0.1;
lower_bound__tax = zeros(HorizonIRF+1,1);
upper_bound_tax = zeros(HorizonIRF+1,1);
lower_bound__spending = zeros(HorizonIRF+1,1);
upper_bound_spending = zeros(HorizonIRF+1,1);

%%
for a = 1:HorizonIRF+1
lower_bound_tax(a,:) = quantile (DYNAMICMULTIPLIERTAX_starMatrix(:,a), alpha/2);
upper_bound_tax(a,:) = quantile (DYNAMICMULTIPLIERTAX_starMatrix(:,a), 1-(alpha/2));
lower_bound_spending(a,:) = quantile (DYNAMICMULTIPLIERSPENDING_starMatrix(:,a), alpha/2);
upper_bound_spending(a,:) = quantile (DYNAMICMULTIPLIERSPENDING_starMatrix(:,a), 1-(alpha/2));
end


%% Dynamic Multipliers - TAX
plot(0:HorizonIRF, multTAX_noCoint, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
multTAX_noCoint = plot(0:HorizonIRF, upper_bound_tax, '--b')
xlabel('Quarters')
title('Tax Dynamic Multiplier')
grid on
saveas(multTAX_noCoint, 'multTAX_noCoint', 'png');

%% Dynamic Multipliers - G
plot(0:HorizonIRF, multG_noCoint, 'LineWidth',2)
hold on
plot(0:HorizonIRF, lower_bound_spending, '--b')
hold on
multG_noCoint = plot(0:HorizonIRF, upper_bound_spending, '--b')
xlabel('Quarters')
title('Spending Dynamic Multiplier')
grid on
saveas(multG_noCoint, 'multG_noCoint', 'png');

%% Defining bands IRF
taxongdp_starMatrix = cell2mat (taxongdp_starArray);
spendingongdp_starMatrix = cell2mat (spendingongdp_starArray);
alpha = 0.1;

lower_bound_tax = zeros(HorizonIRF+1,1);
upper_bound_tax = zeros(HorizonIRF+1,1);

lower_bound_spending = zeros(HorizonIRF+1,1);
upper_bound_spending = zeros(HorizonIRF+1,1);

for a = 1:HorizonIRF+1
lower_bound_tax(a,:) = quantile (taxongdp_starMatrix(:,a), alpha/2);
upper_bound_tax(a,:) = quantile (taxongdp_starMatrix(:,a), 1-(alpha/2));
lower_bound_spending(a,:) = quantile (spendingongdp_starMatrix(:,a), alpha/2);
upper_bound_spending(a,:) = quantile (spendingongdp_starMatrix(:,a), 1-(alpha/2));
end

%% Figure IRF tax_dgp
plot(0:HorizonIRF, tax_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
tax_gdp_noCoint = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('IRF tax on GDP')
grid on
saveas(tax_gdp_noCoint, 'tax_gdp_noCoint', 'png')

%% Figure IRF spend_gdp
plot(0:HorizonIRF, spend_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
spend_gdp_noCoint = plot(0:HorizonIRF, upper_bound_spending, '--b');
xlabel('Quarters')
title('IRF spending on GDP')
grid on
saveas(spend_gdp_noCoint, 'spend_gdp_noCoint', 'png')

%% Rank conditions    
TopZero=zeros(M^2,StructuralParamB);
BotZero=zeros(M^2,StructuralParamA);

HSEL=[HSelectionA,TopZero;BotZero,HSelectionB];

Sigma_kron_A1=kron(Sigma_u,A^(-1));
AB_kron=kron(A^(-1)*B,A^(-1));
BigMatrix=[-kron(eye(1),mDD)*Sigma_kron_A1,kron(eye(1),mDD)*AB_kron];
Jacobian=2*BigMatrix*HSEL;
rank(Jacobian)
size(ParamNumber,1) %full rank ==> local identification

