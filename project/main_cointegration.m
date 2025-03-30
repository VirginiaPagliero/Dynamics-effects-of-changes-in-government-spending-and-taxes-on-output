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

%% VARIABLES IN LEVELS 
exp_tax = exp(TAX);
exp_gdp = exp(GDP);
exp_g = exp(G);

ratio_tax = exp_tax./exp_gdp;
ratio_g = exp_g./exp_gdp;

%% ADF for tax 
[h,pValue,~] = adftest(TAX, model="TS", lag=3)

%% ADF for G
[h,pValue,~] = adftest(G, model="TS", lag=3)

%% ADF for GDP
[h,pValue] = adftest(GDP, model="TS", lag=3)

%% Johansen test
[h,pValue,stat] = jcitest(data, lags=3, Test=["trace" "maxeig"], model='H*');

%% COINTEGRATION - NO RESTRICTIONS
Mdl = vecm(3,2,3);
[EstMdl,EstSEMdl,logL,E] = estimate(Mdl,data,model="H*");
Sigma_u = EstMdl.Covariance;

%% VAR coefficients 
Pi1c = (EstMdl.Adjustment*EstMdl.Cointegration'+ eye(3)+ EstMdl.ShortRun{1,1});
Pi2c = -(EstMdl.ShortRun{1,1} - EstMdl.ShortRun{1,2});
Pi3c = -(EstMdl.ShortRun{1,2} - EstMdl.ShortRun{1,3});
Pi4c = -EstMdl.ShortRun{1,3};
const = (EstMdl.Adjustment*EstMdl.CointegrationConstant) + EstMdl.Constant;
trend = (EstMdl.Adjustment*EstMdl.CointegrationTrend) + EstMdl.Trend;

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

for a_par = 1 : size(ParamNumberA,1)
    A(ParamNumberA(a_par,1)) = StructuralParam_Estimation_MATRIX(a_par); 
    SE_A(ParamNumberA(a_par,1))= SE_Hessian_MATRIX(a_par);
    HSelectionA(ParamNumberA(a_par,1),a_par) = 1;                  
end

i = StructuralParam_Estimation_MATRIX(size(ParamNumberA,1)+1:StructuralParam);
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

%% IRF 
HorizonIRF = 20;

C_IRF = A^(-1)*B;                                   % instantaneous impact at h=0
J=[eye(M) zeros(M,M*(p-1))];                        % selection matrix J used in IRF computation 
CompanionMatrix = [Pi1c Pi2c Pi3c Pi4c;                 % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA(:,:,h+1)=J*(CompanionMatrix^h)*J'*C_IRF;
    end
    
TETA 

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

MultTAX_coint = (tax_gdp./tax_tax(1,1))*mratioGDPTAX; 
MultG_coint = (spend_gdp./spend_spend(1,1))*mratioGDPG; 

%% Bootstrap
mean_residuals = mean(E);

Recentered_residuals = E - mean_residuals;

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
   bootstrappedData(t, :) = const + trend*t + Pi1c*(bootstrappedData(t-1, :)') + Pi2c*(bootstrappedData(t-2, :)') + Pi3c*(bootstrappedData(t-3, :)') + Pi4c*(bootstrappedData(t-4, :)') + resampledResiduals(t, :)';
end
bootstrappedDataArray{j} = bootstrappedData;

Mdl = vecm(3,2,3);
[EstMdl,EstSEMdl,logL,E] = estimate(Mdl,bootstrappedData,model="H*");

mPi1c = (EstMdl.Adjustment*EstMdl.Cointegration'+ eye(3)+ EstMdl.ShortRun{1,1});
mPi2c = -(EstMdl.ShortRun{1,1}-EstMdl.ShortRun{1,2});
mPi3c = -(EstMdl.ShortRun{1,2}-EstMdl.ShortRun{1,3});
mPi4c = -EstMdl.ShortRun{1,3};
mconst = (EstMdl.Adjustment*EstMdl.CointegrationConstant) + EstMdl.Constant;
mtrend = (EstMdl.Adjustment*EstMdl.CointegrationTrend);


Trend_starArray_{j}= mtrend;
Const_starArray{j} = mconst;
mP1_starArray{j}= mPi1c;
mP2_starArray{j} = mPi2c;
mP3_starArray{j} = mPi3c;
mP4_starArray{j} = mPi4c;

% Define covariance matrix in the BS sample

Sigma_u = EstMdl.Covariance;
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
figure
plot(0:HorizonIRF, MultTAX_coint, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
multTAX_coint = plot(0:HorizonIRF, upper_bound_tax, '--b')
xlabel('Quarters')
title('Tax Dynamic Multiplier')
grid on
saveas(multTAX_coint, 'multTAX_coint', 'png');

%% Dynamic Multipliers - G
figure
plot(0:HorizonIRF, MultG_coint, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
multG_coint = plot(0:HorizonIRF, upper_bound_spending, '--b');
xlabel('Quarters')
title('Spending Dynamic Multiplier')
grid on
saveas(multG_coint, 'multG_coint', 'png')


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

%% Figure IRF tax_dp
plot(0:HorizonIRF, tax_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
tax_gdp_coint = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('IRF tax on GDP')
grid on
saveas(tax_gdp_coint, 'tax_gdp_coint', 'png');

%% Figure IRF spend_gdp
plot(0:HorizonIRF, spend_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
spend_gdp_coint = plot(0:HorizonIRF, upper_bound_spending, '--b');
xlabel('Quarters')
title('IRF spending on GDP')
grid on
saveas(spend_gdp_coint, 'spend_gdp_coint', 'png')

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



%% COINTEGRATION - ANGELINI ET AL., 2022
%% Estimate VECM 
Cointegration = [1 0; 0 1; -1 -1.1808];
CointTrend = [0; 0.0020826];
Mdl2 = vecm(3,2,3);
Mdl2.Cointegration = Cointegration;
Mdl2.CointegrationTrend = CointTrend;
[EstMdl2,EstSEMdl2,logL2,E2] = estimate(Mdl2,data,model="H*");
Sigma_u = EstMdl2.Covariance;

%% Plot
beta = EstMdl2.Cointegration;
for i = 1: 119
    y_t(i,:) = data(i,:);
    y(i,:) = beta'*y_t(i,:)';
end

figure
plot(y(:,1),'LineWidth',2);
hold on
coint2_series = plot(y(:,2),'LineWidth',2);
hold on
legend('TAX - GDP','G - (1.18)GDP + (0.002)t')
xlabel('Quarters')
title '{\bf Cointegrated series}';
axis tight
grid on
saveas(coint2_series,'coint2_series','png'); 

%% VAR coefficients 
Pi1c = (EstMdl2.Adjustment*EstMdl2.Cointegration'+ eye(3)+ EstMdl2.ShortRun{1,1});
Pi2c = -(EstMdl2.ShortRun{1,1} - EstMdl2.ShortRun{1,2});
Pi3c = -(EstMdl2.ShortRun{1,2} - EstMdl2.ShortRun{1,3});
Pi4c = -EstMdl2.ShortRun{1,3};
const = (EstMdl2.Adjustment*EstMdl2.CointegrationConstant) + EstMdl2.Constant;
trend = (EstMdl2.Adjustment*EstMdl2.CointegrationTrend) + EstMdl2.Trend;

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

for a_par = 1 : size(ParamNumberA,1)
    A(ParamNumberA(a_par,1)) = StructuralParam_Estimation_MATRIX(a_par); % puts the estimated elements of A in the right place       
    SE_A(ParamNumberA(a_par,1))= SE_Hessian_MATRIX(a_par);
    HSelectionA(ParamNumberA(a_par,1),a_par) = 1;                  % puts "1" in the correct place of selection matrix S_A 
end

i = StructuralParam_Estimation_MATRIX(size(ParamNumberA,1)+1:StructuralParam);
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

%% IRF 
HorizonIRF = 20;

C_IRF = A^(-1)*B;                                   % instantaneous impact at h=0
J=[eye(M) zeros(M,M*(p-1))];                        % selection matrix J used in IRF computation 
CompanionMatrix = [Pi1c Pi2c Pi3c Pi4c;                 % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA(:,:,h+1)=J*(CompanionMatrix^h)*J'*C_IRF;
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

MultTAX_coint2 = (tax_gdp/tax_tax(1,1))*mratioGDPTAX; 
MultG_coint2 = (spend_gdp/spend_spend(1,1))*mratioGDPG;

%% Bootstrap
mean_residuals = mean(E2);

Recentered_residuals = E2 - mean_residuals;

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
   bootstrappedData(t, :) = const + trend*t + Pi1c*(bootstrappedData(t-1, :)') + Pi2c*(bootstrappedData(t-2, :)') + Pi3c*(bootstrappedData(t-3, :)') + Pi4c*(bootstrappedData(t-4, :)') + resampledResiduals(t, :)';
end
bootstrappedDataArray{j} = bootstrappedData;

Mdl = vecm(3,2,3);
Mdl.Cointegration = Cointegration;
Mdl.CointegrationTrend = CointTrend;
[EstMdl,EstSEMdl,logL,E] = estimate(Mdl,bootstrappedData,model="H*");

mPi1c = (EstMdl.Adjustment*EstMdl.Cointegration'+ eye(3)+ EstMdl.ShortRun{1,1});
mPi2c = -(EstMdl.ShortRun{1,1}-EstMdl.ShortRun{1,2});
mPi3c = -(EstMdl.ShortRun{1,2}-EstMdl.ShortRun{1,3});
mPi4c = -EstMdl.ShortRun{1,3};
mconst = (EstMdl.Adjustment*EstMdl.CointegrationConstant) + EstMdl.Constant;
mtrend = (EstMdl.Adjustment*EstMdl.CointegrationTrend);


Trend_starArray_{j}= mtrend;
Const_starArray{j} = mconst;
mP1_starArray{j}= mPi1c;
mP2_starArray{j} = mPi2c;
mP3_starArray{j} = mPi3c;
mP4_starArray{j} = mPi4c;

% Define covariance matrix in the BS sample

Sigma_u = EstMdl.Covariance;
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
plot(0:HorizonIRF, MultTAX_coint2, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
multTAX_coint2 = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('Tax Dynamic Multiplier')
grid on
saveas(multTAX_coint2, 'multTAX_coint2', 'png')

%% Dynamic Multipliers - G
plot(0:HorizonIRF, MultG_coint2, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
multG_coint2 = plot(0:HorizonIRF, upper_bound_spending, '--b');
title('Spending Dynamic Multiplier')
grid on
saveas(multG_coint2, 'multG_coint2', 'png')

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

%% Figure IRF tax_dp
plot(0:HorizonIRF, tax_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
tax_gdp_coint2 = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('IRF tax on GDP')
grid on
saveas(tax_gdp_coint2, 'tax_gdp_coint2', 'png');

%% Figure IRF spend_gdp
plot(0:HorizonIRF, spend_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
spend_gdp_coint2 = plot(0:HorizonIRF, upper_bound_spending, '--b');
xlabel('Quarters')
title('IRF spending on GDP')
grid on
saveas(spend_gdp_coint2, 'spend_gdp_coint2', 'png');

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



%% COINTEGRATION - BLANCHARD E PEROTTI, 2002
%% Estimate VECM 
Cointegration = [1; -1; 0];
CointTrend = [-0.0014365];
Mdl2 = vecm(3,1,3);
Mdl2.Cointegration = Cointegration;
Mdl2.CointegrationTrend = CointTrend;
[EstMdl2,EstSEMdl2,logL2,E2] = estimate(Mdl2,data,model="H*");
Sigma_u = EstMdl2.Covariance;

%% Plot
beta = EstMdl2.Cointegration;
for i = 1: 119
    y_t(i,:) = data(i,:);
    y(i,:) = beta'*y_t(i,:)';
end

figure
coint3_series = plot(y(:,1),'LineWidth',2);
legend('TAX - G - (0.0014) t')
xlabel('Quarters')
title '{\bf Cointegrated series}';
axis tight
grid on
saveas(coint3_series,'coint3_series','png'); 

%% VAR coefficients 
Pi1c = (EstMdl2.Adjustment*EstMdl2.Cointegration'+ eye(3)+ EstMdl2.ShortRun{1,1});
Pi2c = -(EstMdl2.ShortRun{1,1} - EstMdl2.ShortRun{1,2});
Pi3c = -(EstMdl2.ShortRun{1,2} - EstMdl2.ShortRun{1,3});
Pi4c = -EstMdl2.ShortRun{1,3};
const = (EstMdl2.Adjustment*EstMdl2.CointegrationConstant) + EstMdl2.Constant;
trend = (EstMdl2.Adjustment*EstMdl2.CointegrationTrend) + EstMdl2.Trend;

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

for a_par = 1 : size(ParamNumberA,1)
    A(ParamNumberA(a_par,1)) = StructuralParam_Estimation_MATRIX(a_par); % puts the estimated elements of A in the right place       
    SE_A(ParamNumberA(a_par,1))= SE_Hessian_MATRIX(a_par);
    HSelectionA(ParamNumberA(a_par,1),a_par) = 1;                  % puts "1" in the correct place of selection matrix S_A 
end

i = StructuralParam_Estimation_MATRIX(size(ParamNumberA,1)+1:StructuralParam);
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

%% IRF 
HorizonIRF = 20;

C_IRF = A^(-1)*B;                                   % instantaneous impact at h=0
J=[eye(M) zeros(M,M*(p-1))];                        % selection matrix J used in IRF computation 
CompanionMatrix = [Pi1c Pi2c Pi3c Pi4c;                 % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA(:,:,h+1)=J*(CompanionMatrix^h)*J'*C_IRF;
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

MultTAX_coint3 = (tax_gdp./tax_tax(1,1))*mratioGDPTAX;  
MultG_coint3 = (spend_gdp./spend_spend(1,1))*mratioGDPG; 

%% Bootstrap
mean_residuals = mean(E2);

Recentered_residuals = E2 - mean_residuals;

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
   bootstrappedData(t, :) = const + trend*t + Pi1c*(bootstrappedData(t-1, :)') + Pi2c*(bootstrappedData(t-2, :)') + Pi3c*(bootstrappedData(t-3, :)') + Pi4c*(bootstrappedData(t-4, :)') + resampledResiduals(t, :)';
end
bootstrappedDataArray{j} = bootstrappedData;

Mdl = vecm(3,1,3);
Mdl.Cointegration = Cointegration;
Mdl.CointegrationTrend = CointTrend;
[EstMdl,EstSEMdl,logL,E] = estimate(Mdl,bootstrappedData,model="H*");

mPi1c = (EstMdl.Adjustment*EstMdl.Cointegration'+ eye(3)+ EstMdl.ShortRun{1,1});
mPi2c = -(EstMdl.ShortRun{1,1}-EstMdl.ShortRun{1,2});
mPi3c = -(EstMdl.ShortRun{1,2}-EstMdl.ShortRun{1,3});
mPi4c = -EstMdl.ShortRun{1,3};
mconst = (EstMdl.Adjustment*EstMdl.CointegrationConstant) + EstMdl.Constant;
mtrend = (EstMdl.Adjustment*EstMdl.CointegrationTrend);

Trend_starArray_{j}= mtrend;
Const_starArray{j} = mconst;
mP1_starArray{j}= mPi1c;
mP2_starArray{j} = mPi2c;
mP3_starArray{j} = mPi3c;
mP4_starArray{j} = mPi4c;

% Define covariance matrix in the BS sample

Sigma_u = EstMdl.Covariance;
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
plot(0:HorizonIRF, MultTAX_coint3, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
multTAX_coint3 = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('Tax Dynamic Multiplier')
grid on
saveas(multTAX_coint3, 'multTAX_coint3', 'png');

%% Dynamic Multipliers - G
plot(0:HorizonIRF, MultG_coint3, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
multG_coint3 = plot(0:HorizonIRF, upper_bound_spending, '--b')
xlabel('Quarters')
title('Spending Dynamic Multiplier')
grid on
saveas(multG_coint3, 'multG_coint3', 'png');

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

%% Figure IRF tax_dp
plot(0:HorizonIRF, tax_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_tax, '--b')
hold on
tax_gdp_coint3 = plot(0:HorizonIRF, upper_bound_tax, '--b');
xlabel('Quarters')
title('IRF tax on GDP')
grid on
saveas(tax_gdp_coint3, 'tax_gdp_coint3', 'png')

%% Figure IRF spend_gdp
plot(0:HorizonIRF, spend_gdp, 'LineWidth',2)
hold on
plot(0:HorizonIRF,lower_bound_spending, '--b')
hold on
spend_gdp_coint3 = plot(0:HorizonIRF, upper_bound_spending, '--b');
xlabel('Quarters')
title('IRF spending on GDP')
grid on
saveas(spend_gdp_coint3, 'spend_gdp_coint3', 'png')

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


