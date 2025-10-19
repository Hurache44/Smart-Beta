


%%%%%%%%%%%%%%%%% 
% 
% Smart Beta Portfolio
%
%%%%%%%%%%%%%%%%%

% Error: beta market one not oatly

%% Data

% Check

% reading data
data = readtable('/Users/robertlange/Desktop/Investment/Matlab/Returns.xlsx');
dates = data{:,1};
assetNames = data.Properties.VariableNames(2:end);
prices = data{:,2:end};
nAssets = size(prices,2);
nObs = size(prices,1);
returns = diff(prices)./prices(1:end-1,:)*100;

C = cov(returns,'omitrows');
% disp('Covariance of returns (sample):');
% disp(C);

CorrMatrix = corrcov(C);
% disp('Correlation from Covariance:');
% disp(CorrMatrix);
%%
myWeights = [0.236 0.216 0.073 0.031 0.045 0.035 0.035 0.033 0.053 0.033 0.035 0.048 0.039 0.010 0.010 0.008 0.01 0.01 0.01 0.01 0.01 0.01];
myWeights = myWeights(:);
myWeights = myWeights/sum(myWeights)
covR = cov(returns)+1e-6*eye(nAssets);

disp('Initial portfolio weights: ');
disp(table(assetNames',myWeights,'VariableNames',{'Asset','Weight'}));
fprintf('Sum of weigths: %.4f\n\n', sum(myWeights));

portVar = myWeights'*covR*myWeights;
portStd = sqrt(portVar);
fprintf('Portfolio variance: %.6e\n ',portVar);
fprintf('Portfolio standard deviation: %.6e\n ', portStd);

%%

figure;
t1 = tiledlayout(4,4,'TileSpacing','Compact');
for i = 1:min(16,nAssets)
    nexttile
    plot(dates,prices(:,i))
    title(strrep(assetNames{i},'_',' '))
    datetick('x','keeplimits')
end

if nAssets > 16
    figure;
    t2 = tiledlayout(4,4,'TileSpacing','Compact');
    for i = 17:nAssets
        nexttile
        plot(dates,prices(:,i))
        title(strrep(assetNames{i},'_',' '))
        datetick('x','keeplimits')
    end
end

%%
marginalRiskContribution = @(wts,covMat) (wts(:).*(covMat*wts(:)))/sqrt(wts(:)'*covMat*wts(:));
riskContribution = marginalRiskContribution(myWeights,covR);
totalRisk = sqrt(myWeights'*covR*myWeights);
pctRiskContribution = 100*riskContribution/totalRisk;

figure;
bar(riskContribution);
title('Marginal Risk Contribution (MRC)');
set(gca,'XTick',1:nAssets,'XTickLabel',assetNames);
xtickangle(45);
grid on;
axis tight;

disp(table(assetNames',pctRiskContribution, ...
    'VariableNames',{'Asset','PctRiskContribution'}));


% Rolling Portfolio Risk
numInSample = 250;
disp('=== Optimization: Minimum Variance Portfolio ===');
nAsset = size(covR,1);
% meanR = mean(returns);
    % covR = cov(returns)+1e-6*eye(nAsset);
% corR = covR+1e6*eye(size(covR));

% fmincom problem
% o = optimoptions('fmincon','Algorithm','interior-point','TolFun',1e-6,'Display','off');

% function handle
% fH = @(x) riskCostFunction(x,covR);

% Initial Condition
x0 = ones(nAsset,1) / nAsset;

% Constraints for sum, weights == 1
Aeq = ones(1,nAsset);
Beq = 1;
LB = zeros(1,nAsset);
UB = ones(1,nAsset);

opts = optimoptions('fmincon','Algorithm','interior-point','Display','off','TolFun',1e-8);

[optWeights,fVal] = ...
    fmincon(@(x) riskCostFunction(x,covR),x0,[],[],Aeq,Beq,LB,UB,[],opts);

fprintf('Min variance objective: %.6e\n',fVal);
disp(table(assetNames',optWeights,'VariableNames',{'Asset','OptimizedWeight'}));
fprintf('Sum of optimized weights: %.4f\n',sum(optWeights));
%%
riskContribution = marginalRiskContribution(myWeights,covR)

% marching this portfolio forwards
window = numInSample;
nSteps = size(returns,1)-window;
riskRolling = zeros(nSteps,1);
mrcOvertTime = zeros(nSteps,nAssets);

for t = 1:nSteps
    retWindow = returns(t:t+window-1,:);
    covWin = cov(retWindow)+1e-6*eye(nAsset);
    riskRolling(t) = sqrt(myWeights'*covWin*myWeights);
    mrcOverTime(t,:) = marginalRiskContribution(myWeights,covWin);
end

% Plot the evolution of the portfolio risk
figure;
plot(dates(window+1:window+nSteps),riskRolling,'LineWidth',1.5);
title('Rolling Portfolio Risk (Std Dev)');
xlabel('Date');
ylabel('Std dev (Port vol)');
grid on;
%%
% Plot the individual risk
figure
surf(1:nAsset,1:nSteps,mrcOverTime);
shading interp;
axis tight;
view([115 60]);
ylabel('Time step');
xlabel('Asset index');
zlabel('MRC');
set(gca,'XTick',1:nAsset,'XTickLabel',assetNames);
xtickangle(4);
%%
% Beta Hedge

marketIdx = 1;
% asset betas relative to market (full sample returns)
 mr = returns(:,marketIdx);
demeanedReturns = returns - mean(returns);
demeanedMarket = mr-mean(mr);
assetCovWithMarket = (demeanedReturns'*demeanedMarket)/(length(mr)-1);
varMarket = var(mr);
beta = assetCovWithMarket / varMarket;

% for i = 1:nAsset
%     betas(i) = cov(returns(:,i), mr) /  var_m;
% end
% Portfolio Beta with current weights
portfolioBeta = myWeights'*beta;

% Position in the market
% total combined beta
hedgeWeight = -portfolioBeta;

% Beta hedge additional to existing portfolio value
fprintf('Portfolio beta (pre-hedge): %.4f\n',portfolioBeta);
fprintf('Hedge weight in market (seperate position): %.4f\n',hedgeWeight);

% Risk pre and post hedge (with extenden covariance)
% hedged variance

    cov_assets = cov(returns);
    covP_M = myWeights'*assetCovWithMarket;
    varP = myWeights'*cov_assets*myWeights;
    varHedged = varP + (hedgeWeight^2)*varMarket+2*hedgeWeight*covP_M;
    fprintf('Portfolio std dev (pre-hedge): %.6f\n',sqrt(varP));
    fprintf('Portfolio std dev (post-hedge): %.6f\n',sqrt(varHedged));

    result = table(assetNames',myWeights,riskContribution,beta(:), ...
        'VariableNames',{'Asset','Weight','MRC','Beta'});
    disp(result);

    % == post weights result == %
    result_post = result;
    resul_post.MRC = result_post.MRC*(1+hedgeWeight);

    f1 = figure('Name','Pre Hedge Portfolio','NumberTitle','off');
    displayDataPre = [result.Asset,num2cell([result.Weight,result.MRC,result.Beta])];
    t1 = uitable('Parent',f1,'Data',displayDataPre, ...
        'ColumnName',result.Properties.VariableNames, ...
        'Units','normalized','Position',[0 0 1 1]);

    % Post Hedge
    % result_post = result;
    % result_post.MRC = result_post.MRC*(1+hedgeWeight);
    % displayDataPost1 = [result_post.Asset,num2cell([result_post.Weight,result_post.MRC,result_post.Beta])];

 f2 = figure('Name','Post-Hedge Portfolio','NumberTitle','off');
    displayDataPost = [result.Asset,num2cell([result_post.Weight,result_post.MRC,result_post.Beta])];
    t2 = uitable('Parent',f2,'Data',displayDataPost, ...
        'ColumnName',result_post.Properties.VariableNames, ...
        'Units','normalized','Position',[0 0 1 1]);
%%
    %subplot(2,1,2);
    %uit2 = uitable('Data',displayDataPost1, ...
    %    'ColumnName',result_post.Properties.VariableNames, ...
    %    'Units','normalized', ...
    %    'Position',[0 0.1 1 0.8]);
    %title('Post-Hedge Portfolio');

%%
% disp(size(assetNames'));
% disp(size(myWeights));
% disp(size(riskContribution));
% disp(size(beta));

%%
%%%
% Stock Selection Models
%%%

% CalculateDividendYield % function

% CalculatePriceMomentum % function

%%
%%%%%%%%%%%%%%%%% 
% Portfolio Analysis with Hedge Effects
%%%%%%%%%%%%%%%%%

%% Load Data
data = readtable('/Users/robertlange/Desktop/Investment/Matlab/Returns.xlsx');
dates = data{:,1};
assetNames = data.Properties.VariableNames(2:end);
prices = data{:,2:end};
nAssets = size(prices,2);

returns = diff(prices)./prices(1:end-1,:)*100;  % daily returns %

covR = cov(returns,'omitrows');

%% Portfolio Weights
myWeights = [0.236 0.216 0.073 0.031 0.045 0.035 0.035 0.033 0.053 0.033 ...
             0.035 0.048 0.039 0.010 0.010 0.008 0.01 0.01 0.01 0.01 0.01 0.01]';
myWeights = myWeights / sum(myWeights);

%% Portfolio Metrics Pre-Hedge
riskContribution = (myWeights .* (covR*myWeights)) / sqrt(myWeights'*covR*myWeights);
totalRisk = sqrt(myWeights'*covR*myWeights);
pctRiskContribution = 100*riskContribution/totalRisk;

%% Beta Hedge
marketIdx = 1; % assume market is first asset
mr = returns(:,marketIdx);
demeanedReturns = returns - mean(returns);
demeanedMarket = mr - mean(mr);
assetCovWithMarket = (demeanedReturns'*demeanedMarket)/(length(mr)-1);
varMarket = var(mr);
beta = assetCovWithMarket / varMarket;

portfolioBeta = myWeights'*beta;
hedgeWeight = -portfolioBeta; % simple market hedge

% Post-Hedge Portfolio
covP_M = myWeights'*assetCovWithMarket;
varP = myWeights'*covR*myWeights;
varHedged = varP + (hedgeWeight^2)*varMarket + 2*hedgeWeight*covP_M;
totalRiskPost = sqrt(varHedged);

% Post-Hedge MRC scaling (simplistic, linear approximation)
mrcPost = riskContribution * (1 + hedgeWeight);

%% Table of Pre/Post Hedge
result = table(assetNames', myWeights, riskContribution, beta(:), ...
    'VariableNames',{'Asset','Weight','MRC','Beta'});
result_post = result;
result_post.MRC = mrcPost;

%% Figure: Bar charts Pre vs Post Hedge
f = figure('Name','MRC Pre vs Post Hedge','NumberTitle','off');
subplot(2,1,1);
bar([result.MRC, result_post.MRC]);
title('Marginal Risk Contribution: Pre vs Post Hedge');
set(gca,'XTick',1:nAssets,'XTickLabel',assetNames);
xtickangle(90);
ylabel('MRC');
legend('Pre-Hedge','Post-Hedge');

subplot(2,1,2);
bar([result.Beta, result_post.Beta]);
title('Beta per Asset: Pre vs Post Hedge');
set(gca,'XTick',1:nAssets,'XTickLabel',assetNames);
xtickangle(90);
ylabel('Beta');
legend('Pre-Hedge','Post-Hedge');

%% Table of Portfolio Hedge Recommendation
hedgeTable = table;
hedgeTable.Market = assetNames(marketIdx);
hedgeTable.PortfolioBeta = portfolioBeta;
hedgeTable.HedgeWeight = hedgeWeight;
hedgeTable.RiskReduction = totalRisk - totalRiskPost;
hedgeTable.PostHedgeStd = totalRiskPost;

disp('===== Suggested Hedge for Portfolio =====');
disp(hedgeTable);

%% Optional: Dynamic Rolling Risk (Pre/Post Hedge)
numInSample = 250;
nSteps = size(returns,1) - numInSample;
riskRollingPre = zeros(nSteps,1);
riskRollingPost = zeros(nSteps,1);

for t = 1:nSteps
    retWindow = returns(t:t+numInSample-1,:);
    covWin = cov(retWindow) + 1e-6*eye(nAssets);
    riskRollingPre(t) = sqrt(myWeights'*covWin*myWeights);
    
    % Post-hedge variance
    covP_M = myWeights'*((retWindow - mean(retWindow))'*(retWindow(:,marketIdx)-mean(retWindow(:,marketIdx)))/(numInSample-1));
    varWin = myWeights'*covWin*myWeights;
    varHedgedWin = varWin + (hedgeWeight^2)*var(covWin(:,marketIdx)) + 2*hedgeWeight*covP_M;
    riskRollingPost(t) = sqrt(varHedgedWin);
end

% Plot rolling risk
% Ensure equal length of time and risk series
plotDates = dates(numInSample+1:numInSample+nSteps);

figure;
plot(plotDates, riskRollingPre, 'LineWidth', 1.5); hold on;
plot(plotDates, riskRollingPost, 'LineWidth', 1.5);
title('Rolling Portfolio Risk: Pre vs Post Hedge');
xlabel('Date'); 
ylabel('Standard Deviation');
legend('Pre-Hedge','Post-Hedge','Location','best');
grid on;

%%
% === Additional Metrics ===

disp('=== Additional Metrics  ===');
nAssets = size(returns,2);

% Risk Contribution 
w = myWeights(:);
Sigma = covR;
Sigma_w = Sigma*w;
compVar = w.*Sigma_w;
portVar = w'*Sigma*w;
portStd = sqrt(portVar);

pctVar = 100*compVar/portVar;
mrc = (w.*Sigma_w)/portStd;
pctRisk = 100*mrc/portStd;
%%
riskTable = table(assetNames',w,compVar,pctVar,mrc,pctRisk, ...
    'VariableNames',{'Asset','Weight','VarContribution','PctOfVar','MRC','PctOfRisk'});
disp('--- Risk contribution table ---');
disp(riskTable);

figure('Name','Risk contribution per asset (variance & percent)');
subplot(2,1,1);
bar(compVar);
title('Absolute contribution to portfolio variance (units: %^2)');
set(gca,'XTick',1:nAssets,'XTickLabel',assetNames);
xtickangle(44);
grid on;

subplot(2,1,2);
bar(pctVar);
title('% Contribution to total portfolio variance');
set(gca,'XTick',1:nAssets,'XTickLabel',assetNames);
xtickangle(44);
grid on;

%%
% Automatic hedge-scan and analytical optimal hedge
cov_vec = assetCovWithMarket(:);
covP_M = w'*cov_vec;
varP = w'*Sigma*w;

h_star = -covP_M/varMarket;
% Portfolio variance as function of hedge weight
var_at_hstar = varP+h_star^2*varMarket+2*h_star*covP_M;
std_at_hstar = sqrt(var_at_hstar);

h_grid = linspace(h_star*2,h_star*(-2),501);
var_grid = varP+(h_grid.^2)*varMarket+2.*h_grid*covP_M;
std_grid = sqrt(var_grid);

% plot the result
figure('Name','Hedge scan: portfolio std contrary hedge weight');
plot(h_grid,std_grid,'LineWidth',1.5);
hold on;
plot(h_star,sqrt(var_at_hstar),'ro','MarkerFaceColor','r');
xlabel('Hedge weight (position in market; negative = hedge/short)');
ylabel('Portfolio Std Dev (%)');
title('Portfolio Std Dev as function of height weight');
grid on;
legend('std(h)','analytic optimal');

fprintf('\nAnalytic optimal hedge h*= %.6f (reduces std to %.6f',h_star,sqrt(var_at_hstar));

% recommended hedge
hedgeSuggestion = table(assetNames(marketIdx),portfolioBeta,h_star,sqrt(varP),sqrt(var_at_hstar), ...
    'VariableNames',{'MarketAsset','PortfolioBeta','HedgeWeight','StdPre','StdAtOptHedge'});
disp('--- Hedge suggestion ---');
disp(hedgeSuggestion);

% Tracking Error contr. benchmark (static and rolling)
R_p = returns*w;
R_p_hedged = R_p + hedgeWeight * mr;

activePre = R_p-mr;
activePost = R_p_hedged - mr;
TE_pre = std(activePre);
TE_post = std(activePost);

fprintf('\nStatic Tracking Error (pre) = %.6f; (post) = %.6f\n', TE_pre, TE_post);

% Rolling window tracking error
rollWindow = 250;
nT  = size(returns,1);
nStepsTE = nT - rollWindow;
TE_roll_pre = nan(nStepsTE,1);
TE_roll_post = nan(nStepsTE,1);

for t = 1:nStepsTE
    idx = t:(t+rollWindow-1);
    TE_roll_pre(t) = std(R_p(idx)-mr(idx));
    TE_roll_post(t) = std(R_p_hedged(idx)-mr(idx));
end

% Plot tracking error
plotDates_TE = dates(rollWindow+1:rollWindow+nStepsTE);
figure('Name','Tracking Error contr. Market');
plot(plotDates_TE,TE_roll_pre,'LineWidth',1.2);
hold on;
plot(plotDates_TE,TE_roll_post,'LineWidth',1.2);
legend('TE Pre','TE Post','Location','best');
grid on;
title('Rolling Tracking Error (window = 250)');
xlabel('Date');
ylabel('Tracking Error (%)');

% Visualization riskRollingPre contr. riskRollingPost and shaded Difference
if ~exist('riskRollingPre','var') || numel(riskRollingPre) ~= nSteps
    nSteps = size(returns,1) - numInSample;
    riskRollingPre = zeros(nSteps,1);
    riskRollingPost = zeros(nSteps,1);
    for t = 1:nSteps
        idx = t:(t+numInSample-1);
        covWin = cov(returns(idx,:))+1e-6*eye(nAssets);
        riskRollingPre(t) = sqrt(w'*covWin*w);
        cov_vec_win = (returns(idx,:)-mean(returns(idx,:)))'*(returns(idx,marketIdx)-mean(returns(idx,marketIdx)))/(numInSample-1);
        covP_M_win = w'*cov_vec_win;
        varWin = w'*covWin*w;
        varHedgedWin = varWin+(h_star^2)*var(returns(idx,marketIdx))+2*h_star*covP_M_win;
        riskRollingPost(t) = sqrt(varHedgedWin);
    end
end

plotDate = dates(numInSample+1:numInSample+nSteps);
figure('Name','Rolling Risk Pre contr. Post and Reduction');
plot(plotDates,riskRollingPre,'-','LineWidth',1.5);
hold on;
plot(plotDates,riskRollingPost,'-','LineWidth',1.5);
% shaded area of reduction
x = [plotDates; flipud(plotDates)];
y = [riskRollingPre - riskRollingPost; flipud(zeros(size(riskRollingPost)))];
% y is the reduction (pre-post), ensure non-negativity
y_area = max(riskRollingPre - riskRollingPost,0);
area(plotDates,y_area,'FaceAlpha',0.2,'EdgeColor','none');
legend('Risk Pre','Risk Post','Risk Reduction (area)','Location','best');
title('Rolling Risk Pre contr. Post Hedge and Reduction');
xlabel('Date');
ylabel('Std Dev (%)');
grid on;

%%
%% KO PUT Hedge Module (Strike above for protection)
nAssets = size(prices,2);

% Parameters
leverage = 10;
portfolioValue = 1;
deltaATM_Put = -0.5;  % delta for put

% Daily volatility in price units
sigmaDaily = std(returns);                 % daily % returns
sigmaPrice = prices(end,:) .* sigmaDaily / 100;  % convert % to price units

% ----------------------
% KO PUT STRIKES (above price, for hedging downside)
% ----------------------
strikePut_1sigma = prices(end,:)' + sigmaPrice';      % column vector
strikePut_2sigma = prices(end,:)' + 2*sigmaPrice';    % column vector

% Hedge weights & notionals
hedgeWeights = w * h_star;          % nAssets x 1
notionals = abs(hedgeWeights * portfolioValue) ./ (abs(deltaATM_Put) * leverage);

KO_Hedges_Put = table(assetNames(:), ...
                  repmat({'KO Put'}, nAssets,1), ...
                  strikePut_1sigma, ...
                  strikePut_2sigma, ...
                  notionals, ...
                  repmat(leverage, nAssets,1), ...
                  hedgeWeights, ...
                  'VariableNames',{'Assets','HedgeType','Strike_1sigma','Strike_2sigma','Notional','Leverage','OptimalWeight'});

% Display table as figure
dataForTable_P = [KO_Hedges_Put.Assets, KO_Hedges_Put.HedgeType, ...
                num2cell(KO_Hedges_Put.Strike_1sigma), num2cell(KO_Hedges_Put.Strike_2sigma), ...
                num2cell(KO_Hedges_Put.Notional), num2cell(KO_Hedges_Put.Leverage), ...
                num2cell(KO_Hedges_Put.OptimalWeight)];

fHedge_Put = figure('Name','KO Put Hedge Suggestions','NumberTitle','off','Color','w');
tHedge_Put = uitable('Parent',fHedge_Put, ...
                 'Data', dataForTable_P, ...
                 'ColumnName', KO_Hedges_Put.Properties.VariableNames, ...
                 'Units','normalized', ...
                 'Position',[0 0 1 1]);

disp('--- Knock-Out Put Hedge Suggestion ---');
disp(KO_Hedges_Put);

% ----------------------
% KO CALL STRIKES (optional: below price)
% ----------------------
deltaATM_Call = 0.5;
strikeCall_1sigma = max(prices(end,:)' - sigmaPrice', 0);      % column vector
strikeCall_2sigma = max(prices(end,:)' - 2*sigmaPrice', 0);    % column vector

% Hedge weights & notionals
notionals_call = abs(hedgeWeights * portfolioValue) ./ (abs(deltaATM_Call) * leverage);

KO_Hedges_Call = table(assetNames(:), ...
                  repmat({'KO Call'}, nAssets,1), ...
                  strikeCall_1sigma, ...
                  strikeCall_2sigma, ...
                  notionals_call, ...
                  repmat(leverage, nAssets,1), ...
                  hedgeWeights, ...
                  'VariableNames',{'Assets','HedgeType','Strike_1sigma','Strike_2sigma','Notional','Leverage','OptimalWeight'});

% Display table as figure
dataForTable_C = [KO_Hedges_Call.Assets, KO_Hedges_Call.HedgeType, ...
                num2cell(KO_Hedges_Call.Strike_1sigma), num2cell(KO_Hedges_Call.Strike_2sigma), ...
                num2cell(KO_Hedges_Call.Notional), num2cell(KO_Hedges_Call.Leverage), ...
                num2cell(KO_Hedges_Call.OptimalWeight)];

fHedge_Call = figure('Name','KO Call Hedge Suggestions','NumberTitle','off','Color','w');
tHedge_Call = uitable('Parent',fHedge_Call, ...
                 'Data', dataForTable_C, ...
                 'ColumnName', KO_Hedges_Call.Properties.VariableNames, ...
                 'Units','normalized', ...
                 'Position',[0 0 1 1]);

disp('--- Knock-Out Call Hedge Suggestion ---');
disp(KO_Hedges_Call);
