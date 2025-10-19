function r = riskCostFunction(wts,covMat)
% function out = riskCostFunction(wts,covMat)
%
% cost function for the attempt to get equal risk-contribution portfolios

% Build up all the marginal risk contributions
wts = wts(:)/sum(wts);
r = wts'*covMat*wts;
end
