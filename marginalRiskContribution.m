function mrc = marginalRiskContribution(wts,covMat)

w = wts(:);
portVar = w'*covMat*w;
portStd = sqrt(portVar);
if portStd < 1e-12
    warning('Portfolio standard deviation is near zero. Adding epsilon.');
    portStd = 1e-12;
end
disp(['Portfolio variance: ',num2str(w'*covMat*w)])
mrc = (w.*(covMat*w))/portStd;
end