clear all;
cd outputData/

sipOutput_raw_chain
sipOutput_filt_chain

% Priors ------------------------------------------------------------------
% Parameter K
minPrior0 = 0.0;
maxPrior0 = 5.0;
numPrior0 = 200;
numHorizPts0 = 20;
deltaPrior0 = (maxPrior0 - minPrior0)/numPrior0;
ip_prior_0_grid   = [minPrior0-numHorizPts0*deltaPrior0 : deltaPrior0 : maxPrior0+numHorizPts0*deltaPrior0];
ip_prior_0_values = ones(1,numPrior0+2*numHorizPts0+1)./(maxPrior0-minPrior0);
ip_prior_0_values(1,1:numHorizPts0)     = 0.;
ip_prior_0_values(1,numPrior0+numHorizPts0+2:numPrior0+2*numHorizPts0+1) = 0.;
plot(ip_prior_0_grid,ip_prior_0_values,'b*');
hold

left_vertical_line_x0 = ones(1,11)*minPrior0;
left_vertical_line_y0 = [0 : ip_prior_0_values(1,numHorizPts0+1)/10 : ip_prior_0_values(1,numHorizPts0+1)];
plot(left_vertical_line_x0,left_vertical_line_y0,'b--','linewidth',1);

right_vertical_line_x0 = ones(1,11)*maxPrior0;
right_vertical_line_y0 = [0 : ip_prior_0_values(1,numHorizPts0+1)/10 : ip_prior_0_values(1,numHorizPts0+1)];
plot(right_vertical_line_x0,right_vertical_line_y0,'b--','linewidth',1);

ylabel('Prior marginal PFD','fontsize',20);
xlabel('K (min^{-1})','fontsize',20);
title('Parameter K: prior marginal PFD','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng cal_parameterK_prior.png
waitforbuttonpress;
clf

% Parameter M
minPrior1 = 10.0;
maxPrior1 = 20.0;
numPrior1 = 200;
numHorizPts1 = 20;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
ip_prior_1_grid   = [minPrior1-numHorizPts1*deltaPrior1 : deltaPrior1 : maxPrior1+numHorizPts1*deltaPrior1];
ip_prior_1_values = ones(1,numPrior1+2*numHorizPts1+1)./(maxPrior1-minPrior1);
ip_prior_1_values(1,1:numHorizPts1)     = 0.;
ip_prior_1_values(1,numPrior1+numHorizPts1+2:numPrior1+2*numHorizPts1+1) = 0.;
plot(ip_prior_1_grid,ip_prior_1_values,'b*');
hold

left_vertical_line_x1 = ones(1,11)*minPrior1;
left_vertical_line_y1 = [0 : ip_prior_1_values(1,numHorizPts1+1)/10 : ip_prior_1_values(1,numHorizPts1+1)];
plot(left_vertical_line_x1,left_vertical_line_y1,'b--','linewidth',1);

right_vertical_line_x1 = ones(1,11)*maxPrior1;
right_vertical_line_y1 = [0 : ip_prior_1_values(1,numHorizPts1+1)/10 : ip_prior_1_values(1,numHorizPts1+1)];
plot(right_vertical_line_x1,right_vertical_line_y1,'b--','linewidth',1);

ylabel('Prior marginal PFD','fontsize',20);
xlabel('M (J/mol)','fontsize',20);
title('Parameter M: prior marginal PFD','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng cal_parameterM_prior.png
waitforbuttonpress;
clf

% Parameter 2 - M ------------
% KDM plots
[f1,x1] = ksdensity(ip_mh_rawChain_unified(2,:),'function','pdf');
[f2,x2] = ksdensity(ip_mh_rawChain_unified(2,:),'function','pdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('Posterior marginal PFDs of parameter M','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng cal_val_parameterM_PDF.png
waitforbuttonpress;
clf;


% CDF plots
[f1,x1] = ksdensity(ip_mh_rawChain_unified(2,:),'function','cdf');
[f2,x2] = ksdensity(ip_mh_rawChain_unified(2,:),'function','cdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('CFDs of parameter M','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng cal_val_parameterM_CDF.png
waitforbuttonpress;

