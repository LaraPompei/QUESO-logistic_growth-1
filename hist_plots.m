clear all;
a_data
model_data
% Histogram plots ---------------------------------------------------------
% a_data
nbins=100;
hist(a_data,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');%
title("Histogram Values for 'a' Data",'fontname', 'Times', 'fontsize',20);
xlabel('','fontname', 'Times', 'fontsize',20);
ylabel('','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng histogram_Values_a_10000.png
pause;
clf;

%model_data
hist(model_data,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');%
title("Histogram Values for 'y(a,3)' for 'a' Data",'fontname', 'Times', 'fontsize',20);
xlabel('','fontname', 'Times', 'fontsize',20);
ylabel('','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng histogram_Values_model_10000.png
pause;
clf;
