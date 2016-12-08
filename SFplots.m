clc
clear all
close all
load('siteAndMeasurementsData.mat')
% To plot the receiver positions and positions for recordings

for x=[1 5 7 14]
    figure
    subplot(2,1,1)
    plot(dataSet(x).measurements.AMSLHeight)
%     for y=1:length(dataSet(x).measurements.AMSLHeight)
        subplot(2,1,2)
        plot(dataSet(x).measurements.lon,dataSet(x).measurements.lat,'b.','MarkerSize',0.5)
        hold on
        plot(dataSet(x).siteInfo.lon,dataSet(x).siteInfo.lat,'r.','MarkerSize',15)
        text(dataSet(x).siteInfo.lon,dataSet(x).siteInfo.lat,num2str(x))
        %             title(['Site -', num2str(x)])
        hold on
        plot_google_map
        %              z(y)=dataSet(x).measurements.AMSLHeight(y)-dataSet(x).siteInfo.AMSLHeight(1);
%     end
    %     A{1,x}=z;
    %     B(1,x)=mean(A{1,x});
    %
end

% B=transpose(B);
% stem(B)
% ylim([-2 2])
