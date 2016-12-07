clc
clear all
load('siteAndMeasurementsData_LessThan15m.mat')
% To plot the receiver positions and positions for recordings

for x=[1 6 85 111]
%     figure
%     plot(dataSet(x).measurements.AMSLHeight)
        for y=1:length(dataSet(x).measurements.AMSLHeight)
            figure
            plot(dataSet(x).measurements.lon,dataSet(x).measurements.lat,'b.','MarkerSize',0.5)
            hold on
            plot(dataSet(x).siteInfo.lon,dataSet(x).siteInfo.lat,'r.','MarkerSize',15)
            text(dataSet(x).siteInfo.lon,dataSet(x).siteInfo.lat,num2str(x))
            title(['Site -', num2str(x)])
%         z(y)=dataSet(x).measurements.AMSLHeight(y)-dataSet(x).siteInfo.AMSLHeight(1);
        end
%     A{1,x}=z;
%     B(1,x)=mean(A{1,x});
   
end
% B=transpose(B);
% stem(B)
% ylim([-2 2])
