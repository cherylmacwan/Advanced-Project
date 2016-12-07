clear all
clc
load('siteAndMeasurementsData_LessThan15m.mat')
for mesndx=113
%     len=length(dataSet(mesndx).measurements.pwr);
%     C=zeros(len,2);
%     figure(mesndx)
%     plot(dataSet(mesndx).measurements.pwr)
%     title(num2str(mesndx))
%     xlabel('Distance along route (m)  \rightarrow')
%     ylabel('Received Power (db) \rightarrow')
%     hold on
%     C(1:len,1)=(dataSet(mesndx).measurements.lon);
%     C(1:len,2)=(dataSet(mesndx).measurements.lat);
%     A{1,mesndx}=C;
    C(1,1)=dataSet(mesndx).siteInfo.lon;
    C(1,2)=dataSet(mesndx).siteInfo.lat;
    xlswrite('TP113_newestR',C,1,'B2');
    B={'Sr.No','Long','Lat'};
    xlswrite('TP113_newestR',B,1,'A1');
    D=1; %1:len;
    D=transpose(D);
    xlswrite('TP113_newestR',D,1,'A2');
    warning('off','MATLAB:xlswrite:AddSheet')
%     nn=1:len;
%     fitt=fit(transpose(1:len),dataSet(mesndx).measurements.pwr,'poly2');
%     plot(1:len,fitt.p1.*nn.^2+fitt.p2.*nn+fitt.p3,'k-')
%     hold off
end
   warning('off','MATLAB:xlswrite:AddSheet')
    