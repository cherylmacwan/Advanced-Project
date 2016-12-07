clear VARIABLES
clear plot_google_map
close all
clc
load('siteAndMeasurementsData.mat')

%% Indexes
ndx=[804 830 802 829 828 827 826 825 824 851 823 850 848 847 874 846 873 845 872 871 870 869 900 868 867 920 942 964 978 993 1004 1023 1039 1056 1097 1120 1121 1142 1141 1167 1168 1143 1144 1169 1145 1146 1147 1148 1149 1150 1126 1127 1128 1129 1130 1104 1105 1106 1107 1078 1043 1026 1011 995 981 967 946 922 903 877 832 787 767 747 725 706 685 684 660 639 617 597 577 555 518 517 516 515 514 513 537 536 535 554 576 596 616 638 659 658 657 656 681 680 679 678 677 702 722 745 763 785 799 820 866 ];
%Indexes for values along the path
site=7; %Site index among the data

%% Pre-allocation for speed
m=zeros(1,(length(ndx)-1));
n=zeros(1,length(ndx));
D=m;
k=n;
Lat=m;
Long=m;
Lat1=m;
Long1=m;
Power=m;
segLen=m;
seglength=m;
simData=m;

%% Initialize values
a=1:(length(ndx)-1); %Due to segLength the length of Power values will be 1 less than all the index values
b=1;
R=6.371e6; % Radius of earth to calculate distance using Latitude and Longitude values
c = 3*10^8; % speed of light
f=dataSet(site).antennaInfo.MaxFreq; %Frequency of transmission
wavelength = c./(f*10e6); % wavelength
ht=dataSet(site).siteInfo.AGLHeight; %Above Ground height of Antenna
hr=mean(dataSet(site).measurements.AMSLHeight(ndx))-6; %Average above ground height of Receivers
eps = 15 -1i*0.1;
polarization = 0;
exact = 1;

for b=1:(length(ndx)-1)
    k=ndx;
    Power(b)=dataSet(site).measurements.pwr(ndx(b));
    Lat(b)=dataSet(site).measurements.lat(ndx(b));
    Long(b)=dataSet(site).measurements.lon(ndx(b));
    Lat1(b)=dataSet(site).measurements.lat(ndx(b+1));
    Long1(b)=dataSet(site).measurements.lon(ndx(b+1));
    %Lat(a+1)=dataSet(site).measurements.lat(ndx(a+1));
    %Long(a+1)=dataSet(site).measurements.lon(ndx(a+1));
    [arclen,~] = distance(Lat(b),Long(b),Lat1(b),Long1(b));
    seglength(b)=arclen*(R/180)*pi; %Distance between Receiver location and next receiver location or distance between two consecutive points along a path
    %a=a+1;
end

Sitetype(1:25)={'LOS'};
Sitetype(26:38)={'1T1'};
Sitetype(39:58)={'1T2'};
Sitetype(59:83)={'1T3'};
Sitetype(84:93)={'1T4'};
Sitetype(94:99)={'LOS'};
Sitetype(100:107)={'1T5'};
Sitetype(108:114)={'1T6'};%Visually determined type of path for the given sets of points

Sitelat=dataSet(site).siteInfo.lat;
Sitelon=dataSet(site).siteInfo.lon; %Transmitter position

segLen(1)=0;
for w=2:length(ndx)-1;
    segLen(w)=segLen(w-1)+seglength(w-1);
end %Creating distance vector for each point denoting distance from the start


%% Prediction using 2-ray model

for z=1:length(Power)
    if (strcmp(Sitetype(z), 'LOS'))
        [arclen,~]=distance(Lat(z),Long(z),Sitelat,Sitelon);
        R1=arclen*R*pi/180;
        simData(z)=10*log10(exact2RayModel(ht,hr,R1, polarization, eps, wavelength));
        D(z)=0;
        
        
    
    elseif (strcmp(Sitetype(z), '1T1'));
        reflat=dataSet(site).measurements.lat(897);
        reflon=dataSet(site).measurements.lon(897); % 897 is selected since it is the point exactly at the cross-section of the 1-turn so the distance of transmitter to the point 191 will be R1 and from 191 to the 1-turn points will be R2 which keeps varying.
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
      
        
    elseif (strcmp(Sitetype(z), '1T3'));
        reflat=dataSet(site).measurements.lat(804); % 804 is the center of the turn taken to reach the one turn path
        reflon=dataSet(site).measurements.lon(804);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
     elseif (strcmp(Sitetype(z), '1T5'));
        reflat=dataSet(site).measurements.lat(659); % 659 is the center of the turn taken to reach the one turn path
        reflon=dataSet(site).measurements.lon(659);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
       
    elseif (strcmp(Sitetype(z), '1T2'));
        reflat=dataSet(site).measurements.lat(1125);% The turn occurs around point with index 1125
        reflon=dataSet(site).measurements.lon(1125);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
    elseif (strcmp(Sitetype(z), '1T4'));
        reflat=dataSet(site).measurements.lat(553);% The turn occurs around point with index 1125
        reflon=dataSet(site).measurements.lon(553);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
    elseif (strcmp(Sitetype(z), '1T6'));
        reflat=dataSet(site).measurements.lat(897);% The turn occurs around point with index 1125
        reflon=dataSet(site).measurements.lon(897);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
         
    else 
        simData(z)=Power(z); 
%       
%       
    end
%   
end

%% Calibration
s=1;
Calib=zeros(1,25-1+1); %Calibrating along LOS path
for m=1:25;
    Calib(s)=Power(m)-simData(m);
    s=s+1;
end
m=mean(Calib);

simData=simData+m;

%% Error Fit curve

error=Power-simData;
nn=segLen;
ff=db(std(error));
fitt=fit(transpose(nn),transpose(error),'poly2');
er=var(error);

%% Finding the value of S
% D=10.^D;


%% Plot Measured and Simulated Power and Error Values

plot(segLen,Power,'.-')
xlim([0 1170])
hold on
plot(segLen,simData,'.-')
limP=ylim;
xlabel('Distance along route (m) \rightarrow')
ylabel('Power(dB) \rightarrow')
legend('Measured Power','Simulated Power')

hold on
for segndx=[25 38 58 83 93 99 107 114]
    plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
    text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
    hold on
end
hold off     


figure
plot(Long,Lat,'b.-')
plot_google_map
hold on
plot(dataSet(site).siteInfo.lon,dataSet(site).siteInfo.lat,'r.','MarkerSize',15)
plot_google_map
hold off

figure
plot(segLen,error)
hold on
plot(nn,fitt.p1.*nn.^2+fitt.p2.*nn+fitt.p3)
hold on
limPG=xlim;
plot([limPG(1), limPG(2)],[mean(error),mean(error)],'--m')
text(600,mean(error)+1,['Mean=',num2str(mean(error))])
text(850,15,['Variance=',num2str(er)])
hold off

% figure
% plot(segLen,D)
% hold on
% limP=ylim;
% for segndx=[11 24 28 60 69 75 81 85 88 96]
%     plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
%     text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
%     hold on
% end
% hold off
% title('Calculated value of S')
% xlabel('Distance along route \rightarrow')
% ylabel('S (dB) \rightarrow')
