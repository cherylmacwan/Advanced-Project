clear VARIABLES
clear plot_google_map
close all
clc
load('siteAndMeasurementsData_LessThan15m.mat')

%% Indexes
ndx=[170 189 188 169 187 186 185 184 183 182 196 205 208 211 215 216 225 229 234 243 250 255 256 259 266 267 260 261 262 257 252 251 245 244 235 230 226 217 214 212 210 209 207 206 198 199 192 191 154 142 132 126 120 118 114 112 108 106 103 99 87 86 85 84 98 97 96 95 94 105 107 110 117 119 124 123 125 130 141 153 180 193 203 202 201 213 232 237 247 248 249 238 239 240 241 242 233];
%Indexes for values along the path
site=78; %Site index among the data

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
f=dataSet(site).siteInfo.FrequencyMHz; %Frequency of transmission
wavelength = c./(f*10e6); % wavelength
ht=dataSet(site).siteInfo.AGLHeight; %Above Ground height of Antenna
hr=9; %Average above ground height of Receivers
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

Sitetype(1:11)={'1T1'};
Sitetype(12:24)={'2T'};
Sitetype(25:28)={'1T'};
Sitetype(29:60)={'LOS'};
Sitetype(61:69)={'1T2'};
Sitetype(70:75)={'2T1a'};
Sitetype(76:81)={'2T1b'};
Sitetype(82:85)={'1T1'};
Sitetype(86:88)={'2T2'};
Sitetype(89:96)={'2T3'}; %Visually determined type of path for the given sets of points

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
        reflat=dataSet(site).measurements.lat(191);
        reflon=dataSet(site).measurements.lon(191); % 191 is selected since it is the point exactly at the cross-section of the 1-turn so the distance of transmitter to the point 191 will be R1 and from 191 to the 1-turn points will be R2 which keeps varying.
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
      
        
    elseif (strcmp(Sitetype(z), '1T'));
        reflat=dataSet(site).measurements.lat(262); % 262 is the center of the turn taken to reach the one turn path
        reflon=dataSet(site).measurements.lon(262);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
    elseif (strcmp(Sitetype(z), '2T'));
        reflat=dataSet(site).measurements.lat(191);
        reflon=dataSet(site).measurements.lon(191);% The first turn occurs at index 191 hence giving R1
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        reflat1=dataSet(site).measurements.lat(194);% Distance between 191 and 194 is R2
        reflon1=dataSet(site).measurements.lon(194);% 2nd turn occurs at index 194 hence distance from 194 to set of points gives R3
        [arclen,~]=distance(reflat1,reflon1,reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(reflat1,reflon1,Lat(z),Long(z));
        R1=arclen*R*pi/180;
        simData(z)=10*log10(twoTurnPG(wavelength, ht, hr, R1, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/20;
        
    elseif (strcmp(Sitetype(z), '1T2'));
        reflat=dataSet(site).measurements.lat(88);% Th eturn occurs around point with index 88
        reflon=dataSet(site).measurements.lon(88);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
        
    elseif (strcmp(Sitetype(z), '2T1a'));
        reflat=dataSet(site).measurements.lat(88);%First turn occurs at index 88 hence distance from transmitter to 88 is R1
        reflon=dataSet(site).measurements.lon(88);
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        reflat1=dataSet(site).measurements.lat(94);%2nd turn occurs at 94 so distance between 88 and 94 is R2
        reflon1=dataSet(site).measurements.lon(94);% Distane between points and 94 is R3
        [arclen,~]=distance(reflat1,reflon1,reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(reflat1,reflon1,Lat(z),Long(z));
        R1=arclen*R*pi/180;
        simData(z)=10*log10(twoTurnPG(wavelength, ht, hr, R1, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/20;
        
    elseif (strcmp(Sitetype(z), '2T1b'));
        reflat=dataSet(site).measurements.lat(171);% First turn occurs at 171 hence distance from transmitter and 171 is R1
        reflon=dataSet(site).measurements.lon(171);
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        reflat1=dataSet(site).measurements.lat(194);%2nd turn occurs at 194 hence distance between 171 and 194 is R2
        reflon1=dataSet(site).measurements.lon(194);% Distance between each point and 194 is R3
        [arclen,~]=distance(reflat1,reflon1,reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(reflat1,reflon1,Lat(z),Long(z));
        R1=arclen*R*pi/180;
        simData(z)=10*log10(twoTurnPG(wavelength, ht, hr, R1, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/20;
        
    elseif (strcmp(Sitetype(z), '2T2'));
        reflat=dataSet(site).measurements.lat(191);% First turn occurs at 191 hence distance from transmitter to 191 is R1
        reflon=dataSet(site).measurements.lon(191);
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        reflat1=dataSet(site).measurements.lat(200);% 2-nd turn occurs at 200 hence distance between 191 and 200 is R2
        reflon1=dataSet(site).measurements.lon(200);% Distance between each point and 200 is R3
        [arclen,~]=distance(reflat1,reflon1,reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(reflat1,reflon1,Lat(z),Long(z));
        R1=arclen*R*pi/180;
        simData(z)=10*log10(twoTurnPG(wavelength, ht, hr, R1, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/20;
        
     elseif (strcmp(Sitetype(z), '2T3'));
        reflat=dataSet(site).measurements.lat(194);% First turn occurs at 194 hence distacne from transmitter to 194 is R1
        reflon=dataSet(site).measurements.lon(194);
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        reflat1=dataSet(site).measurements.lat(233);% 2-nd turn occurs at 233 hence distance between 233 and 194 is R2
        reflon1=dataSet(site).measurements.lon(233);% Distance between each point and 233 is R3
        [arclen,~]=distance(reflat1,reflon1,reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(reflat1,reflon1,Lat(z),Long(z));
        R1=arclen*R*pi/180;
        simData(z)=10*log10(twoTurnPG(wavelength, ht, hr, R1, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/20;
        
    else 
        simData(z)=Power(z); 
      
      
    end
  
end

%% Calibration
s=1;
Calib=zeros(1,60-29+1); %Calibrating along LOS path
for m=29:60;
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

plot(a,Power,'.-')
% xlim([0 1160])
hold on
plot(a,simData,'.-')
limP=ylim;
xlabel('Distance along route (m) \rightarrow')
ylabel('Power(dB) \rightarrow')
legend('Measured Power','Simulated Power')

% hold on
% for segndx=[11 24 28 60 69 75 81 85 88 96]
%     plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
%     text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
%     hold on
% end
% hold off     


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
text(450,mean(error)+1,['Mean=',num2str(mean(error))])
text(850,10,['Variance=',num2str(er)])
hold off

figure
plot(segLen,D)
hold on
limP=ylim;
for segndx=[11 24 28 60 69 75 81 85 88 96]
    plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
    text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
    hold on
end
hold off
title('Calculated value of S')
xlabel('Distance along route \rightarrow')
ylabel('S (dB) \rightarrow')
