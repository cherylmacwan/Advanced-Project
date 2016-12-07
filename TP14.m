clear VARIABLES
clear plot_google_map
close all
clc
load('siteAndMeasurementsData.mat')

%% Indexes
ndx=[311 312 313 292 294 295 296 297 298 259 242 241 233 223 213 198 186 185 168 155 143 126 125 141 140 139 138 137 136 154 166 184 197 212 222 231 232 239 254 269 288 307 306 305 304 303 302 301 284 268 253 230 220 210 193 181 161 162 163 164 148 149 150 151 152 122 108 102 95 86 77 73 70 69];
%Indexes for values along the path
site=14; %Site index among the data

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
hr=mean(dataSet(site).measurements.AMSLHeight(ndx)); %Average above ground height of Receivers
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

Sitetype(1:9)={'1T1'};
Sitetype(10:21)={'1T2'};
Sitetype(22:40)={'LOS'};
Sitetype(41:47)={'1T1'};
Sitetype(48:56)={'1T3'};
Sitetype(48:length(Power))={'LOS'};%Visually determined type of path for the given sets of points

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
        reflat=dataSet(site).measurements.lat(289);
        reflon=dataSet(site).measurements.lon(289); % 289 is selected since it is the point exactly at the cross-section of the 1-turn so the distance of transmitter to the point 191 will be R1 and from 191 to the 1-turn points will be R2 which keeps varying.
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
     elseif (strcmp(Sitetype(z), '1T2'));
        reflat=dataSet(site).measurements.lat(127);% The turn occurs around point with index 1125
        reflon=dataSet(site).measurements.lon(127);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        
    elseif (strcmp(Sitetype(z), '1T3'));
        reflat=dataSet(site).measurements.lat(160); % 804 is the center of the turn taken to reach the one turn path
        reflon=dataSet(site).measurements.lon(160);
        [arclen,~]=distance(Lat(z),Long(z),reflat,reflon);
        R2=arclen*R*pi/180;
        [arclen,~]=distance(Sitelat,Sitelon,reflat,reflon);
        R3=arclen*R*pi/180;
        R4=sqrt(R2^2+R3^2);
        simData(z)=10*log10(oneTurnPG(wavelength, ht, hr, R2, R3, 1, exact, polarization, eps));
        D(z)=(Power(z)-simData(z))/10;
        

    else 
        simData(z)=Power(z); 
      
    end
  
end

%% Calibration
s=1;
Calib=zeros(1,40-22+1); %Calibrating along LOS path
for m=22:40;
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
xlim([0 770])
hold on
plot(segLen,simData,'.-')
limP=ylim;
xlabel('Distance along route (m) \rightarrow')
ylabel('Power(dB) \rightarrow')
legend('Measured Power','Simulated Power')

hold on
for segndx=[9 21 40 47 56 73]
    plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
    text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
    hold on
end
hold off     


figure
plot(Long,Lat,'bo-')
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
text(600,15,['Variance=',num2str(er)])
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
