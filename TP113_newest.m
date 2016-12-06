clear VARIABLES
clear plot_google_map
clc
close all
load('siteAndMeasurementsData_LessThan15m.mat')

%% Indexes
%ndx=[609 610 611 612 613 593 594 595 596 597 575 578 579 580 581 598 616 630 660 671 705 793 791 790 809 808 807 788 787 786 768 754 743 730 729 718 717 711 702 701 694 679 686 670 669 668 663 659 658 651 646 645 636 625 626 627 608 ];
%ndx=[605 626 627 609 610 611 613 614 592 594 595 596 575 576 577 578 579 580 561 562 564 565 543 544 545 546 548 522 523 524 525 506 507 509 510 511 489 490 491 492 493 474 475 476 454 477 455 457 435 ]; % Slim width straight path 
%ndx=[519 518 502 485 484 467 466 431 421 420 406 419 404 392 391 369 360 337 321 314 313 302 295 285 275 268 262 257 251 250 248 237 230 229 224 218 217 204 193 182 176 170 164 163 158 151 144 143 142 131 126 116 112 ]; % Wider road straight path top half
 ndx=[ 570 553 532 516 515 499 498 482 481 464 463 446 445 429 418 417 403 390 377 359 347 335 328 327 320 312 300 294 293 284 273 266 260 256 246 246 240 233 227 220 213 200 190 179 173 167 149 148 139 135 133 129 119];
site=113;

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
a=1:(length(ndx)-1);
b=1;
R=6.371e6;
c = 3*10^8; % speed of light
f=dataSet(site).siteInfo.FrequencyMHz;
wavelength = c./(f*10e6); % wavelength
ht=dataSet(site).siteInfo.AGLHeight;
hr=2.077;
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
    seglength(b)=arclen*(R/180)*pi;
    %a=a+1;
end
segLen(1)=seglength(1);
for w=2:length(ndx)-1;
    segLen(w)=segLen(w-1)+seglength(w);
end


% Sitetype(1:56)={'NA'};
% Sitetype(31:53)={'LOS'};
Sitelat=dataSet(site).siteInfo.lat;
Sitelon=dataSet(site).siteInfo.lon;

for z=1:length(ndx)-1;
%     if (strcmp(Sitetype(z), 'LOS'))
            [arclen,~]=distance(Lat(z),Long(z),Sitelat,Sitelon);
            R1=arclen*R*pi/180;
            simData(z)=10*log10(exact2RayModel(ht,hr,R1, polarization, eps, wavelength));
            D(z)=0;
    
%     else
%         simData(z)=-94;
%     end
end

%% Calibration
s=1;
Calib=zeros(1,53-31+1);
for m=1:length(ndx)-1;
    Calib(s)=Power(m)-simData(m);
    s=s+1;
end
m=mean(Calib);

simData=simData+m;
%% Plotting

plot(segLen,Power,'.-')
hold on
limX=xlim;
plot([limX(1) limX(2)],[mean(Power) mean(Power)],'-')
text(limX(2)/2,mean(Power)+1,['Mean=',num2str(mean(Power))])
hold on
plot(segLen,simData,'.-')
hold on
% limXB=xlim;
% plot([limXB(1) limXB(2)],[mean(simData) mean(simData)],'-')
% text(limXB(2)/2,mean(simData)+5,['Mean=',num2str(mean(simData))])
limP=ylim;
xlabel('Distance along route (m) \rightarrow')
ylabel('Power(dB) \rightarrow')
legend('Measured Power','Mean','Simulated Power')
% hold on
% for segndx=[31 53]
%     plot([sum(seglength(1:segndx)),sum(seglength(1:segndx))],[limP(1), limP(2)],'--m')
%     text(sum(seglength(1:segndx))+10, limP(2)-10, Sitetype(segndx),'FontSize',10,'Rotation',270)
%     hold on
% end
% hold off

figure
plot(Long,Lat,'b.-')
hold on
plot(dataSet(site).siteInfo.lon,dataSet(site).siteInfo.lat,'r.','MarkerSize',15)
plot_google_map
hold off
% figure
% for i=31:53
%     plot(Long(i),Lat(i),'b.')
%     hold on
% end
% plot_google_map