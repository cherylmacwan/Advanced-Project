clear VARIABLES
clear plot_google_map
close all
clc
%% Load Terrain Data
load('SF_RasterData.mat')
figure; 
imagesc(terrain.data); 
%% Plot terrain data
axis xy; 
xlabel('x pixel'); 
ylabel('y pixel'); 
title('SF terrain')

figure; 
imagesc(buildings.buildingNums > 0); 
axis xy; 
xlabel('x pixel'); 
ylabel('y pixel'); 
title('SF Building Mask')

figure; imagesc(buildings.data - terrain.data); 
axis xy; 
xlabel('x pixel'); 
ylabel('y pixel'); 
title('SF Building AGL Height')

figure;
imagesc(buildings.data);
axis xy;
xlabel('x pixel');
ylabel('y pixel');
title('SF Building AMSL Height')
