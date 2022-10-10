%% Ashleigh Solano - rpCOMB project - August 10 ,2020
% This script will conduct rpcomb analysis by which an intensity frame
% stack will be processed into a brightness movie. Following some
% preliminary NB analysis this movie is split and a 2D-pCF function applied
% to all pixels in each oligomer channel.

%% Open.bin file 
[data,FileName,PathName]=Open_bin(32); % specify file dimensions

%% Detrend raw data 
[datad, datasm]=detrend_image(data,10,1); % only real data 

%% Pileup correction 
[Corrected]=Correctcounts(datad,20,125000); % only real data

%% Brightness movie parameters
n=100;
deln=1;     %  moving average del N=3 or 1 is acceptable 
sFactor = 1;
[B]=B_movie(data,n,deln); % specify corrected if want to calculate B with pileup
dataB=B/sFactor ;

%% Gaussian smoothing filter applied to bmovie 
smooth_dataB=imgaussfilt(dataB,0.7);

%% Brightness histogram
figure
histogram(smooth_dataB(:,:,1:100))
mean(smooth_dataB(:))

figure
imagesc(smooth_dataB(:,:,1))
colormap 'jet'
%% Average the intensity time series to match the B-movie
[intensity_av]=mav_intensity(data,n,deln);
close all
%% N&B to determine cursors for splitting 
% Download required of Robert Henson (2022). Flow Cytometry Data Reader and Visualization (https://www.mathworks.com/matlabcentral/fileexchange/8430-flow-cytometry-data-reader-and-visualization), MATLAB Central File Exchange. Retrieved April 20, 2022.
[NB_Bvals,NB_Ivals,x,y] = plot_IntensityNB(smooth_dataB,intensity_av); %plot this intensity image 
refreshdata 
plot_scatterNB(NB_Ivals,NB_Bvals);  % plot the scatter plot
ylim([0.8 3])
box on
%ensure to link the z values
%% Cursor box selection 
r1=rectangle('Position',[3 1.1 7 0.5],'EdgeColor',[0 0.5 0.5],'LineWidth',3);
r2=rectangle('Position',[3 1.61 7 0.39],'EdgeColor','g','LineWidth',3);
r3=rectangle('Position',[3 2.01 7 0.99],'EdgeColor','r','LineWidth',3);

%% Spatial map NB of all oligomers present  
[Monomer,Dimer,Oligomer]=Split_Bmovie(smooth_dataB,1.1,1.6,1.61,2.1,2.11,3);
[NB_map]=NB_smap(intensity_av,Monomer, Dimer, Oligomer,1);

%% Brightness movie filtering
[Monomer,Dimer,Oligomer]=Split_Bmovie(smooth_dataB,1.1,1.6,1.61,2.1,2.11,3);
%% Pair correlation profiles d = 2

[MPCF_2, angles, rbins] = pcomb2d(Monomer,2,10,[]);
[DPCF_2, angles, rbins] = pcomb2d(Dimer,2,10,[]);
[OPCF_2, angles, rbins] = pcomb2d(Oligomer,2,10,[]);
%% Pair correlation profiles d = 8 

[MPCF_8, angles, rbins] = pcomb2d(Monomer,8,10,[]);
[DPCF_8, angles, rbins] = pcomb2d(Dimer,8,10,[]);
[OPCF_8, angles, rbins] = pcomb2d(Oligomer,8,10,[]);


%% Plotting average pair correlation profiles for d=2 

[mpcf2_curve,mpcf2_ma]=av_pcf(MPCF_2);  % plot monomer d=2
[dpcf2_curve,dpcf2_ma]=av_pcf(DPCF_2);  % plot dimer = d=2
[opcf2_curve,opcf2_ma]=av_pcf(OPCF_2);  % plot oligomer d=2

figure
plot(mpcf2_ma,'color',[0 0.5 0.5],'Linewidth',2)
hold on
plot(dpcf2_ma,'g' ,'Linewidth',2)
hold on
plot(opcf2_ma,'r','Linewidth',2)
xlabel('Log time(s)');ylabel('G_B_(_t_,_d_r_=_2_)');
hold off

%% Plotting average pair correlation profiles for d=8 

frame_time = 0.105; % in seconds
b_sampling = frame_time.*delN; 
time_axisb = rbins .* b_sampling;

[mpcf8_curve,mpcf8_ma]=av_pcf(MPCF_8);  % plot monomer d=8
[dpcf8_curve,dpcf8_ma]=av_pcf(DPCF_8);  % plot dimer = d=8
[opcf8_curve,opcf8_ma]=av_pcf(OPCF_8);  % plot oligomer d=8

figure
plot(mpcf8_ma,'color',[0 0.5 0.5],'Linewidth',2)
hold on
plot(dpcf8_ma,'g','Linewidth',2)
hold on
plot(opcf8_ma,'r','Linewidth',2)
xlabel('Log time(s)');ylabel('G_B_(_t_,_d_r_=_8_)');
hold off 

%% Pair correlation profiles - moments connectivity map  

[MPCF_2, angles, rbins] = pcomb2d(Monomer,6,0,[]);
[DPCF_2, angles, rbins] = pcomb2d(Dimer,6,0,[]);
[OPCF_2, angles, rbins] = pcomb2d(Oligomer,6,0,[]);


%% Calculate Anisotropy values for all oligomer channels 
[Total_lAxisM,Total_sAxis,Total_Eccentricity,Total_AngleM,Total_AnisotropyM]=Moments(MPCF_2);
[Total_lAxisD,TotalM_sAxis,Total_Eccentricity,Total_AngleD,Total_AnisotropyD]=Moments(DPCF_2);
[Total_lAxisO,Total_sAxis,Total_Eccentricity,Total_AngleO,Total_AnisotropyO]=Moments(OPCF_2);

%% Plot connectivity map 
distance = 6; % specify pcf distance 

plot_Connectivity(Total_AnisotropyM,Total_AngleM,distance,[0 0.5 0.5]);
plot_Connectivity(Total_AnisotropyD,Total_AngleD,distance,[0 1 0]);
plot_Connectivity(Total_AnisotropyO,Total_AngleO,distance,[1 0 0]);
