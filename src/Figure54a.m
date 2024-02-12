clc;
clear;
close all;

iterations = 20;
channel_realization_count = 10;

L = 400;
N = 1;
K = 40;

coherence_time = 0.002;
coherence_bandwidth = 100000;
tau_c = coherence_time*coherence_bandwidth;%coherence block length
tau_p = 10;%pilot length

BW = 20000000;

%Angular Standard Deviation - ASD - Nominal Values for azimuth and elevation as suggested in book
Azimuth_ASD = deg2rad(30);  %azimuth angle
Elevation_ASD = deg2rad(15);   %elevation angle

%% Propagation parameters

%Maxmimum uplink power for uplink as mentioned in Page 167
p = 100;

MMSE_Original = zeros(K,iterations); %MMSE (All)
MMSE_DCC = zeros(K,iterations); %MMSE (DCC)
MR_DCC = zeros(K,iterations); %MR (DCC)

for n = 1:iterations
    disp(['Iteration - ' num2str(n) '/' num2str(iterations)]);
    
    [gainOverNoisedB,R,pilot_stream,D,D_small] = generateSetup(L,K,N,tau_p,1,0,Azimuth_ASD,Elevation_ASD);
    
    [Hhat,H,B,C] = functionChannelEstimates(R,channel_realization_count,L,K,N,tau_p,pilot_stream,p);
    
    D_all = ones(L,K);
    
    [SE_MMSE_all] = functionMMSEScheme(Hhat,H,D_all,C,tau_c,tau_p,channel_realization_count,N,K,p);
    [SE_MMSE] = functionMMSEScheme(Hhat,H,D,C,tau_c,tau_p,channel_realization_count,N,K,p);
    [SE_MR]  = functionMRScheme(Hhat,H,D,B,C,tau_c,tau_p,channel_realization_count,N,K,L,p,R,pilot_stream);
    
    MMSE_Original(:,n) = SE_MMSE_all;
    MMSE_DCC(:,n) =  SE_MMSE;
    MR_DCC(:,n) =  SE_MR;
    
    clear Hhat H B C R;
    disp('Completed');
end


% Figure 5.4a
figure;

hold on; 
box on;
grid on;

plot(sort(MMSE_Original(:)),linspace(0,1,K*iterations),'g-','LineWidth',3);
plot(sort(MMSE_DCC(:)),linspace(0,1,K*iterations),'r-.','LineWidth',3);
plot(sort(MR_DCC(:)),linspace(0,1,K*iterations),'k-','LineWidth',3);

title('CDF of Spectral Efficiency');
xlabel('Spectral Efficiency','Interpreter','Latex');
legend({'MMSE','MMSE-DCC','MR-DCC'},'Interpreter','Latex','Location','NorthWest');
xlim([0 12]);

