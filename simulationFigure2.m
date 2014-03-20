%This Matlab script can be used to generate Figure 2 in the paper:
%
%Emil Björnson, Per Zetterberg, Mats Bengtsson, Björn Ottersten, "Capacity
%Limits and Multiplexing Gains of MIMO Channels with Transceiver
%Impairments," IEEE Communications Letters, vol. 17, no. 1, pp. 91-94,
%January 2013.
%
%The paper is available for free download: http://arxiv.org/pdf/1209.4093
%
%This is version 1.0.
%
%License: If you use parts of this code in your research, please cite the
%paper above in your resulting research publications.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

%Initialization
close all;
clear all;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Number of realizations in Monte-Carlo simulations (this number must be
%very large to get a perfect agreement in the asymptotic results, which
%rely on the ratio between two capacities).
nbrOfRealizations = 10000;

Nt = 4; %Number of transmit antennas
Nr = 4; %Number of receive antennas
M = min([Nt,Nr]); %Minimum of Nt and Nr

%Range of SNRs in simulation
SNRdB = -10:1:70; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale

%Level of impairments in transmitter hardware
kappa = [0.05 0.1];


%Part 1: Synthetic channels

%Generate channel realizations based on uncorrelated Rayleigh fading
H = (randn(Nr,max(Nt),nbrOfRealizations)+1i*randn(Nr,max(Nt),nbrOfRealizations))/sqrt(2);

%Placeholders for storing the capacity results
capacityIdealHardware_synthetic = zeros(length(SNR),nbrOfRealizations);
capacityImpairments_synthetic = zeros(length(SNR),nbrOfRealizations,length(kappa));


%Go through all channel realizations
for k = 1:nbrOfRealizations
    
    %Squared singular values of current channel realization
    %(this is the same as the eigenvalues of the Gram matrix)
    eigenvalues = svd(H(:,:,k)).^2;
    
    %Inversion of the squared singular values, which is used in the
    %waterfilling algorithm for power allocation
    waterfillingFactorsIdealHardware = 1./(eigenvalues+1e-10);
    
    %Go through all SNR values
    for m = 1:length(SNR)
        
        %Compute capacity-achieving power allocation for deterministic
        %channels, using the waterfilling algorithm from [1]. Corollary
        %2 shows that it is also the capacity-achieving under
        %transceiver hardware impairments with alpha=1.
        powerAllocationDeterministic = functionWaterfilling(SNR(m),waterfillingFactorsIdealHardware);
        
        %Compute capacity with ideal hardware and a deterministic channel
        %realization, using the power allocation from above.
        capacityIdealHardware_synthetic(m,k) = sum(log2(1+powerAllocationDeterministic./waterfillingFactorsIdealHardware));
        
        %Go through different levels of hardware impairments
        for l = 1:length(kappa)
            
            %Compute capacity with transceiver hardware impairments for alpha=1
            %and a deterministic channel realization. It is based on
            %Eq. (6) and uses the power allocation from above.
            capacityImpairments_synthetic(m,k,l) = sum(log2(1+powerAllocationDeterministic.*eigenvalues ./ ( kappa(l)^2/Nt*SNR(m)*eigenvalues+1)));
        end
        
    end
    
end



%Part 2: Synthetic channels

%Load file with measured channels
load Hmeasured; %Gives the matrix Hmeasured of dimensions 4 x 4 x 1773

nbrOfMeasurements = size(Hmeasured,3); %Extract number of measurements (1773)

%Placeholders for storing the capacity results
capacityIdealHardware_measured = zeros(length(SNR),nbrOfMeasurements);
capacityImpairments_measured = zeros(length(SNR),nbrOfMeasurements,length(kappa));


%Go through all channel measurements
for k = 1:nbrOfMeasurements
    
    %Squared singular values of current channel realization
    %(this is the same as the eigenvalues of the Gram matrix)
    eigenvalues = svd(Hmeasured(:,:,k)).^2;
    
    %Inversion of the squared singular values, which is used in the
    %waterfilling algorithm for power allocation
    waterfillingFactorsIdealHardware = 1./(eigenvalues+1e-10);

    %Go through all SNR values
    for m=1:length(SNR)
        
        %Compute capacity-achieving power allocation for deterministic
        %channels, using the waterfilling algorithm from [1]. Corollary
        %2 shows that it is also the capacity-achieving under
        %transceiver hardware impairments with alpha=1.
        powerAllocationDeterministic = functionWaterfilling(SNR(m),waterfillingFactorsIdealHardware);
        
        %Compute capacity with ideal hardware and a deterministic channel
        %realization, using the power allocation from above.
        capacityIdealHardware_measured(m,k) = sum(log2(1+powerAllocationDeterministic./waterfillingFactorsIdealHardware));
        
        %Go through different levels of hardware impairments
        for l = 1:length(kappa)
            
            %Compute capacity with transceiver hardware impairments for alpha=1
            %and a deterministic channel realization. It is based on
            %Eq. (6) and uses the power allocation from above.
            capacityImpairments_measured(m,k,l) = sum(log2(1+powerAllocationDeterministic.*eigenvalues ./ ( kappa(l)^2/Nt*SNR(m)*eigenvalues+1)));
        end

    end
    
end



%Compute high-SNR limits from Eq. (5). Note that it is the same
%irrespective of the channel distributions
capacityLimits = M*log2(1+Nt./(M*kappa.^2)); 


%Plot Figure 2 from the paper
figure; hold on; box on;

plot(SNRdB,mean(capacityIdealHardware_synthetic,2),'k');
plot(SNRdB,mean(capacityImpairments_synthetic(:,:,1),2),'k--');
plot(SNRdB,mean(capacityImpairments_synthetic(:,:,2),2),'k-.');

plot(SNRdB,mean(capacityIdealHardware_measured,2),'r');
plot(SNRdB,mean(capacityImpairments_measured(:,:,1),2),'r--');
plot(SNRdB,mean(capacityImpairments_measured(:,:,2),2),'r-.');

plot(SNRdB,capacityLimits(1)*ones(size(SNRdB)),'k:');
plot(SNRdB,capacityLimits(2)*ones(size(SNRdB)),'k:');

text(-8,30,'Capacity limits for different \kappa');
text(-8,15,'Synthetic channels');
text(20,10,'Measured channels');

xlabel('SNR [dB]')
ylabel('Average Capacity [bit/s/Hz]')
axis([SNRdB(1) SNRdB(end) 0 70]);
legend('Ideal Transceiver Hardware','Transceiver Impairments: \kappa=0.05','Transceiver Impairments: \kappa=0.1','Location','NorthWest')
