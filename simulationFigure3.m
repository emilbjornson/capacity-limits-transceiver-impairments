%This Matlab script can be used to generate Figure 3 in the paper:
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

Nt = [4 12]; %Number of transmit antennas (two different values)
Nr = 4; %Number of receive antennas
M = Nr; %Minimum of Nt and Nr (same for all Nt in this simulation)

%Generate channel fading realizations
H = (randn(Nr,max(Nt),nbrOfRealizations)+1i*randn(Nr,max(Nt),nbrOfRealizations))/sqrt(2);

%Range of SNRs in simulation
SNRdB = -10:1:50; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale

%Level of impairments in transmitter hardware
kappa = 0.05;

%Number of iterations in the search algorithm that finds the optimal power
%allocation in case of 
nbrOfPowerIterations = 10;


%Placeholders for storing the capacity results
capacityImpairments_determ_alpha0 = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Deterministic MIMO channels for alpha=0 (which will be averaged)
capacityImpairments_determ_alpha1 = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Deterministic MIMO channels for alpha=1 (which will be averaged)
capacityImpairments_random = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Random uncorrelated Rayleigh fading MIMO channels


%Go through all channel realizations
for k = 1:nbrOfRealizations
    
    %Go through all numbers of transmit antennas
    for l = 1:length(Nt)
        
        Ntcurrent = Nt(l); %Current number of transmit antennas
        
        %Compute singular value decomposition of current channel
        %realization
        [U,D,V] = svd(H(:,1:Ntcurrent,k));
        
        %Squared singular values of current channel realization 
        %(this is the same as the eigenvalues of the Gram matrix)
        eigenvalues = diag(D).^2;
        
        %Inversion of the squared singular values, which is used in the
        %waterfilling algorithm for power allocation
        waterfillingFactorsIdealHardware = 1./(eigenvalues + 1e-10);
        
        
        %Go through all SNR values
        for m = 1:length(SNR)
            
            %Compute capacity-achieving power allocation for deterministic
            %channels. Corollary 2 shows that it is also the
            %capacity-achieving (under transceiver hardware impairments with
            %alpha=1) is given by the waterfilling algorithm from [1].
            powerAllocationDeterministicAlpha1 = functionWaterfilling(SNR(m),waterfillingFactorsIdealHardware);
            
            %Compute capacity with transceiver hardware impairments for alpha=1
            %and a deterministic channel realization. It is based on 
            %Eq. (6) and uses the power allocation from above.
            capacityImpairments_determ_alpha1(m,k,l) = sum(log2(1+powerAllocationDeterministicAlpha1.*eigenvalues ./ ( kappa^2/Ntcurrent*SNR(m)*eigenvalues+1)));
            
            %Compute capacity with transceiver hardware impairments and one
            %realization of uncorrelated Rayleigh fading. Equal power
            %allocation is optimal, according to Corollary 1. The averaging
            %is performed when plotting the results.
            capacityImpairments_random(m,k,l) = real(log2(det(eye(Nr)+(SNR(m)/Ntcurrent+SNR(m)*kappa^2/Ntcurrent)*H(:,1:Ntcurrent,k)*H(:,1:Ntcurrent,k)'))-log2(det(eye(Nr)+(SNR(m)*kappa^2/Ntcurrent)*H(:,1:Ntcurrent,k)*H(:,1:Ntcurrent,k)')));

            %Compute capacity with transceiver hardware impairments for alpha=0
            %and a deterministic channel realization. The
            %capacity-achieving power allocation is not known in this case.
            %It is computed iteratively by 1) Fixing Upsilont in Lemma 1;
            %2) Computing the rate-maximizing signal covariance matrix Q; 
            %3) Updating Upsilont based on the Q-matrix; 4) Iterate these
            %steps. This algorithm converges to the local optimum, which
            %appears to also be global optimum in this case.
            Q = V(:,1:Nr)*diag(powerAllocationDeterministicAlpha1)*V(:,1:Nr)';
            Upsilont = kappa^2*diag(diag(Q)); %Compute Upsilon_t as defined in Eq. (3)
            
            for j = 1:nbrOfPowerIterations %Update the power allocation iteratively
                effectiveChannel = sqrtm(eye(Nr)+H(:,1:Ntcurrent,k)*Upsilont*H(:,1:Ntcurrent,k)')*H(:,1:Ntcurrent,k); %Compute the effective (normalized) channel for a given Upsilon_t
                [U,D,V] = svd(effectiveChannel); %Compute singular value decomposition of effective channel
                waterfillingFactorsIteration = 1./(diag(D).^2 + 1e-10); %Inversion of the squared singular values, used in waterfilling algorithm
                powerAllocationDeterministicIteration = functionWaterfilling(SNR(m),waterfillingFactorsIteration); %New waterfilling
                Q = V(:,1:Nr)*diag(powerAllocationDeterministicIteration)*V(:,1:Nr)'; %Compute the new signal covariance matrix
                Upsilont = kappa^2*diag(diag(Q)); %Compute Upsilon_t as defined in Eq. (3)
            end
            
            %Compute capacity with transceiver hardware impairments for alpha=0
            capacityImpairments_determ_alpha0(m,k,l) = real(log2(det(eye(Nr)+H(:,1:Ntcurrent,k)*(Q+Upsilont)*H(:,1:Ntcurrent,k)'))-log2(det(eye(Nr)+H(:,1:Ntcurrent,k)*Upsilont*H(:,1:Ntcurrent,k)')));
            
        end
        
    end
end


%Compute high-SNR limits from Eq. (5)
capacityLimits = M*log2(1+Nt/(M*kappa^2)); 


%Plot Figure 3 from the paper
figure; hold on; box on;

plot(SNRdB,mean(capacityImpairments_determ_alpha1(:,:,1),2),'k');
plot(SNRdB,mean(capacityImpairments_determ_alpha0(:,:,1),2),'r--');
plot(SNRdB,mean(capacityImpairments_random(:,:,1),2),'b-.');

plot(SNRdB,mean(capacityImpairments_determ_alpha1(:,:,2),2),'k');
plot(SNRdB,mean(capacityImpairments_determ_alpha0(:,:,2),2),'r--');
plot(SNRdB,mean(capacityImpairments_random(:,:,2),2),'b-.');

plot(SNRdB,capacityLimits(1)*ones(size(SNRdB)),'k:');
plot(SNRdB,capacityLimits(2)*ones(size(SNRdB)),'k:');

text(SNRdB(41),mean(capacityImpairments_determ_alpha1(41,:,1),2)-1.5,['N_t = ' num2str(Nt(1))]);
text(SNRdB(41),mean(capacityImpairments_determ_alpha1(41,:,2),2)-1.5,['N_t = ' num2str(Nt(2))]);

xlabel('SNR [dB]')
ylabel('Capacity [bit/s/Hz]')
axis([SNRdB(1) SNRdB(end) 0 50]);
legend('Deterministic, \alpha=1 (Average)','Deterministic, \alpha=0 (Average)','Uncorr Rayleigh Fading, any \alpha','Location','SouthEast')
