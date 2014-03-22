%This Matlab script can be used to generate Figure 4 and Figure 5 in the paper:
%
%Emil Björnson, Per Zetterberg, Mats Bengtsson, Björn Ottersten, "Capacity
%Limits and Multiplexing Gains of MIMO Channels with Transceiver
%Impairments," IEEE Communications Letters, vol. 17, no. 1, pp. 91-94,
%January 2013.
%
%The paper is available for free download: http://arxiv.org/pdf/1209.4093
%
%This is version 1.0 (Last edited: 2014-03-22)
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

Nt = [4 8 12]; %Number of transmit antennas in MIMO case (three different values)
Nr = 4; %Number of receive antennas in MIMO case
M = Nr; %Minimum of Nt and Nr (same for all Nt in this simulation)

%Generate channel fading realizations
H = (randn(Nr,max(Nt),nbrOfRealizations)+1i*randn(Nr,max(Nt),nbrOfRealizations))/sqrt(2); %MIMO case
h = (randn(nbrOfRealizations,1)+1i*randn(nbrOfRealizations,1))/sqrt(2); %SISO case

%Range of SNRs in simulation
SNRdB = -20:1:70; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale

%Level of impairments in transmitter hardware
kappa = 0.05;


%Placeholders for storing the capacity results: Ideal hardware
MIMOcapacityIdealHardware_deterministic = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Deterministic MIMO channels (which will be averaged)
MIMOcapacityIdealHardware_random = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Random uncorrelated Rayleigh fading MIMO channels
SISOcapacityIdealHardware = zeros(length(SNR),nbrOfRealizations); %SISO channels

%Placeholders for storing the capacity results: Ideal hardware
MIMOcapacityImpairments_deterministic = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Deterministic MIMO channels (which will be averaged)
MIMOcapacityImpairments_random = zeros(length(SNR),nbrOfRealizations,length(Nt)); %Random uncorrelated Rayleigh fading MIMO channels
SISOcapacityImpairments = zeros(length(SNR),nbrOfRealizations,length(Nt)); %SISO channels


%Go through all channel realizations
for k = 1:nbrOfRealizations
    
    %Go through all numbers of transmit antennas
    for l = 1:length(Nt)
        Ntcurrent = Nt(l); %Current number of transmit antennas
        
        %Squared singular values of current channel realization 
        %(this is the same as the eigenvalues of the Gram matrix)
        eigenvaluesMIMO = svd(H(:,1:Ntcurrent,k)).^2;
        
        %Inversion of the squared singular values, which is used in the
        %waterfilling algorithm for power allocation
        waterfillingFactorsIdealHardware = 1./(eigenvaluesMIMO+1e-10);
        
        %Compute capacity with ideal hardware in the SISO case. It is the
        %same for both deterministic and random channels, since there is no
        %power allocation. Note that all SNR values are handled.
        SISOcapacityIdealHardware = log2(1 + SNR'*abs(h').^2 );
        
        %Compute capacity with transceiver hardware impairments in the SISO
        %case. It is the same for both deterministic and random channels,
        %since there is no power allocation. Note that all SNR values are handled.
        SISOcapacityImpairments = log2(1 + (SNR'*abs(h').^2)./(kappa^2*SNR'*abs(h').^2+1) );
        
        
        %Go through all SNR values
        for m = 1:length(SNR)
            
            %Compute capacity-achieving power allocation for deterministic
            %channels, using the waterfilling algorithm from [1]. Corollary
            %2 shows that it is also the capacity-achieving under
            %transceiver hardware impairments with alpha=1.
            powerAllocationDeterministic = functionWaterfilling(SNR(m),waterfillingFactorsIdealHardware);
 
            %Compute capacity with ideal hardware and a deterministic channel
            %realization, using the power allocation from above.
            MIMOcapacityIdealHardware_deterministic(m,k,l) = sum(log2(1+powerAllocationDeterministic./waterfillingFactorsIdealHardware));
            
            %Compute capacity with transceiver hardware impairments for alpha=1
            %and a deterministic channel realization. It is based on 
            %Eq. (6) and uses the power allocation from above.
            MIMOcapacityImpairments_deterministic(m,k,l) = sum(log2(1+powerAllocationDeterministic.*eigenvaluesMIMO ./ ( kappa^2/Ntcurrent*SNR(m)*eigenvaluesMIMO+1)));
            
            %Compute capacity with ideal hardware and one realization of
            %uncorrelated Rayleigh fading. The averaging is performed when
            %plotting the results.
            MIMOcapacityIdealHardware_random(m,k,l) = real(log2(det(eye(Nr)+(SNR(m)/Ntcurrent)*H(:,1:Ntcurrent,k)*H(:,1:Ntcurrent,k)')));
            
            %Compute capacity with transceiver hardware impairments and one
            %realization of uncorrelated Rayleigh fading. Equal power
            %allocation is optimal, according to Corollary 1. The averaging
            %is performed when plotting the results.
            MIMOcapacityImpairments_random(m,k,l) = real(log2(det(eye(Nr)+(SNR(m)/Ntcurrent+SNR(m)*kappa^2/Ntcurrent)*H(:,1:Ntcurrent,k)*H(:,1:Ntcurrent,k)'))-log2(det(eye(Nr)+(SNR(m)*kappa^2/Ntcurrent)*H(:,1:Ntcurrent,k)*H(:,1:Ntcurrent,k)')));
            
        end
    end
end


%Plot Figure 4 from the paper
figure(1); hold on; box on;

for m = 1:length(Nt)
    plot(SNRdB,mean(MIMOcapacityIdealHardware_random(:,:,m),2)./mean(SISOcapacityIdealHardware,2),'k'); %Compute and plot finite-SNR multiplexing gain with ideal hardware
    plot(SNRdB,mean(MIMOcapacityImpairments_random(:,:,m),2)./mean(SISOcapacityImpairments,2),'r--'); %Compute and plot finite-SNR multiplexing gain with transceiver hardware impairments
end

plot(SNRdB,Nr*ones(size(SNRdB)),'k:'); %Plot high-SNR limit of finite-SNR multiplexing gain for transceiver hardware impairments and ideal hardware (see Eq. (10))

text(SNRdB(21),mean(MIMOcapacityIdealHardware_random(21,:,1),2)./mean(SISOcapacityIdealHardware(21,:),2)+0.02,['N_t = ' num2str(Nt(1))]);
text(SNRdB(21),mean(MIMOcapacityIdealHardware_random(21,:,2),2)./mean(SISOcapacityIdealHardware(21,:),2)+0.08,['N_t = ' num2str(Nt(2))]);
text(SNRdB(21),mean(MIMOcapacityIdealHardware_random(21,:,3),2)./mean(SISOcapacityIdealHardware(21,:),2)+0.15,['N_t = ' num2str(Nt(3))]);

xlabel('SNR [dB]')
ylabel('Finite-SNR Multiplexing Gain')
axis([SNRdB(1) SNRdB(end) M-0.5 M+0.5]);
legend('Ideal Transceiver Hardware','Transceiver Impairments: \kappa=0.05','Location','SouthEast')



%Plot Figure 5 from the paper
figure(2); hold on; box on;

for m = 1:length(Nt)
    plot(SNRdB,mean(MIMOcapacityIdealHardware_deterministic(:,:,m),2)./mean(SISOcapacityIdealHardware,2),'k'); %Compute and plot finite-SNR multiplexing gain with ideal hardware
    plot(SNRdB,mean(MIMOcapacityImpairments_deterministic(:,:,m),2)./mean(SISOcapacityImpairments,2),'b-.'); %Compute and plot finite-SNR multiplexing gain with transceiver hardware impairments
    
    highSNRlimit = M*log2(1+Nt(m)/M/kappa^2)/log2(1+1/kappa^2); %High SNR limit from Eq. (10) for hardware impairments
    plot(SNRdB,highSNRlimit*ones(size(SNRdB)),'k:'); %Plot high-SNR limit of finite-SNR multiplexing gain
end

plot(SNRdB,M*ones(size(SNRdB)),'k:'); %Plot high-SNR limit of finite-SNR multiplexing gain for ideal hardware

text(SNRdB(41),mean(MIMOcapacityIdealHardware_deterministic(41,:,1),2)./mean(SISOcapacityIdealHardware(41,:),2)+0.08,['N_t = ' num2str(Nt(1))]);
text(SNRdB(41),mean(MIMOcapacityIdealHardware_deterministic(41,:,2),2)./mean(SISOcapacityIdealHardware(41,:),2)+0.05,['N_t = ' num2str(Nt(2))]);
text(SNRdB(41),mean(MIMOcapacityIdealHardware_deterministic(41,:,3),2)./mean(SISOcapacityIdealHardware(41,:),2)+0.08,['N_t = ' num2str(Nt(3))]);

xlabel('SNR [dB]')
ylabel('Average Finite-SNR Multiplexing Gain')
axis([SNRdB(1) SNRdB(end) 3.5 M+4]);
legend('Ideal Transceiver Hardware','Transceiver Impairments: \kappa=0.05','Location','NorthEast')
