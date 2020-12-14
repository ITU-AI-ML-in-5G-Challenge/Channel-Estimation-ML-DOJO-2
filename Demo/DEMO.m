%% Clears
clc
clear

%% Training parameters
Nc = 50; % Number of channels available in the channel data set up to 10000
Ntrain = 100; % Number of training symbols to be received for each one of the available channels
data_set = 1; %1, 2 or 3

%% System parameters
Nt = 16; % Number of TX antennas
Nr = 64; % Number of RX antennas
Nbits= 4; % Number of bits available to represent a phase shift in the analog precoder/combiner.
Lt = 2;  % Number of TX RF chains 
Lr = 4;  % Number of RX RF chains 
Ns = 2;  % Number of data streams to be transmitted
Nfft=256; % Number of subcarriers in the MIMO-OFDM system
Pt=1; % Transmit power(mw)
Nfilter = 20;
Mfilter = 1; %no oversampling
rolloff = 0.8;
MHz = 1e6; 
fs = 1760*MHz; %Sampling frequency
Ts = 1/fs;

%% Display
fprintf("Generating %i samples for data-set %i.", Nc, data_set)
tic

%% SNR
switch data_set
    case 1
        SNR = -15;
    case 2
        SNR = -10;
    case 3
        SNR = -5;
    otherwise 
        error('data_set not found')
end %switch
snr = 10.^(SNR/10);

%% Code
Nres = 2^Nbits; %resolution ofr the phase shifters
Phi=zeros(Ntrain*Lr,Nt*Nr);%Initialize measurement matrix Phi in [1,(10)] of size LrxNtNr.
                           % We have a different measurement matrix for every training symbol.
                           % Here we are initilizaing the measurement
                           % matrices for all the training symbols, which
                           % correspond to matrix Phi in (13)
rng(1);
tt=randi(Nres,[Nt Ntrain*Lt]);
for i = 1:Nres
   tt(tt==i) = exp(1i*2*pi*(i-1)/Nres);
end
Ftr=tt;  
        
tt=randi(Nres,[Nr Ntrain*Lr]);
for i = 1:Nres
   tt(tt==i) = exp(1i*2*pi*(i-1)/Nres);
end
Wtr=tt;    

Ftr = Ftr/sqrt(Nt);% pseudorandom training precoders (generating Ntrain, precoders for all the training symbols)
Wtr = Wtr/sqrt(Nr);% pseudoranmdom training combiners

save TrainingPrecoders.mat Ftr
save TrainingCombiners.mat Wtr
        
for i=1:Ntrain
   signal = sqrt(1/2/Lt)*(sign(randn(Lt,1))+1i*sign(randn(Lt,1))); %training signal q (frequency flat)
   Phi((i-1)*Lr+(1:Lr),:)=kron(signal.'*Ftr(:,(i-1)*Lt+(1:Lt)).',Wtr(:,(i-1)*Lr+(1:Lr))');% Generate Phi in (13)
end



r = zeros(Ntrain*Lr,Nfft);%Initialize RX training symbols for one channel
R = zeros(Nc,Ntrain*Lr,Nfft);% Initialize RX training synbols for all the channels
Channels = zeros(Nc,Nr,Nt,Nfft);
nn = zeros(Lr*Ntrain,Nfft);% Initialize noise matrix at the RF combiner output for all the training symbols

for j=1:Nc %Nc is number of channels 
    [Hk,H_time,At,Ar] = gen_channel_ray_tracing(j,Nr,Nt,Nfft,Ts,rolloff,Mfilter); 
    Channels(j,:,:,:) =  Hk;
    % Load channel parameters for channel j and build channel matrix including filtering
    % effects
   
   
    var_n = Pt/snr;
    Noise = sqrt(var_n/2)*(randn(Nr,Ntrain,Nfft)+1i*randn(Nr,Ntrain,Nfft));
    SNRaux = zeros(Nfft,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRLB calc
%     upsilon = Phi*kron(conj(At),Ar); %% new
%     C = noise_cov_matrix(Wtr,Ntrain,Lr,var_n); % new
%     crlb = CRLB(C,At,Ar,upsilon) % new
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:Nfft % Generate RX pilots for every subcarrier
            for t=1:Ntrain
                Wrf_t = Wtr(:,(t-1)*Lr+(1:Lr));
                nn((1:Lr)+Lr*(t-1),k) = Wrf_t'*Noise(:,t,k);
            end
            signal_k = Phi*reshape(Hk(:,:,k),[],1);
            noise_k = nn(:,k);
            r(:,k) = Phi*reshape(Hk(:,:,k),[],1) + nn(:,k);
            SNRaux(k) = signal_k'*signal_k/(noise_k'*noise_k);
    end
    Average_SNR = 10*log10(mean(SNRaux));
    R(j,:,:)=r;
end
fprintf("   |    Finished (%.1f s)\n", toc);

%% Training
fprintf("Training.")
tic

%% Calculate NMSE of the channel
TH_domain = linspace(0.7, 1, 2^12);
nmse = zeros(size(TH_domain));
iterations = zeros(Nc,1);
D_w = Whitening(Wtr,Ntrain,Lr);
size_R = size(R);
for ex=1:Nc
    Pilots = reshape(R(ex,:,:),size_R(2),size_R(3));
    H = reshape(Channels(ex,:,:,:), Nr,Nt,Nfft);
    [CH, TH] = ReconstructChannel_train(Pilots, Phi, Nr, Nt, D_w);
    TH = [TH, 0.7];
    for ii = 1:length(CH)
        th1 = TH(ii);
        th2 = TH(ii+1);
        Ch = CH{ii};
        nmse_Ch=NMSE_channel(Ch,H, Nfft);
        nmse(TH_domain >= th2 & TH_domain < th1) = nmse(TH_domain >= th2 & TH_domain < th1) + 10*log10(nmse_Ch);
    end
end
nmse = nmse/Nc;

figure(10);
plot(TH_domain, nmse,'r')
switch data_set
    case 1
        ylim([-19, -17.5])
    case 2
        ylim([-23, -20])
    case 3
        ylim([-25.5, -24])
end
xlabel("Detection threshold \theta")
ylabel("NMSE [dB]")

[amax, imax] = min(nmse);
fprintf("   |    Finished (%.1f s)\n", toc);
fprintf("Best threshold for data_set=%i : %1.4f with %3.2fdB NMSE\n", data_set, TH_domain(imax), amax)

print('-f10', sprintf('Training_%i_new', data_set), '-dpng')

%% calculate NMSE
function nmse = NMSE_channel(H_hat,H, Nffc)
% INPUTS    
% H_hat - Channel estimate
% H - real channel
%OUTPUTS NMSE
    nmse = 0;
    den = 0;
    for subcarrier=1:Nffc
        sub = H_hat(:,:,subcarrier)-H(:,:,subcarrier);
        H_k = H(:,:,subcarrier);
        nmse = nmse + norm(sub(:),'fro')^2;
        den = den + norm(H_k(:),'fro')^2;
    end
    nmse = nmse/den;
end