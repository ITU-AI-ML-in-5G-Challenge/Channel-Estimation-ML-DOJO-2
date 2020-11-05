
%System parameters
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

% Training parameters
%Nc = 30; % Number of channels available in the channel data set up to 10000


Ntrain = 20; % Number of training symbols to be received for each one of the available channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate training precoders and combiners and matrix Phi in (13) in [1]
% They are generated as pseudorandom phase shifts. 
% We use a different training precoder and combiner for every training
% symbol to be transmitted through the same channel.
% We use the same set of training precoders and combiners for every channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nres = 2^Nbits; %resolution ofr the phase shifters
Phi=zeros(Ntrain*Lr,Nt*Nr);%Initialize measurement matrix Phi in [1,(10)] of size LrxNtNr.
                           % We have a different measurement matrix for every training symbol.
                           % Here we are initilizaing the measurement
                           % matrices for all the training symbols, which
                           % correspond to matrix Phi in (13)
Wtr = Wtr_save;
Ftr = Ftr_save;
for i=1:Ntrain
   signal = signal_save(:,i); %training signal q (frequency flat)
   Phi((i-1)*Lr+(1:Lr),:)=kron(signal.'*Ftr(:,(i-1)*Lt+(1:Lt)).',Wtr(:,(i-1)*Lr+(1:Lr))');% Generate Phi in (13)
end

Pilots_imag = h5read('test_dataset_v3_20_pilots_1_data_set.hdf5','/training_data_imag');
Pilots_real = h5read('test_dataset_v3_20_pilots_1_data_set.hdf5','/training_data_real');
Pilots = Pilots_real + 1i*Pilots_imag;
N_RX = Nr;
N_TX = Nt;
P_s = size(Pilots);
Nc = P_s(1); 
Ch_all = zeros(Nc,64,16,256);
for i =1:Nc
    Pilotsi = squeeze(Pilots(i,:,:));
    [Ch, AoA, AoD, ToF, Alpha, NUM] = ReconstructChannel_v2(Pilotsi, Phi, N_RX, N_TX);
    Ch_all(i,:,:,:)=Ch;
end   

training_data_real = real(Ch_all);
training_data_imag = imag(Ch_all);

chan_save_file_hdf5 = 'Channels_test_dataset_v3_20_pilots_1_data_set.hdf5'

if ~isfile(chan_save_file_hdf5)
    h5create(chan_save_file_hdf5, '/training_data_real', size(Ch_all));
    h5create(chan_save_file_hdf5, '/training_data_imag', size(Ch_all));
end

h5write(chan_save_file_hdf5, '/training_data_real', training_data_real);
h5write(chan_save_file_hdf5, '/training_data_imag', training_data_imag);

%%
% % 
% % Pilots_imag1 = h5read('Channels_test_dataset_v3_20_pilots_1_data_set.hdf5','/training_data_imag');
% % Pilots_real1 = h5read('Channels_test_dataset_v3_20_pilots_1_data_set.hdf5','/training_data_real');
% % 
% %     
   