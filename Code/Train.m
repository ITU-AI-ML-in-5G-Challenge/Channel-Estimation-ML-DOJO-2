%Calculate NMSE of the channel
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
fprintf("Best threshold for data_set=%i : %1.4f with %3.2fdB NMSE\n", data_set, TH_domain(imax), amax)

print('-f10', sprintf('Training_%i', data_set), '-dpng')

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