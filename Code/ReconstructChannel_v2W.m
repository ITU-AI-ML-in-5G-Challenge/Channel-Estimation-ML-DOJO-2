function [Ch, AoA, AoD, ToF, Alpha] = ReconstructChannel_v2W(Pilots, Phi, N_RX, N_TX,D_w)
% function R = ReconstructChannel(Pilots, Phi, N_RX, N_TX)
%
% INPUTS
% Pilots - Measured antenna configurations with dimensions (N_Phi, N_fft)
% Phi - Used antenna configurations with dimensions (N_Phi, N_RX*N_TX)
% N_RX - Number of antenna elements in the receiver
% N_TX - Number of antenna elements in the transmitter
% Note: N_Phi is the number of used antenna configurations and N_fft the
% number of subcarriers
%
% OUTPUTS
% Ch - Reconstructed channel
% AoA - Angles of arrival for the computed paths
% AoD - Angles of departure for the computed paths
% ToF - Times of flight for the computed paths
% Alpha - Complex gains for the computed paths

%% Parameters
Time_res = 2^12;         % Time resolution
Angle_res = 2^12;        % Angle resolution
Angle_oversampling = 2;  % Angle oversampling for detection
Confidence = 0.995;      % Detection confidence
PingPong_it = 3;         % Number of iteration to refine the angles
bool_debug = false;       % Debug condition for plots
%% Compute N_Phi and N_fft
Pilots = inv(D_w)*Pilots;
Phi = inv(D_w)*Phi;
[N_Phi, N_fft] = size(Pilots);
%% Compute angle responses
Angle = linspace(-pi, pi, Angle_res+1); Angle(end) = [];
A_RX = exp((0:N_RX-1).'*Angle*1i);
A_TX = exp((0:N_TX-1).'*Angle*1i);
%% Compute simple angle transform
Angle_simple = linspace(-pi, pi, N_RX*Angle_oversampling+1); Angle_simple(end) = [];
A_RX_simple = exp((0:N_RX-1).'*Angle_simple*1i);
Angle_simple = linspace(-pi, pi, N_TX*Angle_oversampling+1); Angle_simple(end) = [];
A_TX_simple = exp((0:N_TX-1).'*Angle_simple*1i);
A_RTX = kron(A_TX_simple.', A_RX_simple.');
%% Dump pilots info into the residual pilots and convert it to time information
Pilots_res = Pilots;
Pilots_res_t = ifft(Pilots_res, Time_res, 2);
Pilots_res_t_norm = sum(abs(Pilots_res_t).^2, 1);
Pilots_res_t_esp = (A_RTX*Phi')*Pilots_res_t;
%% Compute noise level and its standard deviation using Rayleigh distribution properties
Noise_level = sqrt(pi/(4*log(2)))*median(abs(Pilots_res_t_esp(:)), 1);
Noise_std = sqrt((2-pi/2)/(pi/2))*Noise_level;
%% Plot noise level and std margins
if bool_debug
    figure(1)
    plot(max(abs(Pilots_res_t_esp)).', 'b', 'LineWidth', 1.5), hold on
    plot([1, Time_res], repmat(Noise_level, 1, 2), 'k', 'LineWidth', 1);
    plot([1, Time_res], repmat(Noise_level+Noise_std, 1, 2), 'r', 'LineWidth', 0.5);
    plot([1, Time_res], repmat(Noise_level-Noise_std, 1, 2), 'r', 'LineWidth', 0.5); hold off
end
%% Path substraction loop
[amax, imax] = max(max(abs(Pilots_res_t_esp), [], 1), [], 2);
LoS_P = amax;
AoA = [];
AoD = [];
ToF = [];
H   = [];
P   = [];
Confidence_TH = Noise_level*sqrt(-2*log(1-Confidence^(1/size(Pilots_res_t_esp, 1))));
while amax > Confidence_TH
    %% Extract time and spatial measurements
    t = (imax-1)/Time_res; % Time computed between 0 and 1
    Mt = Pilots_res_t(:, imax);
    %% Pseudochannel computation
    % These two options are similar but have different properties
    H_pseudo = reshape(Phi'*Mt, [N_RX, N_TX]);
%     H_pseudo = reshape(Phi\Mt, [N_RX, N_TX]);
    %% aoa or aod imaging
    H_Angle = abs(A_RX'*H_pseudo*A_TX);
    [amax, imax] = max(H_Angle(:));
    [ii_aoa, ii_aod] = ind2sub([Angle_res, Angle_res], imax);
    if bool_debug
        figure(4)
        imshow(H_Angle/amax)
    end
    %% Angular values
    aoa = Angle(ii_aoa);
    aod = Angle(ii_aod);
    %% Reconstruct pilots measurement
    h = reshape(A_RX(:, ii_aoa)*A_TX(:, ii_aod)', [], 1)*exp(-(0:N_fft-1)*(t*2i*pi));
    p = (Phi*reshape(A_RX(:, ii_aoa)*A_TX(:, ii_aod)', [], 1))*exp(-(0:N_fft-1)*(t*2i*pi));
    %% Dump data
    AoA = [AoA, aoa];
    AoD = [AoD, aod];
    ToF = [ToF, t];
    H   = [H, h(:)];
    P   = [P, p(:)];
    Alpha = P \ Pilots(:);
    %% Update pilot residual
    Pilots_res = reshape(Pilots(:) - P*Alpha, size(Pilots));
    Pilots_res_t = ifft(Pilots_res, Time_res, 2);
    Pilots_res_t_norm = sum(abs(Pilots_res_t).^2, 1);
    Pilots_res_t_esp = (A_RTX*Phi')*Pilots_res_t;
    %% Plot reconstruction and residual
    if bool_debug
        figure(2)
        plot(max(abs(ifft(reshape(P*Alpha, size(Pilots)), Time_res, 2)), [], 1), 'b', 'LineWidth', 1.5)
        figure(3)
        plot(max(abs(Pilots_res_t_esp), [], 1), 'b', 'LineWidth', 1.5), hold on
        plot([1, Time_res], repmat(Confidence_TH, 1, 2), 'r', 'LineWidth', 1); hold off
    end
%     %% While building function
%     break
    %% Recompute max
    [amax, imax] = max(max(abs(Pilots_res_t_esp), [], 1), [], 2);
end
%% Reconstruct channel
Ch = reshape(H*Alpha, [N_RX, N_TX, N_fft]);

end
