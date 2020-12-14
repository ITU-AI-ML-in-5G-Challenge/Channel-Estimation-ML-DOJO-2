function [CH, TH] = ReconstructChannel_train(Pilots, Phi, N_RX, N_TX, D_w)
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
Time_sr_d = 4;             % Time sampling ratio
Time_sr_r = 32;            % Time sampling ratio for refinement
Angle_sr_d = 2;            % Angle sampling ratio for detection
Angle_sr_r = 32;           % Angle sampling ratio for refinement
Confidence_min = 0.7;      % Detection confidence minimum
PingPong_it = 3;           % Number of iteration to refine the angles
bool_debug = false;        % Debug condition for plots
%% Whitening
if nargin > 4
    Pilots = D_w\Pilots;
    Phi = D_w\Phi;
end
%% Compute N_Phi and N_fft
[N_Phi, N_fft] = size(Pilots);
%% Derived parameters
Time_res_d = Time_sr_d*N_fft;
Time_res_r = Time_sr_r*N_fft;
%% Compute angle responses
Angle_RX_r = linspace(-pi, pi, N_RX*Angle_sr_r+1); Angle_RX_r(end) = [];
Angle_TX_r = linspace(-pi, pi, N_TX*Angle_sr_r+1); Angle_TX_r(end) = [];
A_RX_r = exp((0:N_RX-1).'*Angle_RX_r*1i)/sqrt(N_RX);
A_TX_r = exp((0:N_TX-1).'*Angle_TX_r*1i)/sqrt(N_TX);
%% Compute simple angle transform
if Angle_sr_d < 2
    Angle_d = linspace(-pi, pi, N_RX*Angle_sr_d+1); Angle_d(end) = [];
    A_RX_d = SincBeam(N_RX, 1.1*2*pi/(N_RX*Angle_sr_d)).*exp((0:N_RX-1).'*Angle_d*1i);
    Angle_d = linspace(-pi, pi, N_TX*Angle_sr_d+1); Angle_d(end) = [];
    A_TX_d = SincBeam(N_TX, 1.1*2*pi/(N_TX*Angle_sr_d)).*exp((0:N_TX-1).'*Angle_d*1i);
else
    Angle_d = linspace(-pi, pi, N_RX*Angle_sr_d+1); Angle_d(end) = [];
    A_RX_d = exp((0:N_RX-1).'*Angle_d*1i);
    Angle_d = linspace(-pi, pi, N_TX*Angle_sr_d+1); Angle_d(end) = [];
    A_TX_d = exp((0:N_TX-1).'*Angle_d*1i);
end
A_RTX_d = kron(A_TX_d.', A_RX_d');
%% Dump pilots info into the residual pilots and convert it to time information
Pilots_res = Pilots;
Pilots_res_t = ifft(Pilots_res, Time_res_d, 2);
Pilots_res_t_esp = (A_RTX_d*Phi')*Pilots_res_t;
%% Compute noise level and the confidence threshold
Noise_level = sqrt(1/(2*log(2)))*median(abs(Pilots_res_t_esp(:)), 1);
Confidence_TH = Noise_level*sqrt(-2*log(1-Confidence_min^(1/numel(Pilots_res_t_esp))));
%% Plot noise level and std margins
if bool_debug
    figure(1)
    plot(max(abs(Pilots_res_t_esp)).', 'b', 'LineWidth', 1.5), hold on
    plot([1, Time_res_d], repmat(Confidence_TH, 1, 2), 'r', 'LineWidth', 1);hold off
end
%% Initialization
CH = cell(0);
TH = [];
%% Path substraction loop
AoA = [];
AoD = [];
ToF = [];
H   = [];
P   = [];
Alpha = [];
[amax, imax] = max(abs(Pilots_res_t_esp(:)));
while amax > Confidence_TH
    %% Imaging
    [ii_aoa, ii_aod, ii_t] = ind2sub([N_RX*Angle_sr_d, N_TX*Angle_sr_d, Time_res_d], imax);
    %% Extract time and spatial measurements
    t = (ii_t-1)/Time_res_d; % Time computed between 0 and 1
    Mt = Pilots_res_t(:, ii_t);
    %% Pseudochannel computation
    % These two options are similar but have different properties
    H_pseudo = reshape(Phi'*Mt, [N_RX, N_TX]);
%     H_pseudo = reshape(Phi\Mt, [N_RX, N_TX]);
    %% First iteration estimation
    if N_RX > N_TX
        [~, ii_aoa] = max(abs((H_pseudo*A_TX_d(:, ii_aod))'*A_RX_r));
        [~, ii_aod] = max(abs((A_RX_r(:, ii_aoa)'*H_pseudo)*A_TX_r));
    else
        [~, ii_aod] = max(abs((A_RX_d(:, ii_aoa)'*H_pseudo)*A_TX_r));
        [~, ii_aoa] = max(abs((H_pseudo*A_TX_r(:, ii_aod))'*A_RX_r));
    end
    [~, ii_t] = max(abs(ifft((kron(A_TX_r(:, ii_aod)', A_RX_r(:, ii_aoa).')*Phi')*Pilots_res, Time_res_r, 2)));
    %% Ping pong
    for iter = 1:PingPong_it
        if N_RX > N_TX
            [~, ii_aoa] = max(abs((H_pseudo*A_TX_r(:, ii_aod))'*A_RX_r));
            [~, ii_aod] = max(abs(A_RX_r(:, ii_aoa)'*H_pseudo*A_TX_r));
        else
            [~, ii_aod] = max(abs(A_RX_r(:, ii_aoa)'*H_pseudo*A_TX_r));
            [~, ii_aoa] = max(abs((H_pseudo*A_TX_r(:, ii_aod))'*A_RX_r));
        end
        [~, ii_t] = max(abs(ifft((kron(A_TX_r(:, ii_aod).', A_RX_r(:, ii_aoa)')*Phi')*Pilots_res, Time_res_r, 2)));
    end
    %% Angular values and normallized time in [0, 1[
    aoa = Angle_RX_r(ii_aoa);
    aod = Angle_TX_r(ii_aod);
    t = (ii_t-1)/Time_res_r; % Time computed between 0 and 1
    %% Reconstruct pilots measurement
    h = reshape(A_RX_r(:, ii_aoa)*A_TX_r(:, ii_aod)', [], 1)*exp(-(0:N_fft-1)*(t*2i*pi));
    p = (Phi*reshape(A_RX_r(:, ii_aoa)*A_TX_r(:, ii_aod)', [], 1))*exp(-(0:N_fft-1)*(t*2i*pi));
    %% Dump data
    AoA = [AoA, aoa];
    AoD = [AoD, aod];
    ToF = [ToF, t];
    H   = [H, h(:)];
    P   = [P, p(:)];
    Alpha = P \ Pilots(:);
    %% Reconstruct channel
    Ch = reshape(H*Alpha, [N_RX, N_TX, N_fft]);
    CH{end+1} = Ch;
    TH(end+1) = (1-exp(-(amax/Noise_level)^2/2))^numel(Pilots_res_t_esp);
    %% Update pilot residual
    Pilots_res = reshape(Pilots(:) - P*Alpha, size(Pilots));
    Pilots_res_t = ifft(Pilots_res, Time_res_d, 2);
    Pilots_res_t_esp = (A_RTX_d*Phi')*Pilots_res_t;
    %% Compute noise level and the confidence threshold
    Noise_level = sqrt(1/(2*log(2)))*median(abs(Pilots_res_t_esp(:)), 1);
    Confidence_TH = Noise_level*sqrt(-2*log(1-Confidence_min^(1/numel(Pilots_res_t_esp))));
    %% Plot reconstruction and residual
    if bool_debug
        figure(2)
        plot(max(abs(ifft(reshape(P*Alpha, size(Pilots)), Time_res_d, 2)), [], 1), 'b', 'LineWidth', 1.5)
        figure(3)
        plot(max(abs(Pilots_res_t_esp), [], 1), 'b', 'LineWidth', 1.5), hold on
        plot([1, Time_res_d], repmat(Confidence_TH, 1, 2), 'r', 'LineWidth', 1); hold off
    end
    %% Recompute max
    [amax, imax] = max(abs(Pilots_res_t_esp(:)));
end

end
