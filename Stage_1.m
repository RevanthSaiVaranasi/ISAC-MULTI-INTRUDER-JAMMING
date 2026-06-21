%% ==========================================================
% Stage-1 : Multi-Intruder Detection and Angle Estimation
%
% Purpose:
%   1. Detect multiple intruders using Capon Beamforming.
%   2. Estimate Direction of Arrival (DoA).
%   3. Compute CRB-based uncertainty intervals.
%   4. Generate sampled angle sets for Algorithm-2.
%
% Inputs:
%   Nt, Nr, Ps, sigma2, Lmax, true_theta, true_gain, dist
%
% Outputs:
%   theta_samples
%   Angle_Info
%   Algo1_Output.mat
%
% Author: Revanth Sai Varanasi
% ==========================================================
clc;
clear;
close all;

%% User Defined Values
N = input('Enter the number of Intruders: ');
Nt = input('Enter the Number of Transmit Antennas: ');
Nr = input('Enter the Number of Recieve Antennas: ');
Lmax = input('Enter the Maximum number of Sensing Beams: ');
Pmax = input('Enter the Maximum Power available at Base Station: ');
true_theta = zeros(N,1);
for n = 1 : N
    true_theta(n) = input(sprintf('Enter the true angle of intruder %d:',n));
end
true_gain = zeros(N,1);
for n = 1: N
    true_gain(n) = input(sprintf('Enter the true gain of intruder %d:',n));
end
Ps = input('Enter the transmit Power for Sensing at Base Station: ');
Sigma_2 = input('Enter the Noise Power of AWGN: ');
Sth = input('Enter the Maximum Threshold Power of the intruders: ');
delta = input('Enter the value of Minimum Criteria for Binary Search Method: ');
Z = input('Enter the Number of Noise Samples: ');
lambda = input('Enter the Wavelength of the Signal: ');
dist = zeros(N,1);
for n = 1: N
    dist(n) = input(sprintf('Enter the distance of the intruder %d:',n));
end
d = lambda/2;
Qthetha = 3;

%% Binary Search for Optimal Number of Sensing Beams
feasible_solution_found = false;
Llow = Nt;
Lup = Lmax;
P_final = [];
thetha_test = -90:0.1:90;

%% Generate DFT Sensing Beam Directions
while(abs(Llow - Lup) > delta)
    L = floor((Llow + Lup) / 2);
    fprintf('Testing for L = %d\n', L);
     
    K_last = 0;
    theta_hat_last = [];
    P_last = [];
    for z = 1:Z
        thetha_scan = zeros(1,L);
        %% Construct Received Signal Matrix
        for l = 1:L
            thetha_scan(l) = asind(-1 + (2*l-1)/L); 
        end
        %% Covariance Matrix Computation
        Y = zeros(Nr, L);
        X = zeros(Nt, L);
        
        for l = 1:L
            ns = sqrt(sigma2/2)*(randn(Nr,1) + 1j*randn(Nr,1)); % Recieved AWG noise every time 
            v = (1/sqrt(Nt))*steering_vector(Nt, thetha_scan(l)); 
            xs = 1;
            Y(:,l) = ns;
            for n = 1:N
                at = steering_vector(Nt, true_theta(n));
                ar = steering_vector(Nr, true_theta(n));
                Y(:,l) = Y(:,l) + sqrt(Ps)*ar*true_gain(n)*(at.')*v*xs;
            end
            X(:,l) = sqrt(Ps/Nt)*steering_vector(Nt, thetha_scan(l));
        end

        %Co-Varaince Matrices
        Rxx = (1/L)*Y*(Y');
        Rx = (1/L)*X*(X');

        % Capon spectrum
        P = zeros(size(thetha_test));
        for p = 1:length(thetha_test)
            ar_test = steering_vector(Nr, thetha_test(p));
            Rxx_reg = Rxx + 1e-6*eye(size(Rxx));
            P(p) = 1/real(ar_test'*(Rxx_reg\ar_test));
        end
        
        Pnorm = P/max(P);
        [pks,locs] = findpeaks(Pnorm,'MinPeakProminence',1e-3,'MinPeakDistance',20);
        theta_hat = thetha_test(locs);
        K = length(theta_hat);
        
        %detection check
        if K == 0
            fprintf('z=%d: No targets detected, skipping.\n', z);
            break;
        end
        if K < N
            fprintf('Only %d out of %d intruders detected for L = %d\n',K,N,L);
            break;
        end
        % Steering matrices
        At = zeros(Nt, K);
        Br = zeros(Nr, K);
        for k = 1:K
            At(:,k) = steering_vector(Nt, theta_hat(k));
            Br(:,k) = steering_vector(Nr, theta_hat(k));
        end

        % Gain estimation with regularization
        inner     = At.' * Rx * conj(At);
        reg_inner = inner + 1e-8*eye(size(inner));
        gamma     = Y*Y' - (1/K)*Y*X'*conj(At)*(reg_inner\(At.'*X*Y'));
        gamma_reg = gamma + 1e-6*norm(gamma,'fro')*eye(size(gamma));

        num        = (Br'*(gamma_reg\Br)) .* (At.'*conj(Rx)*At);
        den_matrix = Br'*(gamma_reg\(Y*X'*conj(At)));
        den        = diag(den_matrix);
        Lambda_hat = (1/L)*(num\den);

        %CRB and uncertainty intervals of every intruder
        CRB = zeros(K,1);
        xi1 = zeros(K,1);
        xi2 = zeros(K,1);
        for k = 1:K
            CRB(k) = sigma2/(2*L*Ps*(abs(Lambda_hat(k))^2));
            xi1(k) = theta_hat(k) - 3*sqrt(CRB(k));
            xi2(k) = theta_hat(k) + 3*sqrt(CRB(k));
        end

        %Last valid result for plotting
        K_last = K;
        theta_hat_last = theta_hat;
        P_last = P;

    end 

    % Bisection decision
    if K_last < N
        Llow = L;
    else    
        feasible_solution_found = true;
        Lup = L;
        P_final = P_last;
    end

end 

if feasible_solution_found
    fprintf('Optimal Number of Sensing Beams = %d\n',Lup);
    Angle_Info = zeros(K,3);
    for k = 1:K
        Angle_Info(k,1) = theta_hat(k);
        Angle_Info(k,2) = xi1(k);
        Angle_Info(k,3) = xi2(k);
    end
    disp('Angle Information Matrix')
    disp('   Theta_hat      Xi1         Xi2')
    disp(Angle_Info)
    theta_samples = zeros(K,Qtheta);
    for k = 1:K
        theta_samples(k,:) = linspace(xi1(k),xi2(k),Qtheta);
    end
    %% Save Outputs for Stage-2
    save('Algo1_Output.mat','Nt','Sth','dist','lambda','sigma2','Pmax',...
        'Angle_Info','theta_samples');
end

% Capon Spectrum plot
if ~isempty(P_final)
    figure;
    plot(thetha_test,P_final,'LineWidth',2);
    xlabel('Angle (degrees)');
    ylabel('Capon Spectrum');
    title(sprintf('Capon Spectrum at Optimal L = %d',Lup));
    grid on;

    % Power graph
    figure;
    plot(thetha_test,10*log10(P/max(P)),'LineWidth',2);
    grid on;
    xlabel('Angle (deg)');
    ylabel('Normalized Power (dB)');
else
    disp(theta_hat)
    fprintf('No valid Capon spectrum available for plotting.\n');
end

