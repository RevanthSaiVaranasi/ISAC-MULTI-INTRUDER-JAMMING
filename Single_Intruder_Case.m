%% ==========================================================
% Algorithm 1 : Robust Beamforming for Single Intruder
%
% Description:
%   This MATLAB implementation realizes Algorithm 1 proposed in
%   the Robust Sensing-Aware Jamming (SAJA) framework.
%
%   The algorithm performs:
%     1. Binary search to obtain the minimum feasible sensing beams.
%     2. DFT beam scanning for sensing.
%     3. Capon beamforming for DoA estimation.
%     4. CAML channel gain estimation.
%     5. CRB computation for angle uncertainty.
%     6. Robust beamformer design.
%     7. Worst-case JNR evaluation.
%
% Output:
%   - Optimal number of sensing beams
%   - Capon Spectrum
%
% Author : Revanth Sai Varanasi
%% ==========================================================

clc;
clear;
close all;

%% User Inputs
% System parameters, sensing parameters and target information.
Nt = input('Enter number of transmit antennas: ');
Nr = input('Enter number of receive antennas: ');
Pmax = input('Enter maximum Base Station power: ');
Ps = input('Enter sensing power: ');
sigma2 = input('Enter noise power: ');
Sth = input('Enter threshold Power: ');
Lmax = input('Enter maximum sensing beams: ');
delta = input('Enter binary search stopping criteria: ');
Z = input('Enter number of Monte Carlo noise samples: ');
true_theta = input('Enter true angle of Intudder : ');
true_gain = input('Enter true gain of Intruder : ');
lambda = input('Enter wavelength: ');
d = input('Enter antenna spacing: ');
dist = input('Enter target distance: ');

%% ---------------- Binary Search Initialization ----------------
% Initialize lower and upper bounds for the number of sensing beams.
feasible=false;
Llow = Nt;
Lup = Lmax;

%% ---------------- Binary Search ----------------
% Search for the minimum number of sensing beams satisfying
% the required worst-case JNR threshold.
while(abs(Llow-Lup)>delta)
    L = floor((Llow+Lup)/2);
    fprintf('Testing L = %d\n',L);
    Q = zeros(1,Z);
    for z = 1:Z

        %% Generate DFT Sensing Beam Directions
        % Construct uniformly distributed sensing beams over the angular region.
        theta_scan = zeros(1,L);
        for l = 1:L
            theta_scan(l) = asind(-1+(2*l-1)/L);
        end
        Y = zeros(Nr,L);
        X = zeros(Nt,L);
        at = steering_vector(Nt,true_theta);
        ar = steering_vector(Nr,true_theta);
        %% Received Signal Generation
        for l = 1 : L
            ns = sqrt(sigma2/2)*(randn(Nr,1) + 1j*randn(Nr,1));
            v = (1/sqrt(Nt))*steering_vector(Nt,theta_scan(l));
            xs = 1;
            Y(:,l) = sqrt(Ps)*ar*true_gain*(at.')*v*xs + ns;
            X(:,l) = sqrt(Ps/Nt)*steering_vector(Nt,theta_scan(l));
        end
        %% Covariance Matrix Computation
        Rxx=(1/L)*Y*(Y');
        Rx = (1/L)*X*(X');
        %% Direction of Arrival Estimation using Capon Beamforming Method
        theta_test = -90:0.1:90;
        P = zeros(size(theta_test));
        for p = 1 : length(theta_test)
            ar_test = steering_vector(Nr,theta_test(p));
            Rxx_reg = Rxx + 1e-6*eye(size(Rxx));
            P(p) = 1/real(ar_test'*(Rxx_reg\ar_test));
        end
        [~,idx] = max(P);
        theta_hat = theta_test(idx); 

        %% Channel Gain Estimation 
        gain = (ar'*Y*X'*conj(at))/(L*(ar'*ar)*(at'*conj(Rx)*at));   
        
        %% CRB Computation
        CRB = sigma2/(2*L*Ps*(abs(gain)^2));   
        % angular uncertainty interval.
        xi1 = theta_hat - 3*sqrt(CRB);
        xi2 = theta_hat + 3*sqrt(CRB);   

        %% Robust Beamformer Design
        thetha_mid = (xi1+xi2)/2;
        w = sqrt(Pmax/Nt)*steering_vector(Nt,thetha_mid);   
        beta = lambda*exp(1j*2*pi*dist/lambda)/(4*pi*dist);
        rho = abs(beta)^2/sigma2;
        
        %% Worst-Case JNR Evaluation
        at1 = steering_vector(Nt,xi1);
        JnrLeft = rho*abs((at1')*w)^2;  
        at2 = steering_vector(Nt,xi2);
        JnrRight = rho*abs((at2')*w)^2;
            
        F_worst = min(JnrLeft,JnrRight);
        Q(z) = F_worst;
    end

    %% Binary Search Decision
    Qmin = min(Q);
    fprintf('Estimated Gain = %.4f\n',abs(gain));
    fprintf('Least JNR = %.2f\n',Qmin);
    if(Qmin>=Sth)
        fprintf('Constraint satisfied.\n\n');
        Lup = L;
        feasible = true;
    else
        fprintf('Constraint NOT satisfied.\n\n');
        Llow = L;
    end
end
%% ---------------- Final Result ----------------
if(feasible == true)
    fprintf('Optimal Number of Sensing Beams = %d\n',Lup);
    figure;
    plot(theta_test,P,'LineWidth',2);
    xlabel('Angle(degrees)');
    ylabel('Capon Value');
    grid on;
else
    fprintf('We cannot find Optimal Number of Sensing Beams.');
end



function a = steering_vector(N,thetha)
thetha_rad = deg2rad(thetha);
a = zeros(N,1);
for m = 1 : N
    idx = 2*m-N-1;
    a(m) = exp(-1j*idx*pi*sin(thetha_rad)/2);
end
end
