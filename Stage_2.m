%% ==========================================================
% Algorithm-2 : Robust Jamming Beamforming Optimization
%
% Purpose:
%   Design a transmit beamformer that satisfies
%   worst-case JNR constraints for all detected intruders.
%
% Outputs:
%   JNR_intruder
%   Power
%   Feasible
%   Algo2_output.mat
% ==========================================================

%% Load Detection Results and User Defined Inputs
load('Algo1_Output.mat');
[N,Qtheta] = size(theta_samples);
Iterations = input('Enter the number of times iteration should be: ');
eps_conv = input('Enter the converging criteria to stop: ');

%% Path Loss Computation
rho_intruder = zeros(N,1);
for k = 1:N
    beta_k = lambda*exp(1j*2*pi*dist(k)/lambda)/(4*pi*dist(k));
    rho_intruder(k) = abs(beta_k)^2/sigma2;
end


%% Main Optimization Loop
rng(3);
w = randn(Nt,1) + 1j*randn(Nt,1);
F_prev = real(w' * w);
n = 0;

%% Construct Steering Matrices    
while n < Iterations
    n = n + 1;
    %% Update Auxiliary Variables eta
    A = cell(N,Qtheta);
    for k = 1:N
        for q = 1:Qtheta
            theta = theta_samples(k,q);
            alpha = steering_vector(Nt,theta);
            A{k,q} = alpha*alpha';
        end
    end
    %%Iterative solvation for eta
    eta = ones(N,Qtheta);
    %% Beamformer Update
    for iter_eta = 1:20
        eta_old = eta;
        C = zeros(N,Qtheta);
        for k = 1:N
            for q = 1:Qtheta
                C(k,q) = (real(w'*A{k,q}*w) + Sth/rho_intruder(k))/2;
            end
        end
        for k = 1:N
            for q = 1:Qtheta
                bk = A{k,q}*w;
                numerator = C(k,q);
                for i = 1:N
                    for j = 1:Qtheta
                        if ~(i==k && j==q)
                            bi = A{i,j}*w;
                            numerator = numerator - eta(i,j)*real(bi'*bk);
                        end
                    end
                end
                denominator = real(bk'*bk);
                eta(k,q) = max(numerator/denominator,0);
            end
        end
        if norm(eta-eta_old,'fro') < 1e-4
            break;
        end
    end
    %% Power Constraint Enforcement
    w_new = zeros(Nt,1);
    for k = 1:N
        for q = 1:Qtheta
            w_new = w_new + eta(k,q)*A{k,q}*w;
        end
    end
    %% JNR Evaluation
    if real(w_new'*w_new) > Pmax
        w_new = sqrt(Pmax/(real(w_new'*w_new))) * w_new;
    end
    F_new = real(w_new' * w_new);
    if abs(F_new - F_prev) < eps_conv
        break;
    end
    w = w_new;
    F_prev = F_new;
end

%% Feasibility Check
JNR = zeros(N,Qtheta);

for k = 1:N

    rho_k = rho_intruder(k);

    for q = 1:Qtheta

        JNR(k,q) = rho_k*real(w'*A{k,q}*w);

    end

end

%% Worst-case JNR of each intruder
JNR_intruder = min(JNR,[],2);

%% Final transmit power
Power = real(w'*w);

%% Feasibility check
tol = 1e-3;


fprintf('Algorithm-2 Result\n');

fprintf('Total Transmit Power = %.4f\n',Power);
fprintf('Required JNR Threshold = %.4f\n',Sth);

fprintf('\nWorst Case JNR of Each Intruder:\n');
disp(JNR_intruder);

if min(JNR_intruder) >= (Sth - tol)

    fprintf('\nBeamforming Solution Feasible\n');
    fprintf('All intruders can be intercepted.\n');

    Feasible = 1;

else

    fprintf('\nBeamforming Solution Infeasible\n');
    fprintf('Scheduling is required.\n');

    Feasible = 0;

end

%% Save Results for Stage-3
save('Algo2_output.mat','Power','JNR_intruder','Feasible');

%% Beam Pattern Visualization
theta_plot = -90:0.1:90;
BeamPattern = zeros(size(theta_plot));

for p = 1:length(theta_plot)
    a = steering_vector(Nt,theta_plot(p));
    BeamPattern(p) = abs(a'*w)^2;
end

BeamPattern = BeamPattern/max(BeamPattern);

figure;
plot(theta_plot,10*log10(BeamPattern),'LineWidth',2);
xlabel('Angle (degrees)');
ylabel('Normalized Beam Power (dB)');
title('Optimized Beam Pattern');
grid on;

