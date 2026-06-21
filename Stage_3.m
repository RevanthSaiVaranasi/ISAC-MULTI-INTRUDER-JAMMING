%% ==========================================================
% Algorithm-3 : Intruder Scheduling Scheme
%
% Purpose:
%   Schedule intruders into feasible jamming groups.
%   If all intruders cannot be jammed simultaneously,
%   defer the weakest intruders and continue scheduling
%   until all intruders are successfully jammed.
%
% Outputs:
%   Scheduling Rounds
%   Jammed Groups
% ==========================================================


%% Load Detection Results
load('Algo1_Output.mat');
theta_all = theta_samples;
dist_all = dist;

%% Initialize Active and Deferred Sets
Remaining_Set = 1:size(theta_samples,1);
Deferred_Set = [];

Round = 1;

%% Scheduling Loop
while ~isempty(Remaining_Set)
    fprintf('\n====================================\n');
    fprintf('Scheduling Round %d\n',Round);
    fprintf('Active Intruders : ');
    disp(Remaining_Set);
    
    %% Build Current Active Intruder Set
    theta_samples = theta_all(Remaining_Set,:);
    dist = dist_all(Remaining_Set);

    while true
        %% Save active set for Algorithm-2
        save('Algo1_Output.mat','Nt','Sth','lambda','sigma2','Pmax',...
            'theta_samples','dist');

        %% Run Algorithm-2
        run('Algorithm_2.m');

        %% Load results
        load('Algo2_output.mat');

        tol = 1e-3;

        %% Feasible Check
        if min(JNR_intruder) >= (Sth - tol)
            fprintf('\nFeasible Jamming Group Found\n');
            fprintf('Jammed Intruders : ');
            disp(Remaining_Set);
            break;
        end

        %% Defer Weakest Intruder
        [~,idx_remove] = min(JNR_intruder);
        fprintf('Deferring Intruder %d\n',Remaining_Set(idx_remove));
        Deferred_Set = [Deferred_Set Remaining_Set(idx_remove)];

        %% Update Active Set
        Remaining_Set(idx_remove) = [];
        theta_samples(idx_remove,:) = [];
        dist(idx_remove) = [];
        if isempty(Remaining_Set)
            break;
        end
    end

    %% Current Group Successfully Jammed
    fprintf('\nCurrent Group Jammed Successfully\n');

    %% Move Deferred Intruders to Next Round
    Remaining_Set = Deferred_Set;
    Deferred_Set = [];
    Round = Round + 1;
end

fprintf('\n====================================\n');
fprintf('All intruders have been successfully jammed.\n');
fprintf('Total Scheduling Rounds = %d\n',Round-1);
fprintf('====================================\n');