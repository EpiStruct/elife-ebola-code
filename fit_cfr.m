%% Script to infer the Case fatality ratio parameters

clear

%% Import data

xx = importdata('ebola_short.csv');
cs=xx(:,2);
de=xx(:,3);

%% Main MCMC loop

L = 1.01e6;% Number of unthinned samples to take
params = zeros(L,2);
% Starting parameters:
params(1,1) = 7;
params(1,2) = 7;
LL = zeros(L,1);
LL(1) = 0;
for i=1:length(cs)
    LL(1) = LL(1) ...
        + log(beta(params(1,1)+de(i),params(1,2)+cs(i)-de(i))) ...
        - log(beta(params(1,1),params(1,2)));
end
vars = [2, 1]; % Proposal variances
for l=2:L
    params(l,:) = params(l-1,:);
    LL(l) = LL(l-1);
    for w=1:2
        pprop = params(l,:);
        pprop(w) = pprop(w) + vars(w)*randn;
        if (pprop(w) >0) && (pprop(w) <100)
            Lprop = 0;
            for i=1:length(cs)
                Lprop = Lprop ...
                    + log(beta(pprop(1)+de(i),pprop(2)+cs(i)-de(i))) ...
                    - log(beta(pprop(1),pprop(2)));
            end
            if (log(rand) < (Lprop - LL(l)))
                LL(l) = Lprop;
                params(l,:) = pprop;
            end
        end
    end
end

%% Save poserior samples

B = 1e4; % Burn-in period
T = 1e3; % Thinning
dlmwrite('cfr_out.txt',params(B:T:end,:));


