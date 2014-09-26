%% Script to infer the arrival time of new Ebloa outbreaks

clear

%% Import data

xx = importdata('ebola_short.csv');
yr = xx(:,1);
t=max(yr)-min(yr)+1;
k = length(yr);

%% Main MCMC loop

L = 1.01e6; % Number of unthinned samples to take
params = zeros(L,1);
params(1,1) = k/t; % Start of chain
LL = zeros(L,1);
la = params(1,1);
LL(1) = (k*log(la*t)) - (la*t);
vars = [0.2]; % Proposal variance
for l=2:L
    params(l,:) = params(l-1,:);
    LL(l) = LL(l-1);
    for w=1
        pprop = params(l,:);
        pprop(w) = pprop(w) + vars(w)*randn;
        if (pprop(w) >0) && (pprop(w) <100)
            la = pprop(1);
            Lprop = (k*log(la*t)) - (la*t);
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
dlmwrite('la_out.txt',params(B:T:end,:));


