%% Script to infer the arrival time of new Ebloa outbreaks

clear

%% Import data

xx = importdata('ebola_short.csv');
cs=xx(:,2);
m = max(cs);
[yyy, xx] = hist(cs,1:m);
yy = yyy/sum(yyy);
jj=find(yyy);

%% Main MCMC loop

L = 1.01e4;
np = 2;
params = zeros(L,np);
params(1,1) = 0.8;
params(1,2) = 0.1;
LL = zeros(L,1);
R0 = params(1,1);
rr = 1;
pp = rr/(R0+rr);
im = params(1,2);
Z = geopdf(0:m,im);
Y = zeros(1,m);
for a=1:(m+1)
    ii = jj-a;
    for k=ii(ii>=0)
        Y(k+a) = Y(k+a)+(Z(a)*(a/(k+a))*nbinpdf(k,(k+a)*rr,pp));
    end
end
LL(1) = sum(yyy(jj).*log(Y(jj)));
vars = [0.05, 0.05];
pl = [0, 0];
pu = [2, 1];

for l=2:L
    params(l,:) = params(l-1,:);
    LL(l) = LL(l-1);
    for w=1:np
        pprop = params(l,:);
        pprop(w) = pprop(w) + vars(w)*randn;
        R0 = pprop(1);
        rr = 1;
        pp = rr/(R0+rr);
        im = pprop(2);
        Z = geopdf(0:m,im);
        Y = zeros(1,m);
        for a=1:(m+1)
            ii = jj-a;
            for k=ii(ii>=0)
                Y(k+a) = Y(k+a)+(Z(a)*(a/(k+a))*nbinpdf(k,(k+a)*rr,pp));
            end
        end
        Lprop = sum(yyy(jj).*log(Y(jj)));
        if (log(rand) < (Lprop - LL(l))) && (pprop(w) > pl(w)) && (pprop(w) < pu(w))
            LL(l) = Lprop;
            params(l,:) = pprop;
        end
    end
end

%% Save posterior samples
B = 1e2;
T = 1e1;
dlmwrite('bp_out.txt',params(B:T:end,:));


