clear

xx = importdata('ebola_short.csv');
cs=xx(:,2);
m = max(cs); % Largest in data
M1 = 1e3; % Minimum size for ongoing outbreak
M2 = 5e3; % Maximum size for ongoing outbreak
[yyy, xx] = hist(cs,1:m);
yy = yyy/sum(yyy);
jj=find(yyy);

%%

L = 1.01e4;
np = 2;
params = zeros(L,np);
params(1,1) = 0.95;
params(1,2) = 0.1;
LL = zeros(L,1);
CC = zeros(L,1);
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
CC(1) = largeC(M1,M2,im,pp);
LL(1) = sum(yyy(jj).*log(Y(jj)));
vars = [0.035, 0.05];
pl = [0, 0];
pu = [2, 1];

for l=2:L
    params(l,:) = params(l-1,:);
    CC(l) = CC(l-1);
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
        C = largeC(M1,M2,im,pp);
        Lprop = sum(yyy(jj).*log(Y(jj))) + log(C);
        if (log(rand) < (Lprop - LL(l))) && (pprop(w) > pl(w)) && (pprop(w) < pu(w))
            LL(l) = Lprop;
            params(l,:) = pprop;
            CC(l) = C;
        end
    end

end

%% Save posterior samples
B = 1e2;
T = 1e1;
dlmwrite('bp_aug_out.txt',params(B:T:end,:));
dlmwrite('CC.txt',CC(B:T:end,:));


