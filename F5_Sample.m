%% Classfication and Localization 2-NF 2-FF overestimate K, K = 5
clear
close all
TrialAll = 500;
for trial = 1:TrialAll;
    Sensor = 9;
    Snap = 1000;
    Doa = [10 20 30 40];
    Rng = [3 6  2000 3000];
    SigAll = length(Doa); NfAll = length(Rng);
    Lambda = 1;  D = 0.5*Lambda;
    Sp = [-2 -1.4 -1.1 -0.5 0 0.5 1.2 1.6 2];
    for num = 1:NfAll
        rngall  = sqrt(Rng(num).^2 + Sp'.^2 - 2*Rng(num).*Sp'*sind(Doa(num)));
        tall = asin((rngall*sind(Doa(num)) - Sp')./(sqrt(rngall.^2 + (Sp').^2 - 2*rngall.*Sp'.*sind(Doa(num)))));
        A(:,num) = Rng(num)./rngall .*exp(-j *2*pi * (Rng(num) - rngall));
    end
    SnrdB = [19 21 22 25];
    Snr=sqrt(10.^(SnrdB/10));
    f = [0.2, 0.4, -0.12, -0.31];
    for num = 1:SigAll
        Sig(num,:) = fmlin(Snap,f(num), f(num));
    end
    Noise = sqrt(0.5)* (randn(Sensor,Snap) + j* randn(Sensor,Snap));
    G = diag(randn(1,Sensor)*0.122 + 1);  %% magnitude error
    X = G* A * diag(Snr)* Sig(1:NfAll,:)   + 1*Noise;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Smd = (Sensor + 1)/2; ssa = Smd - 1; ssb = Smd+ 1; ssc = Smd - 4; ssd = Smd +4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for nc = 1:Sensor
        C(:,:,nc) = cum4mtx(X(nc,:), X(Smd,:), X,X);
    end
    SigAll = SigAll + 1;
    [A1,A2,A3]=comfac(C,SigAll);
    Ahat = A1 + conj(A2) + A3;
    Ahat = norml(A2,Smd);
    R = X*X'/Snap;  [U S V] = svd(R); En = U(:,SigAll + 1:end);
    
    % %%%%%%%   Detection Nf - FF
    
    for num = 1:SigAll
        %%%%%%%%%    Coarse Estimation
        c = Ahat(:,num);
        ahat = c/c((Sensor+1)/2);
        [rref, uref] = Cest(ahat,ssa,ssb,Sp);     rr = abs(rref);
        rcoa(trial, num) = abs(rref);  ucoa(trial, num) = uref;
        %%%%%%%%%    Coarse Estimation
        
        %%%%%%%%%               Fine Estimation
        [rhat, uhat] = Fest(ahat,ssa,ssb,ssc,ssd,Sp,rr);
        rpha(trial, num) = rhat; upha(trial, num) = uhat;
        %%%%%%%%%                Fine Estimation
    end
    rini(trial,:) = rcoa(trial,:);  uini(trial,:) = ucoa(trial,:);
    
    rfin(trial,:) = rpha(trial,:); ufin(trial,:) = upha(trial,:);
    
end %% for trial
rini = abs(rini).*(2*(abs(rini) >10)) + abs(rini);

for trial = 1:TrialAll
    rini(trial, :) = (rini(trial,:));
    uini(trial, :) = asind((uini(trial,:)));
    rfin(trial, :) = (rfin(trial,:));
    ufin(trial, :) = (asind(ufin(trial,:)));
end

for trial = 1:TrialAll
    temp = ufin(trial, :);
    for num = 1:SigAll - 1;
        [b a] = min(abs(temp - Doa(num)));
        ind(num) = a; temp(a) = 1000;
        rinib(trial,num) = rini(trial,a);  ufinb(trial,num) = ufin(trial,a);
    end
    m = 15 - sum(ind);
    rinib(trial,5) = rini(trial,m);  ufinb(trial,5) = ufin(trial,m);
end
DR = 32;
N = 1:TrialAll;
figure; set(gcf,'DefaultLineLineWidth',1.5)
semilogy(N,rinib(:,1),'.', N,rinib(:,2), '.', N,rinib(:,3), '.', N,rinib(:,4), '.', N,rinib(:,5), '*', N,DR*ones(1,TrialAll),'-k'); grid on;
axis([1,TrialAll, 1,1e5])
xlabel('Number of trials')
ylabel('Range, in wavelength')
legend('NF source 1', 'NF source 2', 'FF source 1', 'FF source 2', 'Nonexistent source')


N = 1:TrialAll;
figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(N,ufinb(:,1), '.', N,ufinb(:,2), '.', N,ufinb(:,3), '.', N,ufinb(:,4), '.', N,ufinb(:,5), '*'); grid on;
xlabel('Number of trials')
ylabel('Angle, in degrees')
legend('NF source 1', 'NF source 2', 'FF source 1', 'FF source 2','Nonexistent source')
