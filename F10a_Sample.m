%% Identification results for L = 5, K = 4;
clear
close all
TrialAll = 500;
for trial = 1:TrialAll;
    Sensor = 5;
    Snap = 1000;
    Doa = [10 20 30 40];
    Rng = [3 6  2000 3000];
    SigAll = length(Doa); NfAll = length(Rng);
    Lambda = 1;  D = 0.5*Lambda;
    Sp = [-2   -0.5 0 0.5   2];
    for num = 1:NfAll
        rngall  = sqrt(Rng(num).^2 + Sp'.^2 - 2*Rng(num).*Sp'*sind(Doa(num)));
        tall = asin((rngall*sind(Doa(num)) - Sp')./(sqrt(rngall.^2 + (Sp').^2 - 2*rngall.*Sp'.*sind(Doa(num)))));
        A(:,num) = Rng(num)./rngall .*exp(-j *2*pi * (Rng(num) - rngall));
    end
    SnrdB = 20;
    Snr=sqrt(10.^(SnrdB/10));
    f = [0.2, 0.4, -0.12, -0.31, 0.43, 0.15];
    for num = 1:SigAll
        Sig(num,:) = fmlin(Snap,f(num), f(num));
    end
    Noise = sqrt(0.5)* (randn(Sensor,Snap) + j* randn(Sensor,Snap));
    G = diag(randn(1,Sensor)*0.122 + 1);
    X = Snr *G* A *  Sig(1:NfAll,:)   + 1*Noise;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Smd = (Sensor + 1)/2; ssa = Smd - 1; ssb = Smd+ 1; ssc = Smd - 2; ssd = Smd +2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for nc = 1:Sensor
        C(:,:,nc) = cum4mtx(X(nc,:), X(Smd,:), X,X);
    end
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
for trial = 1:TrialAll
    rini(trial, :) = sort(rini(trial,:));
    uini(trial, :) = asind(sort(uini(trial,:)));
    rfin(trial, :) = sort(rfin(trial,:));
    ufin(trial, :) = sort(asind(ufin(trial,:)));
end

figure; set(gcf,'DefaultLineLineWidth',1.5)
semilogy(ufin(:,1), rfin(:,1), '*', ufin(:,2), rfin(:,2), '*', ufin(:,3), 2000*ones(TrialAll,1), '*', ...
    ufin(:,4), 3000*ones(TrialAll,1), '*'); grid on
title('(a)')
axis([-20,45,1,10000])
xlabel('Angle, in degree')
ylabel('Range, in wavelength')
legend('Source 1', 'Source 2', 'Source 3', 'Source 4', 'Location', 'NorthWest')

Errr= rmse(rfin,Rng);
Erru = rmse(ufin,Doa);

[Errr; Erru]