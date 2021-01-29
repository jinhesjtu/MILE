%% Classfication Success-Rate
clear
close all
TrialAll = 500; count = 0;
for SnrdB = 0:2:40
    disp(sprintf('SnrdB: %.1f', SnrdB));
    count = count + 1;
    suc1 = 0; suc2 = 0; suc3 = 0;
    for trial = 1:TrialAll;
        Sensor = 9;
        Snap = 1000;
        Doa = [10 30 32];
        Rng = [3 2000 3000];
        SigAll = length(Doa); NfAll = length(Rng);
        Lambda = 1;  D = 0.5*Lambda;
        Sp = [-2 -1.4 -1.1 -0.5 0 0.5 1.2 1.6 2]; 
        Sp2 = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
        for num = 1:NfAll
            rngall  = sqrt(Rng(num).^2 + Sp'.^2 - 2*Rng(num).*Sp'*sind(Doa(num)));
            tall = asin((rngall*sind(Doa(num)) - Sp')./(sqrt(rngall.^2 + (Sp').^2 - 2*rngall.*Sp'.*sind(Doa(num)))));
            A(:,num) = Rng(num)./rngall .*exp(-j *2*pi * (Rng(num) - rngall));
        end
        for num = 1:NfAll
            rngall  = sqrt(Rng(num).^2 + Sp2'.^2 - 2*Rng(num).*Sp2'*sind(Doa(num)));
            tall = asin((rngall*sind(Doa(num)) - Sp2')./(sqrt(rngall.^2 + (Sp2').^2 - 2*rngall.*Sp2'.*sind(Doa(num)))));
            A2(:,num) = Rng(num)./rngall .*exp(-j *2*pi * (Rng(num) - rngall));
        end
        
        Snr=sqrt(10.^(SnrdB/10));
        f = [0.2, 0.4, -0.12, -0.31];
        for num = 1:SigAll
            Sig(num,:) = fmlin(Snap,f(num), f(num));
        end
        Noise = sqrt(0.5)* (randn(Sensor,Snap) + j* randn(Sensor,Snap));
        G = diag(randn(1,Sensor)*0.122 + 1);
        X = A * diag(Snr)* Sig(1:NfAll,:)   + 1*Noise;
        X2 = A2 * diag(Snr)* Sig(1:NfAll,:)   + 1*Noise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Smd = (Sensor + 1)/2; ssa = Smd - 1; ssb = Smd+ 1; ssc = Smd - 4; ssd = Smd +4;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        for nc = 1:Sensor
            C(:,:,nc) = cum4mtx(X(nc,:), X(Smd,:), X,X);
        end
        [A1,A2,A3]=comfac(C,SigAll);
        Ahat = A1 + conj(A2) + A3;
        Ahat = norml(A2,Smd);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Classification Before Estimation
        for num = 1:SigAll
            %%%%%%%%%    Coarse Estimation
            c = Ahat(:,num);
            ahat = c/c((Sensor+1)/2);
            %    rcoa(trial, num) = db(abs(ahat(1)));
            [rref, uref] = Cest(ahat,ssa,ssb,Sp);     rr = abs(rref);
            rcoa(trial, num) = abs(rref);  ucoa(trial, num) = uref;
            %%%%%%%%%    Coarse Estimation
            
             [rhat, uhat] = Fest(ahat,ssa,ssb,ssc,ssd,Sp,rr);
             rpha(trial, num) = rhat; upha(trial, num) = uhat;
        end
        rest(trial,:) = sort(rcoa(trial,:));
        rfin(trial,:) = sort(rpha(trial,:));    
        DR = 6;
        a = find(rest(trial,:) < DR);
        if length(a) == 1
            suc1 = suc1 + 1;
        end
        a = find(rfin(trial,:) < DR);
        if length(a) == 1
            suc2 = suc2 + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Classification Before Estimation
       
    end %% for trial
    Suc(:, count) = [suc1;suc2];
end
Suc = Suc/TrialAll;

SnrdB = 0:2:40;
D = SnrdB;
figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(D,Suc(1,:), '-o', D,Suc(2,:), '-*'); grid on;
