%%  Performance of MILE for the variant of correlation coefficient rho
clear
close all
TrialAll = 500;
count = 0;
for co = 0:0.1:0.9
    disp(sprintf('co: %.1f', co));
    SnrdB = 20;
    count = count + 1;
    for trial = 1:TrialAll;
        Sensor = 9;
        Snap = 1000;
        Doa = [10 30];
        Rng = [3 2000];
        SigAll = length(Doa); NfAll = length(Rng);
        Lambda = 1;  D = 0.5*Lambda;
        Sp = [-2 -1.4 -1.1 -0.5 0 0.5 1.2 1.6 2];
        for num = 1:NfAll
            rngall  = sqrt(Rng(num).^2 + Sp'.^2 - 2*Rng(num).*Sp'*sind(Doa(num)));
            tall = asin((rngall*sind(Doa(num)) - Sp')./(sqrt(rngall.^2 + (Sp').^2 - 2*rngall.*Sp'.*sind(Doa(num)))));
            A(:,num) = Rng(num)./rngall .*exp(-j *2*pi * (Rng(num) - rngall));
        end
        %Af = exp(-j*2*pi*Sp'*sind(Doa(2:3)));
        Snr=sqrt(10.^(SnrdB/10));
        f = [0.2, 0.4, -0.12,-0.31];
        for num = 1:SigAll
            Sig(num,:) = fmlin(Snap,f(num), f(num));
        end
        if co ~=0
            Sig(2,1:co*Snap) = Sig(1,1:co*Snap);
        end
        
        Noise = sqrt(0.5)* (randn(Sensor,Snap) + j* randn(Sensor,Snap));
        G = diag(randn(1,Sensor)*0.122 + 1);
        X = Snr* G* A * Sig(1:NfAll,:)   + 1*Noise;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Smd = (Sensor + 1)/2;        ssa = Smd - 1; ssb = Smd+ 1;   ssc = Smd - 4;  ssd = Smd +4;
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
        end   %% for num
    end   %% for trial
    
    if Doa(1) < Doa(2)
        for trial = 1:TrialAll
            rpha(trial,:) = sort(abs(rpha(trial,:)));        upha(trial,:) =  asind(sort(upha(trial,:)));
        end
    else
        for trial = 1:TrialAll
            rpha(trial,:) = sort(abs(rpha(trial,:)));        upha(trial,:) =  asind(sort(upha(trial,:),'descend'));
        end
    end
    
    Errrpha(:, count)= rmse(rpha,Rng);
    Errupha(:, count) = rmse(upha,Doa);
end

D = 0:0.1:0.9;
figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(D,Errupha(1,:), '-o', D,Errupha(2,:), '-s'); grid on;   %%% FF 
xlabel('Correlation coefficient'); ylabel('RMSE of angle estimates, in degrees')
title('(a): Performance of angle estimates')

figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(D,Errrpha(1,:), '-d'); grid on;   %%% FF 
xlabel('Correlation coefficient'); ylabel('RMSE of range r_1 estimates, in wavelength')
title('(b): Performance of range estimates')

