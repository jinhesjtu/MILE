%%  Performance versus sensor location errors
clear
close all
TrialAll = 500;
count = 0;
for Le = [0:0.002:0.1];
    disp(sprintf('Le: %.1f', Le*1000));
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
        Spe = randn(1,Sensor)*Le;
        Sp = Sp + Spe;
        clear A C
        for num = 1:NfAll
            rngall  = sqrt(Rng(num).^2 + Sp'.^2 - 2*Rng(num).*Sp'*sind(Doa(num)));
            tall = asin((rngall*sind(Doa(num)) - Sp')./(sqrt(rngall.^2 + (Sp').^2 - 2*rngall.*Sp'.*sind(Doa(num)))));
            A(:,num) =Rng(num)./rngall .*exp(j *2*pi * ( Rng(num) - rngall));
        end
        %Af = exp(-j*2*pi*Sp'*sind(Doa(2:3)));
        Snr=sqrt(10.^(SnrdB/10));
        f = [0.2, 0.4, -0.12,-0.31];
        for num = 1:SigAll
            Sig(num,:) = fmlin(Snap,f(num), f(num));
        end
        Noise = sqrt(0.5)* (randn(Sensor,Snap) + j* randn(Sensor,Snap));
        G = diag(randn(1,Sensor)*0.122 + 1);
        X = Snr* G*A * Sig(1:NfAll,:)   + 1*Noise;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Smd = (Sensor + 1)/2;        ssa = Smd - 1; ssb = Smd + 1;   ssc = max(Smd - 4,1);  ssd = min(Sensor, Smd + 4);
        ssc = 1; ssd = Sensor;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        Sp = [-2 -1.4 -1.1 -0.5 0 0.5 1.2 1.6 2];
        for nc = 1:Sensor
            C(:,:,nc) = cum4mtx(X(nc,:), X(Smd,:), X,X);
        end
        [A1,A2,A3]=comfac(C,SigAll);
        Ahat = A1 + conj(A2) + A3;
        Ahat = norml(A3,Smd);
        R = X*X'/Snap;  [U S V] = svd(R); En = U(:,SigAll + 1:end);
        
        % %%%%%%%   Detection Nf - FF        
        
        for num = 1:SigAll
            %%%%%%%%%    Coarse Estimation
            c = Ahat(:,num);
            ahat = c/c((Sensor+1)/2);
            [rref, uref] = Cest(ahat,ssa,ssb,Sp);     rr = abs(rref);
            rcoa(trial, num) = abs(rref);  ucoa(trial, num) = uref;
            %%%%%%%%%    Coarse Estimation            
            %%%%%%%%%    Fine Estimation
            [rhat, uhat] = Fest(ahat,ssa,ssb,ssc,ssd,Sp,rr);
            rpha(trial, num) = rhat; upha(trial, num) = uhat;
            %%%%%%%%%                Fine Estimation            
        end   %% for num
    end   %% for trial
    
    for trial = 1:TrialAll
        rpha(trial,:) = sort(abs(rpha(trial,:)));        
        upha(trial,:) =  asind(sort(upha(trial,:)));
    end    
    Errrpha(:, count)= rmse(rpha,Rng);
    Errupha(:, count) = rmse(upha,Doa);
end

syms u r 
count = 0;
for    Le = [0:0.002:0.1];
    SnrdB = 20;
    Snap = 1000;
    Sensor = 9;
    DD = 0.5;        count = count + 1;
    tu = d2r([10 30]);
    rv = [3 2000];           SigAll = length(rv);
    Lambda = 1;  D = 0.5*Lambda;
        Sp = [-2 -1.4 -1.1 -0.5 0 0.5 1.2 1.6 2];
        Spe = randn(1,Sensor)*Le;
        Sp = Sp + Spe;
        clear A C
    Snr=sqrt(10.^(SnrdB/10)); % D = 0.5*Lambda;    
    clear A Ju Jr J JJ
    rall  = sqrt(r.^2 + Sp'.^2 - 2*r.*Sp'*sin(u)); % rl = sqrt(r.^2 + D^2 - 2*r.*D*sin(u));
    A = (r./rall).*exp(-j*2*pi/Lambda * (rall - r));    
    for num = 1:SigAll
        Ju = double(subs(diff(A,'u'), {u,r},{tu(num),rv(num)}));
        Jr = double(subs(diff(A,'r'), {u,r},{tu(num),rv(num)}));
        J = [Ju'*Ju,Ju'*Jr;
            Jr'*Ju,Jr'*Jr];
        JJ = double(J);
        JJ = real(JJ);
        Jf = diag(inv(JJ));
        crbu(num) = asind(sqrt((Jf(1)/(2*Snap*(Snr^2)))));
        crbr(num) =  sqrt((Jf(2)/(2*Snap*(Snr^2))));
    end
    Crbu(:,count) = crbu;
    Crbr(:,count) = crbr;
end

D = 0:0.002:0.1;
figure; set(gcf,'DefaultLineLineWidth',1.5)
semilogy(D,Errupha(1,:), '-o', D,Errupha(2,:), '-s', D, Crbu(1,:), '-k', D, Crbu(2,:), '-b'); grid on;   %%% FF 
xlabel('Standard deviation of location error'); ylabel('RMSE of angle estimates, in degrees')
title('(a): Performance of angle estimates')

figure; set(gcf,'DefaultLineLineWidth',1.5)
semilogy(D,Errrpha(1,:), '-o',  D,Crbr(1,:), '-k'); grid on; 
xlabel('Standard deviation of location error'); ylabel('RMSE of range r_1 estimates, in wavelength')
title('(b): Performance of range estimates')
