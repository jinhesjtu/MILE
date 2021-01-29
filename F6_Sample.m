%%  FF source
clear
close all
TrialAll = 500;
count = 0;
for SnrdB = 0:5:40
    disp(sprintf('SnrdB: %.1f', SnrdB));
    count = count + 1;
    for trial = 1:TrialAll;
        Sensor = 9;
        Snap = 1000;
        Doa = [10 30];
        Rng = [3 2000];
        SigAll = length(Doa); NfAll = length(Rng);
        Lambda = 1;  D = 0.25*Lambda;
        Sp = [-(Sensor-1)/2 : (Sensor-1)/2]*D;
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
        X = Snr* A * Sig(1:NfAll,:)   + 1*Noise;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Smd = (Sensor + 1)/2;        ssa = Smd - 2; ssb = Smd+ 2;   ssc = Smd - 4;  ssd = Smd +4;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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
            
            %%%%%%%%%               Fine Estimation
            [rhat, uhat] = Fest(ahat,ssa,ssb,ssc,ssd,Sp,rr);
            rpha(trial, num) = rhat; upha(trial, num) = uhat;
            %%%%%%%%%                Fine Estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refine
            uini(num) = round(1000*asind(upha(trial,num)))/1000; rini(num) = round(1000*abs(rpha(trial,num)))/1000;
            rstep = 0.01; ustep = 0.01;
            if rini(num) < 10;
                nr = 0;
                for dr = rini(num) - 1 :rstep: rini(num) + 1
                    nr = nr + 1; nu = 0;
                    for du = uini(num) - 1 : ustep: uini(num) + 1
                        nu = nu + 1;
                        r1  = sqrt(dr.^2 + Sp'.^2 - 2*dr.*Sp'*sind(du));
                        ar1 = dr./r1 .*exp(-j *2*pi * (r1 - dr));
                        Pr(nr,nu) = abs(1/(ar1'*En*En'*ar1));
                    end
                end
                [b1 a1] = max(max(Pr)); rmu(num) = (a1 -1)*rstep + rini(num) - 1;
                [b1 a1] = max(max(Pr.')); umu(num) = (a1 -1)*ustep + uini(num) - 1;
            else
                dr = 5000; nu = 0;
                for du = uini(num) - 1 : ustep: uini(num) + 1
                    nu = nu + 1;
                    r1  = sqrt(dr.^2 + Sp'.^2 - 2*dr.*Sp'*sind(du));
                    ar1 = dr./r1 .*exp(-j *2*pi * (r1 - dr));
                    Prf(nu) = abs(1/(ar1'*En*En'*ar1));
                end
                [b1 a1] = max(Prf); rmu(num) = dr;
                umu(num) = (a1 -1)*ustep + uini(num) - 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refine            
        end   %% for num
        rrefi(trial,:) = sort(rmu); urefi(trial,:) = sort(umu);
    end   %% for trial
    for trial = 1:TrialAll
        rpha(trial,:) = sort(abs(rpha(trial,:)));        upha(trial,:) =  asind(sort(upha(trial,:)));
        rrefi(trial,:) = sort(abs(rpha(trial,:)));        urefi(trial,:) =  (sort(urefi(trial,:)));
    end
    Errrpha(:, count)= rmse(rpha,Rng);
    Errupha(:, count) = rmse(upha,Doa);
    Errrfefi(:, count)= rmse(rrefi,Rng);
    Errurefi(:, count) = rmse(urefi,Doa);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRB FF
end
count = 0;
for SnrdB = 0:5:40
    count = count + 1;
    syms soru sorr;
    tu = Doa(2)*pi/180;
    rv = [2000];
    Snr=sqrt(10.^(SnrdB/10));
    rall  = sqrt(sorr.^2 + Sp'.^2 - 2*sorr.*Sp'*sin(soru));  rl = sqrt(sorr.^2 + D^2 - 2*sorr.*D*sin(soru));
    A = sorr./rall.*exp(-j*2*pi/Lambda * (rall - sorr));
    Ju = double(subs(diff(A,'soru'), {soru,sorr},{tu,rv}));
    Jr = double(subs(diff(A,'sorr'), {soru,sorr},{tu,rv}));
    J = [Ju'*Ju,Ju'*Jr;
        Jr'*Ju,Jr'*Jr];
    JJ = double(J);JJ = real(JJ);JJ = inv(JJ);
    crbu = asind(sqrt((JJ(1,1)/(2*Snap*(Snr^2)))));
    crbr =  sqrt((JJ(2,2)/(2*Snap*(Snr^2))));
    Crbu(count) = crbu;
    Crbr(count) = crbr;
end

SnrdB = 0:5:40;
D = SnrdB;
figure; set(gcf,'DefaultLineLineWidth',1.5)
semilogy(D,Errurefi(2,:), '-x', D,Crbu, '-k'); grid on;   %%% FF
% axis([50 10000 1e-3 1])
xlabel('SNR, in dB'); ylabel('RMSE of angle \theta_2 estimates, in degrees')

