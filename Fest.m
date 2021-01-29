function [ro, uo] = Fest(ahat,ssa,ssb,ssc,ssd,Sp,rr)
%Fine estimate
[rref, uref] = Cest(ahat,ssa,ssb,Sp);
rr = 1*(abs(rr)< 1 ) + rr;
phc = ahat(ssc,:);  phd = ahat(ssd,:);
d1 = Sp(ssc); d2 = Sp(ssd);
n1 = -ceil(d2 - d1):floor(d2 - d1); n2 = -ceil(d2 - d1):floor(d2 - d1);
[nn1, nn2] = meshgrid(n1,n2);
p1 = angle(phc)/(2*pi) + nn1; p2 = angle(phd)/(2*pi) + nn2;
ram = ( d1*p2.^2 - d2*p1.^2 + d2*d1^2 -  d1*d2^2 )./(2*p2*d1 - 2*p1*d2);
ram1 = ( d1*p2.^2 - d2*p1.^2 + d2*d1^2 -  d1*d2^2 )./(2*p2*d1 - 2*p1*d2);
uu1 =  (2*rr*p1 - p1.^2 + d1^2)./(2*rr*d1);
uu2 =  (2*rr*p2 - p2.^2 + d2^2)./(2*rr*d2);
row = find(abs(uu1(1,:)) < 1);
col = find(abs(uu2(:,1)) < 1);
ua = uu1(1,row);
ub = uu2(col,1)';
[inda] = (find(abs(ua - uref) == min(abs(ua - uref))));  nn(1) = row(inda) + n1(1) - 1;
[indb] = (find(abs(ub - uref) == min(abs(ub - uref))));  nn(2) = col(indb) + n2(1) - 1;
p1 = angle(phc)/(2*pi) + nn(1); p2 = angle(phd)/(2*pi) + nn(2);
ram = ( d1*p2.^2 - d2*p1.^2 + d2*d1^2 -  d1*d2^2 )./(2*p2*d1 - 2*p1*d2);
ramall = reshape(ram, 1, size(ram,1)*size(ram,2));
rrr = ramall(find(abs(ramall - rr) == min(abs(ramall - rr))));
ro = abs(rrr);
uo = (2*ro*p1 - p1.^2 + d1^2)./(2*ro*d1);
end

