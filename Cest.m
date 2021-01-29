function [ro, uo] = Cest(ahat,ssa,ssb,Sp)
%Coarse estimate
pha = ahat(ssa,:);  phb = ahat(ssb,:);
p1 = angle(pha)/(2*pi); p2 = angle(phb)/(2*pi);
d1 = Sp(ssa); d2 = Sp(ssb);
ro = ( d1*p2.^2 - d2*p1.^2 + d2*d1^2 -  d1*d2^2 )./(2*p2*d1 - 2*p1*d2);
uo = (2*ro*p1 - p1.^2 + d1^2)./(2*ro*d1);
end

