function [dMdt,dMremin,dMfrag] = interactions(t,M,p,xz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M(M<0) = 0;
N = M./p.m(:);
m=p.m(:);
dM = zeros(size(N));
dMremin = zeros(size(dM));
dMfrag = zeros(size(dM));

%% Aggregation


%N = N(:);

for k = 1:length(p.bi)
    
    ii = p.bi(k)+1;
    jj = p.bj(k)+1;
    mi = m(ii);
    mj = m(jj);
    mij = m(ii)+m(jj);
    d00 = p.b300(k) + 1;
    d01 = p.b301(k) + 1;
    d10 = p.b310(k) + 1;
    d11 = p.b311(k) + 1;
    

    dN = p.alpha*p.beta(k)*N(ii)*N(jj);

    if dN > 0
        dM(ii) = dM(ii)-dN*mi;
        dM(jj) = dM(jj)-dN*mj;
        dM(d00) = dM(d00) + p.f00(k)*dN*mij; 
        dM(d01) = dM(d01) + p.f01(k)*dN*mij; 
        dM(d10) = dM(d10) + p.f10(k)*dN*mij; 
        dM(d11) = dM(d11) + p.f11(k)*dN*mij;
    end
end

%% Degradation
if p.remin~=0
for i = 1:length(N)
    tmp = xz(i-1);
    xtmp = tmp(1);
    ztmp = tmp(2);
    if mod(i,p.nD)==0
        dN_remin = p.remin*p.q^xtmp *((1-p.q)/(3-p.a))*(-ztmp*N(i));
        
    else
        dN_remin = p.remin*p.q^xtmp *((1-p.q)/(3-p.a))*((ztmp+1)*N(i+1)-ztmp*N(i));
    end
    dMremin(i) = dN_remin*m(i);
end
end
%% Fragmentation
for i  = 1:length(M)
    if i< p.nD+1
        dMfrag(i) = p.pfrag(i+p.nD).*M(i+p.nD)*(p.frag_div);
    elseif i>length(M)-p.nD
        dMfrag(i) = -p.pfrag(i).*M(i)*(p.frag_div);
    else
         dMfrag(i) = p.pfrag(i+p.nD)*M(i+p.nD)*p.frag_div - p.pfrag(i)*M(i)*p.frag_div;
    end
end
 

%% collecting

dMdt = dM +p.prod(:) - M.*p.wWhites(:)./p.H +dMremin +dMfrag;


end

