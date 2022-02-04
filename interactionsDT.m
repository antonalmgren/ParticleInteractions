function [dMdt,dMremin,dMfrag] = interactionsDT(t,M,m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M(M<0) = 0;
N = M./m(:);
m=m(:);
dM = zeros(size(N));
dMremin = zeros(size(dM));
dMfrag = zeros(size(dM));

%% Aggregation


%N = N(:);

for k = 1:length(bi)
    
    ii = bi(k)+1;
    jj = bj(k)+1;
    mi = m(ii);
    mj = m(jj);
    mij = m(ii)+m(jj);
    d00 = b300(k) + 1;
    d01 = b301(k) + 1;
    d10 = b310(k) + 1;
    d11 = b311(k) + 1;
    
    
    if length(beta) ==1;
        dN = alpha*beta*N(ii)*N(jj);
    else
        dN = alpha*beta(k)*N(ii)*N(jj);
    end
    

    if dN > 0
%         dM(ii) = dM(ii)-dN*m(ii);
%         dM(jj) = dM(jj)-dN*m(jj);
%         dM(d00) = dM(d00) + f00(k)*dN*(m(ii)+m(jj)); 
%         dM(d01) = dM(d01) + f01(k)*dN*(m(ii)+m(jj)); 
%         dM(d10) = dM(d10) + f10(k)*dN*(m(ii)+m(jj)); 
%         dM(d11) = dM(d11) + f11(k)*dN*(m(ii)+m(jj));
        dM(ii) = dM(ii)-dN*mi;
        dM(jj) = dM(jj)-dN*mj;
        dM(d00) = dM(d00) + f00(k)*dN*mij; 
        dM(d01) = dM(d01) + f01(k)*dN*mij; 
        dM(d10) = dM(d10) + f10(k)*dN*mij; 
        dM(d11) = dM(d11) + f11(k)*dN*mij;
    end
end

%% Degradation
if remin~=0
for i = 1:length(N)
    tmp = xz(i-1);
    xtmp = tmp(1);
    ztmp = tmp(2);
    if mod(i,nD)==0
        dN_remin = remin*q^xtmp *((1-q)/(3-a))*(-ztmp*N(i));
        
    else
        dN_remin = remin*q^xtmp *((1-q)/(3-a))*((ztmp+1)*N(i+1)-ztmp*N(i));
    end
    dMremin(i) = dN_remin*m(i);
end
end
%% Fragmentation
for i  = 1:length(M)
    if i< nD+1
        dMfrag(i) = pfrag(i+nD).*M(i+nD)*(frag_div);
    elseif i>length(M)-nD
        dMfrag(i) = -pfrag(i).*M(i)*(frag_div);
    else
         dMfrag(i) = pfrag(i+nD)*M(i+nD)*frag_div - pfrag(i)*M(i)*frag_div;
    end
end
 

%% collecting

if sum(dMfrag) >1E-5
    keyboard
end

dMdt = dM +prod(:) - M.*wWhites(:)./H +dMremin +dMfrag;


end

