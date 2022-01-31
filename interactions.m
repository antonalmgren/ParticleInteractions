function [dMdt] = interactions(M,m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = M./m(:);

dM = zeros(size(N));

dMremin  = zeros(size(dM));

dMfrag = zeros(size(dM));

%% Aggregation
% Keep inside state space
b300(b300>(L-1)) = L-1;
b301(b301>(L-1)) = L-1;
b310(b310>(L-1)) = b310(b310>(L-1))-nD;
b311(b311>(L-1)) = b311(b311>(L-1))-nD;
b311(b311>(L-1)) = b311(b311>(L-1)) -1;

N = N(:);


for k = 1:length(bi)
    %index of parent bins
    ii = bi(k)+1;
    jj = bj(k)+1;
    % index of daughter bins
    d00 = b300(k) + 1;
    d01 = b301(k) + 1;
    d10 = b310(k) + 1;
    d11 = b311(k) + 1;
    
    
    dN = alpha*beta(k)*N(ii)*N(jj);

    %if dN > 0
        dM(ii) = dM(ii)-dN.*m(ii);
        dM(jj) = dM(jj)-dN.*m(jj);
        dM(d00) = dM(d00) + f00(k)*dN.*(m(ii)+m(jj)); 
        dM(d01) = dM(d01) + f01(k)*dN.*(m(ii)+m(jj)); 
        dM(d10) = dM(d10) + f10(k)*dN.*(m(ii)+m(jj)); 
        dM(d11) = dM(d11) + f11(k)*dN.*(m(ii)+m(jj));
   % end
end
%% Degradation
for i = 1:length(N)
    tmp = xz(i-1);
    xtmp = tmp(1);
    ztmp = tmp(2);
    if mod(i,nD)==0
        dN_remin = remin*q^xtmp *((1-q)/(3-a))*(-ztmp*N(i));
        
    else
        dN_remin = remin*q^xtmp *((1-q)/(3-a))*((ztmp+1)*N(i+1)-ztmp*N(i));
    end
    dMremin(i) = dMremin(i) + dN_remin*m(i);
end

%% Fragmentation
dMfrag(nD+1:end) = -pfrag(nD+1:end).*M(nD+1:end)*frag_div;
dMfrag(1:end-nD) = pfrag(nD+1:end).*M(nD+1:end)*frag_div;

%% collecting

dMdt = dM  + prod(:) +dMremin - M.*wWhites(:)/H +dMfrag;

% %% Degradation
% 
% 
% for i = 1:length(N)
%     tmp = xz(i-1);
%     xtmp = tmp(1);
%     ztmp = tmp(2);
%     if mod(i,nD)==0
%         dN_remin = remin*q^xtmp *((1-q)/(3-a))*(-ztmp*N(i));
%         
%     else
%         dN_remin = remin*q^xtmp *((1-q)/(3-a))*((ztmp+1)*N(i+1)-ztmp*N(i));
%     end
%     dMremin(i) = dMremin(i) + dN_remin*m(i);
%     
%     %clear dN_remin
% end
%         
% 
% 
% dM = dM + dMremin- M.*wWhites(:)/H + prod(:) ;
% % dM = dM + - M.*wWhites(:)/H + prod(:) ;

end

