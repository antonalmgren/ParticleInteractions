function plotting(sim)
M = sim.M;
r = sim.r;
y = sim.y;
N = sim.N;
DELTA = sim.DELTA;
RMSE = sim.RMSE;
w = sim.w;
nD = sim.nD;

%calculate spectrum
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
Nspec = Nc./(DELTA*1E-4);

% export flux
export = M.*w/sim.H;
export = sum(export,1,'omitnan')/sim.prod_tot;

slope  = @(x,b)   1E-5*x.^(-b); 

figure('Position',[300,200,600,800])
tiledlayout(4,1,"TileSpacing","tight")

nexttile
loglog(r*2*1E-4,Nspec)
hold on
loglog(r*2*1E-4,slope(r*2*1E-4,3))
loglog(r*2*1E-4,slope(r*2*1E-4,4))
xlabel('diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
legend('Particle spectrum','slope=-3','slope=-4')
ylim([1E-5 1E10])

nexttile
surface(r,1:nD,M)
shading flat
title('M')
c = colorbar;
set(gca,'ColorScale','log','xscale','log')
c.Label.String = '\mu g C m^{-3}';

nexttile
surface(r,1:nD,N)
shading flat
title('N')
c = colorbar;
set(gca,'Colorscale','log','xscale','log')
c.Label.String = '# m^{-3}';

nexttile
semilogx(r,export)
xlabel('radius [\mu m]')
ylabel('export [% of total production]' )


%% Diagnostics

% for i = 1:length(t)
%     Mt = reshape(dM(i,:),nD,nR);
%     [Mdt,Mremin,Mfrag] = interactionsDT(t(i),Mt(:),m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div);
%     frag(i) = sum(Mfrag);
%     COM(i) = sum(Mdt)-sum(Mremin)+sum(Mt(:).*wWhites(:)/H)-sum(prod(:)) ;
% end

% figure
% plot(t,COM)

end


% 
% figure
% surface(x,z,N)
% title('N transient')
% colorbar
% set(gca,'ColorScale','log')
% set(gca,'ZScale','log')
% 
% for i = 1:length(t)
%     MM(:,:,i) = reshape(dM(i,:),nD,nR);
% end
% 
% exportDT = wWhites.*M/H;
% exportDT_x = sum(exportDT,1);
% exportFlux = sum(exportDT,'all','omitnan')*1E-3 %[mgC/m2/d]
% 
% figure
% plot(x,exportDT_x)
% title('export flux DT')
% text(0.01, 0.95, ['Export flux ~ ',num2str(round(exportFlux)),'mgC/m^2/d'] ,'Units','normalized')
% xlabel('size')
% ylabel('\mu g C m^{-2} d^{-1}') 
% 
% 
% 
% %% for annual retreat
% 
% Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
% 
% DELTA(1:nR-1) = 1E-4*(r(2:nR)-r(1:nR-1));
% DELTA(nR) = 1E-4*r(nR+1)-r(nR);
% NNN = Nc./DELTA;
% 
% 
% 
% slope  = @(x,b)   1E-5*x.^(-b); 
% 
% monterey = @(x) 0.027*x.^(-2.965);
% 
% 
% figure
% loglog(2*r(x)*1E-4,NNN, 'LineWidth',2)
% hold on
% plot(2*r(x)*1E-4,slope(2*r(x)*1E-4,4))
% plot(2*r(x)*1E-4,monterey(2*r(x)*1E-4))
% title('Particle size spectrum')
% xlabel('Particle Diameter [cm]')
% ylabel('Number spectrum [# cm^{-4}]')
% legend('Size spectrum','slope = -4','monterey bay','Location','SouthWest')
% set(gca,'FontSize',16)
% 
% 
% %%
% N3 = N;
% N3(N3==0) = NaN;
% N3 = log(N3);
% %N3 = log10(N3);
% N3 = N3+abs(min(N3,[],'all'));
% N3(N3==0) = NaN;
% 
% w_weight = sum(wWhites.*N,1)./sum(N,1);
% 
% figure
% loglog(r(x),w_weight,'-','LineWidth',2) 
% hold on
% scatter(r(xMesh(:)),wWhites(:),N3(:))
% xlabel('Radius (\mu m)')
% ylabel('Sinking velocity, (m/d)')
% title('Mean sinking velocity of aggregates')
% 
% 
% %Weight sinking speed by particle numbers
% 
% 

% 

% %% w comparison
% 
% B117 = 1.10; % [cm^-0.17 s^-1] (Kries & Evans 2000)
% 
% B062 = (1.10^(-1/0.17))^(0.38);
% 
% w_PL = @(r,eta,B) B * (2*r.*1E-4).^eta*3600*24*1E-2;
% 
% W_PL117 = w_PL(r(x),1.17,B117);
% W_PL062 = w_PL(r(x),0.62,B062);
% %%
% % 
% % figure
% % loglog(r(x),W_PL117,'r-',r(x),W_PL062,'g-',r(x),W,'b-',r(x),wWhites,'m-')
% % % hold on
% % % loglog(1E-3*r(1,:),W_PL117,'r-')
% % %xlim([r(1,1) r(1,end)])
% % %ylim([0 1000])
% % legend('\eta = 1.17','\eta = 0.62', 'Stokes','Whites')
% % ylabel('Sinking velocity [m d^{-1}]')
% % xlabel('radius [\mu m]')
% 
% figure
% loglog(r(x),W,'b-',r(x),wWhites,'r-')
%  hold on
%  loglog(r(x),W_PL117,'g-')
%  loglog(r(x),W_PL062,'m-')
% %xlim([1E2 r(x(end))])
% %ylim([0 1000])
% ylabel('Sinking velocity [m d^{-1}]')
% xlabel('radius [\mu m]')
% set(gca,'FontSize', 14)
% 
% %% Velocity figures
% % 
% % 
% % 
% % figure
% % surface(Re)
% % colorbar
% % %set(gca,'xscale','log')
% % set(gca,'ColorScale','log')
% % title('Re')
% % 
% % figure
% % subplot(1,3,1)
% % surface(W)
% % colorbar;
% % title('w')
% % %set(gca,'xscale','log')
% % set(gca,'ColorScale','log')
% % colormap jet
% % 
% % subplot(1,3,2)
% % surface(wWhites)
% % colorbar
% % title('w Whites approximation')
% % %set(gca,'xscale','log')
% % set(gca,'ColorScale','log')
% % colormap jet
% % 
% % subplot(1,3,3)
% % surface(W./wWhites)
% % colorbar
% % title('w:w_whites')
% % %set(gca,'xscale','log')
% %  %set(gca,'ColorScale','log')
% % colormap jet