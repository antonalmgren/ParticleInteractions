function  SAplotting(simA,simAl,simE)

for i = 1:length(simA.sim)
    N = simA.sim(i).N;
    Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
    NspecA(i,:) = Nc./(simA.sim(i).DELTA*1E-4);
    exportA = (simA.sim(i).w.*simA.sim(i).M/simA.sim(i).H);
    exportAx(i,:) = sum(exportA,1,'omitnan');
end

tot_export = sum(exportAx,2)

r = simA.sim(1).r;
figure
subplot(1,2,1)
loglog(r*2E-4,NspecA')
ylabel('number spectrum [# cm^{-4}')
xlabel('diameter [cm]')

subplot(1,2,2)
semilogx(r*2E-4,exportAx')
ylabel('Export [%]')
xlabel('diameter [cm]')
end