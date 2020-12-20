clear; close all;

% data1=load('KGD_fullcouple_Inj1__fracData.txt');
data=load('iso_calibration_Reaction.txt');
% data3=load('KGD_fullcouple_Inj3__fracData.txt');
% data4=load('KGD_fullcouple_Inj4__fracData.txt');
% data5=load('KGD_fullcouple_Inj5__fracData.txt');

figure
plot(data(:,1), data(:,2),'sr','MarkerSize',10)
xlabel('Square root of injection rate, Q^{1/2}','FontSize',20)
ylabel('CMOD @ t = 3.5s (mm)','FontSize',20)
grid on
set(gca,'FontSize',15)
% print(gcf,'-depsc','demon_1.eps')


m=size(data2,1);
time = 0.01*(1:m);

% plot(time,data1(:,1),'-r','LineWidth',2)
% figure
% plot(time,data1(:,2),'-r','LineWidth',2)
% figure
% plot(time,data1(:,3),'-r','LineWidth',2)

% injectionRate = 1E-4.*[2 3 3.5 4 4.5 5 5.5];
% data = [1.142669e+06   1.688939e-04   1.120000e+00;
%     8.836000e+05   2.755945e-04   1.980000e+00;
%     8.473235e+05   3.092039e-04   2.280000e+00;
%     8.240078e+05   3.385167e-04   2.540000e+00;
%     8.102142e+05   3.640095e-04   2.760000e+00;
%     8.002025e+05   3.875604e-04   2.980000e+00;
%     7.966129e+05   4.072101e-04   3.140000e+00;];
 
data = [2.0  8.444270e+05   2.459905e-04   1.860000e+00;
2.5  7.481310e+05   3.083425e-04   2.520000e+00;
3.0  7.105216e+05   3.526504e-04   3.000000e+00;
3.5  6.899368e+05   3.894387e-04   3.400000e+00;
4.0  6.765879e+05   4.222097e-04   3.740000e+00;
4.5  6.677731e+05   4.519628e-04   4.060000e+00;
5.0  6.615755e+05   4.797466e-04   4.340000e+00;
5.5  6.574006e+05   5.061588e-04   4.600000e+00;];

injectionRate = data(:,1)*1E-4;
CMOD_3_5s = data(:,3);
CrackLength_3_5s= data(:,4);

figure
plot(injectionRate.^0.5, CMOD_3_5s.*1000,'sr','MarkerSize',10)
xlabel('Square root of injection rate, Q^{1/2}','FontSize',20)
ylabel('CMOD @ t = 3.5s (mm)','FontSize',20)
grid on
set(gca,'FontSize',15)
print(gcf,'-depsc','demon_1.eps')


figure

plot(injectionRate.^0.5, CrackLength_3_5s,'or','MarkerSize',10)
xlabel('Square root of injection rate, Q^{1/2}','FontSize',20)
ylabel('Fracture length @ t = 3.5s (m)','FontSize',20)
grid on
set(gca,'FontSize',15)
print(gcf,'-depsc','demon_2.eps')



figure;
plot(time,data2(:,1)./10^6,'-r','LineWidth',2)
hold on
plot(time,data3(1:1000,1)./10^6,'--k','LineWidth',2)
plot(time,data4(:,1)./10^6,':b','LineWidth',2)
plot(time,data5(:,1)./10^6,'-.c','LineWidth',2)
% plot(time,data5(:,1)./10^6,'-m','LineWidth',2)

legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','northeast')
xlabel('Injection time (s)','FontSize',20)
ylabel('CMP (MPa)','FontSize',20)
grid on
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_Presure.pdf')


figure;
plot(time(1:20:1000),2.*data2(1:20:1000,2).*10^3,'-*r','LineWidth',2,'MarkerSize',10)
hold on
plot(time(1:20:1000),2.*data3(1:20:1000,2).*10^3,'-+k','LineWidth',2,'MarkerSize',10)
plot(time(1:20:1000),2.*data4(1:20:1000,2).*10^3,'-sb','LineWidth',2,'MarkerSize',10)
plot(time(1:20:1000),2.*data5(1:20:1000,2).*10^3,'->c','LineWidth',2,'MarkerSize',10)
% plot(time(1:20:1000),2.*data5(1:20:1000,2).*10^3,'-om','LineWidth',2,'MarkerSize',10)

xlabel('Injection time (s)','FontSize',20)
ylabel('CMOD (mm)','FontSize',20)
legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','northwest')
grid on
% axis([0 10 0 5])
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_CMOD.pdf')


figure;
plot(time(1:20:1000),data2(1:20:1000,3),'-*r','LineWidth',2,'MarkerSize',10)
hold on
plot(time(1:20:1000),data3(1:20:1000,3),'-+k','LineWidth',2,'MarkerSize',10)
plot(time(1:20:1000),data4(1:20:1000,3),'-sb','LineWidth',2,'MarkerSize',10)
plot(time(1:20:1000),data5(1:20:1000,3),'->c','LineWidth',2,'MarkerSize',10)
% plot(time(1:20:1000),data5(1:20:1000,3),'-om','LineWidth',2,'MarkerSize',10)
xlabel('Injection time (s)','FontSize',20)
ylabel('Fracture Length (m)','FontSize',20)
% legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','northwest')
axis([0 10 0 9])
grid on
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_Crack_length.pdf')


% data1=load('KGD_fullcouple_Inj1__fracData.txt');
data2=load('KGD_fullcouple_Inj2__PP_Jump.txt');
data3=load('KGD_fullcouple_Inj3__PP_Jump.txt');
data4=load('KGD_fullcouple_Inj4__PP_Jump.txt');
data5=load('KGD_fullcouple_Inj5__PP_Jump.txt');



%pore presure distribution along fracture
figure;
plot(data2(1:10:size(data2,1),1)+10,data2(1:10:size(data2,1),2)./10^6,'-*r','LineWidth',2,'MarkerSize',10)
hold on
plot(data3(1:15:size(data3,1),1)+10,data3(1:15:size(data3,1),2)./10^6,'-+k','LineWidth',2,'MarkerSize',10)
plot(data4(1:20:size(data4,1),1)+10,data4(1:20:size(data4,1),2)./10^6,'-sb','LineWidth',2,'MarkerSize',10)
plot(data5(1:25:size(data5,1),1)+10,data5(1:25:size(data5,1),2)./10^6,'->c','LineWidth',2,'MarkerSize',10)
% plot(time(1:20:1000),data5(1:20:1000,3),'-om','LineWidth',2,'MarkerSize',10)
xlabel('Length (m)','FontSize',20)
ylabel('Water Poressure (MPa)','FontSize',20)
% legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','southwest')
axis([0 10 -0.4 0.6])
grid on
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_PorePressure.pdf')


%fracture opening 
figure;
plot(data2(1:10:size(data2,1),1)+10,data2(1:10:size(data2,1),3).*10^3,'-*r','LineWidth',2,'MarkerSize',10)
hold on
plot(data3(1:15:size(data3,1),1)+10,data3(1:15:size(data3,1),3).*10^3,'-+k','LineWidth',2,'MarkerSize',10)
plot(data4(1:20:size(data4,1),1)+10,data4(1:20:size(data4,1),3).*10^3,'-sb','LineWidth',2,'MarkerSize',10)
plot(data5(1:25:size(data5,1),1)+10,data5(1:25:size(data5,1),3).*10^3,'->c','LineWidth',2,'MarkerSize',10)
% plot(time(1:20:1000),data5(1:20:1000,3),'-om','LineWidth',2,'MarkerSize',10)
xlabel('Length (m)','FontSize',20)
ylabel('COD (mm)','FontSize',20)
% legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','southwest')
axis([0 10 0 0.8])
grid on
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_COD.pdf')



global PROP

PROP.GI = 90;                  %N/m
PROP.GII = PROP.GI;                  %N/m
PROP.sigmaMax = 1e6;              %N/m^2  Pa
PROP.tauMax = 1e6;              %N/m^2  Pa
PROP.lambdaN = 0.05;
PROP.lambdaT = 0.05;
PROP.alpha = 4;
PROP.beta = 4;
PROP.m = PROP.alpha*(PROP.alpha-1)*PROP.lambdaN^2/(1-PROP.alpha*PROP.lambdaN^2);
PROP.n = PROP.beta*(PROP.beta-1)*PROP.lambdaT^2/(1-PROP.beta*PROP.lambdaT^2);
PROP.deltaN = PROP.GI/PROP.sigmaMax*PROP.alpha*PROP.lambdaN*(1-PROP.lambdaN)^(PROP.alpha-1)*(PROP.alpha/PROP.m+1)*(PROP.lambdaN*PROP.alpha/PROP.m+1)^(PROP.m-1);
PROP.deltaT = PROP.GII/PROP.tauMax*PROP.alpha*PROP.lambdaT*(1-PROP.lambdaT)^(PROP.beta-1)*(PROP.beta/PROP.n+1)*(PROP.lambdaT*PROP.beta/PROP.n+1)^(PROP.n-1);
PROP.PenaltyStiffness = 1e8*PROP.sigmaMax/PROP.deltaN;
PROP.dGnt = 0;
PROP.dGtn = 0;
PROP.deltaN_conj = PROP.deltaN-PROP.deltaN*(PROP.dGnt/PROP.GI)^(1/PROP.alpha);
PROP.deltaT_conj = PROP.deltaT-PROP.deltaT*(PROP.dGtn/PROP.GI)^(1/PROP.beta);
PROP.GammaN = -PROP.GI*(PROP.alpha/PROP.m)^PROP.m;
PROP.GammaT = (PROP.beta/PROP.n)^PROP.n ;


traction2=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-data2(:,3)./PROP.deltaN).^PROP.alpha).*(PROP.m/PROP.alpha+data2(:,3)./PROP.deltaN).^(PROP.m-1)-...
                                PROP.alpha*(1-data2(:,3)./PROP.deltaN).^(PROP.alpha-1).*(PROP.m/PROP.alpha+data2(:,3)./PROP.deltaN).^PROP.m);
traction3=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-data3(:,3)./PROP.deltaN).^PROP.alpha).*(PROP.m/PROP.alpha+data3(:,3)./PROP.deltaN).^(PROP.m-1)-...
                                PROP.alpha*(1-data3(:,3)./PROP.deltaN).^(PROP.alpha-1).*(PROP.m/PROP.alpha+data3(:,3)./PROP.deltaN).^PROP.m);
traction4=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-data4(:,3)./PROP.deltaN).^PROP.alpha).*(PROP.m/PROP.alpha+data4(:,3)./PROP.deltaN).^(PROP.m-1)-...
                                PROP.alpha*(1-data4(:,3)./PROP.deltaN).^(PROP.alpha-1).*(PROP.m/PROP.alpha+data4(:,3)./PROP.deltaN).^PROP.m);
traction5=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-data5(:,3)./PROP.deltaN).^PROP.alpha).*(PROP.m/PROP.alpha+data5(:,3)./PROP.deltaN).^(PROP.m-1)-...
                                PROP.alpha*(1-data5(:,3)./PROP.deltaN).^(PROP.alpha-1).*(PROP.m/PROP.alpha+data5(:,3)./PROP.deltaN).^PROP.m);


ind = find(traction2<0);
traction2(ind) = 0;
ind = find(traction3<0);
traction3(ind) = 0;
ind = find(traction4<0);
traction4(ind) = 0;
ind = find(traction5<0);
traction5(ind) = 0;
                            
%pore presure distribution along fracture
figure;
plot(data2(1:6:size(data2,1),1)+10,traction2(1:6:size(data2,1),1)./10^6,'-*r','LineWidth',2,'MarkerSize',10)
hold on
plot(data3(1:6:size(data3,1),1)+10,traction3(1:6:size(data3,1),1)./10^6,'-+k','LineWidth',2,'MarkerSize',10)
plot(data4(1:6:size(data4,1),1)+10,traction4(1:6:size(data4,1),1)./10^6,'-sb','LineWidth',2,'MarkerSize',10)
plot(data5(1:6:size(data5,1),1)+10,traction5(1:6:size(data5,1),1)./10^6,'->c','LineWidth',2,'MarkerSize',10)
% plot(time(1:20:1000),data5(1:20:1000,3),'-om','LineWidth',2,'MarkerSize',10)
xlabel('Length (m)','FontSize',20)
ylabel('Cohesive attraction (MPa)','FontSize',20)
% legend({'Q=0.0002 m^2/s','Q=0.0003 m^2/s','Q=0.0004 m^2/s','Q=0.0005 m^2/s'},'FontSize',15,'Location','Northwest')
axis([0 10 0 1])
grid on
set(gca,'FontSize',15)
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); 
print(gcf,'-dpdf','Injection_cohesiveTraction.pdf')
