%% post plot
clear; close all; clc

%% pf v.s. n
% ss-da
ex1_2 = load('Ex1_2.mat');
ex1_4 = load('Ex1_4.mat');
ex1_10 = load('Ex1_10.mat');
ex1_50 = load('Ex1_50.mat');
ex1_100 = load('Ex1_100.mat');
ex1_200 = load('Ex1_200.mat');

Pf_da = [mean(ex1_2.Pf_da) mean(ex1_4.Pf_da) mean(ex1_10.Pf_da)...
    mean(ex1_50.Pf_da) mean(ex1_100.Pf_da) mean(ex1_200.Pf_da)];
cov_da = [nanmean(ex1_2.cov_da) nanmean(ex1_4.cov_da) nanmean(ex1_10.cov_da)...
    nanmean(ex1_50.cov_da) nanmean(ex1_100.cov_da) nanmean(ex1_200.cov_da)];
nfun_da = [mean(ex1_2.nfun_da) mean(ex1_4.nfun_da) mean(ex1_10.nfun_da)...
    mean(ex1_50.nfun_da) mean(ex1_100.nfun_da) mean(ex1_200.nfun_da)]-500;

% ss
ex1s_2 = load('Ex1s_2.mat');
ex1s_4 = load('Ex1s_4.mat');
ex1s_10 = load('Ex1s_10.mat');
ex1s_50 = load('Ex1s_50.mat');
ex1s_100 = load('Ex1s_100.mat');
ex1s_200 = load('Ex1s_200.mat');

Pf_ss = [mean(ex1s_2.Pf_ss) mean(ex1s_4.Pf_ss) mean(ex1s_10.Pf_ss)...
    mean(ex1s_50.Pf_ss) mean(ex1s_100.Pf_ss) mean(ex1s_200.Pf_ss)];
cov_ss = [nanmean(ex1s_2.cov_ss) nanmean(ex1s_4.cov_ss) nanmean(ex1s_10.cov_ss)...
    nanmean(ex1s_50.cov_ss) nanmean(ex1s_100.cov_ss) nanmean(ex1s_200.cov_ss)];
nfun_ss = [mean(ex1s_2.nfun_ss) mean(ex1s_4.nfun_ss) mean(ex1s_10.nfun_ss)...
    mean(ex1s_50.nfun_ss) mean(ex1s_100.nfun_ss) mean(ex1s_200.nfun_ss)];

(nfun_ss-nfun_da)./(nfun_ss-1000)

x1 = [2 4 10 50 100 200];
Pf_da_r = (Pf_da-1e-3)*1e3;
Pf_ss_r = (Pf_ss-1e-3)*1e3;
figure;
semilogx(x1,Pf_da_r,'k-o')
hold on
semilogx(x1,Pf_ss_r,'b-*')
legend('ss-da','ss')
axis([1 500 -0.06 0.08])
grid on
xlabel('Dimension n','FontName','Times New Roman')
ylabel('Relative bias','FontName','Times New Roman')

figure;
semilogx(x1,nfun_da,'k-o')
hold on
semilogx(x1,nfun_ss,'b-*')
legend('ss-da','ss')
axis([1 500 1000 5000])
grid on
xlabel('Dimension n','FontName','Times New Roman')
ylabel('N_{call}','FontName','Times New Roman')

figure;
semilogx(x1,cov_da,'k-o')
hold on
semilogx(x1,cov_ss,'b-*')
legend('ss-da','ss')
axis([1 500 0.2 0.28])
grid on
xlabel('Dimension n','FontName','Times New Roman')
ylabel('Empirical c.o.v.','FontName','Times New Roman')


%% Pf v.s.  pf
% ss-da
ex1_3 = load('Ex1_100.mat');
ex1_4 = load('Ex1_100_4.mat');
ex1_5 = load('Ex1_100_5.mat');
ex1_6 = load('Ex1_100_6.mat');
ex1_7 = load('Ex1_100_7.mat');
ex1_8 = load('Ex1_100_8.mat');

Pf_da = [mean(ex1_3.Pf_da) mean(ex1_4.Pf_da) mean(ex1_5.Pf_da)...
    mean(ex1_6.Pf_da) mean(ex1_7.Pf_da) mean(ex1_8.Pf_da)];
cov_da = [nanmean(ex1_3.cov_da) nanmean(ex1_4.cov_da) nanmean(ex1_5.cov_da)...
    nanmean(ex1_6.cov_da) nanmean(ex1_7.cov_da) nanmean(ex1_8.cov_da)];
nfun_da = [mean(ex1_3.nfun_da) mean(ex1_4.nfun_da) mean(ex1_5.nfun_da)...
    mean(ex1_6.nfun_da) mean(ex1_7.nfun_da) mean(ex1_8.nfun_da)]-500;

% ss
ex1s_3 = load('Ex1s_100.mat');
ex1s_4 = load('Ex1s_100_4.mat');
ex1s_5 = load('Ex1s_100_5.mat');
ex1s_6 = load('Ex1s_100_6.mat');
ex1s_7 = load('Ex1s_100_7.mat');
ex1s_8 = load('Ex1s_100_8.mat');

Pf_ss = [mean(ex1s_3.Pf_ss) mean(ex1s_4.Pf_ss) mean(ex1s_5.Pf_ss)...
    mean(ex1s_6.Pf_ss) mean(ex1s_7.Pf_ss) mean(ex1s_8.Pf_ss)];
cov_ss = [nanmean(ex1s_3.cov_ss) nanmean(ex1s_4.cov_ss) nanmean(ex1s_5.cov_ss)...
    nanmean(ex1s_6.cov_ss) nanmean(ex1s_7.cov_ss) nanmean(ex1s_8.cov_ss)];
nfun_ss = [mean(ex1s_3.nfun_ss) mean(ex1s_4.nfun_ss) mean(ex1s_5.nfun_ss)...
    mean(ex1s_6.nfun_ss) mean(ex1s_7.nfun_ss) mean(ex1s_8.nfun_ss)];

x2 = [1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
Pf_da_r = (Pf_da-x2)./x2;
Pf_ss_r = (Pf_ss-x2)./x2;

figure;
semilogx(x2,Pf_da_r,'k-o')
hold on
semilogx(x2,Pf_ss_r,'b-*')
legend('ss-da','ss')
axis([1e-9 1e-2 0 0.15])
grid on
xlabel('Magnitude of failure probability','FontName','Times New Roman')
ylabel('Relative bias','FontName','Times New Roman')

figure;
semilogx(x2,nfun_da,'k-o')
hold on
semilogx(x2,nfun_ss,'b-*')
legend('ss-da','ss')
axis([1e-9 1e-2 1000 10000])
grid on
xlabel('Magnitude of failure probability','FontName','Times New Roman')
ylabel('N_{call}','FontName','Times New Roman')

figure;
semilogx(x2,cov_da,'k-o')
hold on
semilogx(x2,cov_ss,'b-*')
legend('ss-da','ss')
axis([1e-9 1e-2 0.2 0.6])
grid on
xlabel('Magnitude of failure probability','FontName','Times New Roman')
ylabel('Empirical c.o.v.','FontName','Times New Roman')


%% Ex 2

ex2_1 = load('Ex2_1.mat');
ex2_2 = load('Ex2_2.mat');
ex2_3 = load('Ex2_3.mat');

ex2s_1 = load('Ex2s_1.mat');
ex2s_2 = load('Ex2s_2.mat');
ex2s_3 = load('Ex2s_3.mat');

Pf = [mean(ex2s_1.Pf_ss) mean(ex2s_2.Pf_ss) mean(ex2s_3.Pf_ss)...
    mean(ex2_1.Pf_da) mean(ex2_2.Pf_da) mean(ex2_3.Pf_da)];
cov = [nanmean(ex2s_1.cov_ss) nanmean(ex2s_2.cov_ss) nanmean(ex2s_3.cov_ss)...
    nanmean(ex2_1.cov_da) nanmean(ex2_2.cov_da) nanmean(ex2_3.cov_da)];
nfun = [mean(ex2s_1.nfun_ss) mean(ex2s_2.nfun_ss) mean(ex2s_3.nfun_ss)...
    mean(ex2_1.nfun_da)-300 mean(ex2_2.nfun_da)-300 mean(ex2_3.nfun_da)-300];




for i = 1:100;
    y1_da(i) = y_da{i}(4);
    y1_ss(i) = y_ss{i}(4);
    p1_da(i) = p_da{i}(4);
    p1_ss(i) = p_ss{i}(4);
end

figure;
plot(y1_da);
hold on
plot(y1_ss,'r--')

mean(y1_da)
mean(y1_ss)
mean(p1_da)
mean(p1_ss)

std(Pf_da)/mean(Pf_da)*mean(nfun_da-300)
std(Pf_ss)/mean(Pf_ss)*mean(nfun_ss)

mean(cov_da)
cov_ss(isnan(cov_ss))=[];
mean(cov_ss)

save('Ex2_2.mat','Pf_da','nfun_da','cov_da','accep','y_da','p_da')

save('Ex2s_2.mat','Pf_ss','nfun_ss','cov_ss','y_ss','p_ss')
