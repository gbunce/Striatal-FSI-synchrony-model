%% FSI Run Script
% This file generates panels in fig 4 & 6 for n =  10. To generate other
% figures of various ends just go through each section and change n.
% To run this file you will need the following file:
% FSI_modelV2.m
n = 10;
pfsi = 0.3;
p0 = 0.005;
%% Laser on;
xrm = zeros(41, 20);
xrm2 = zeros(41, 20);
s = zeros(n, 20);
s2 = zeros(n, 20);
for i = 1:20
    [xr, lags, spikes] = FSI_modelV2(n, n, 0, p0, 1); %Fully Convergent
    [xr2, lags, spikes2] = FSI_modelV2(n, n^2, pfsi, p0, 1); %Fully independent
    xrm(:, i) = mean(xr, 2); % Calculates the average of the cross correlation, for the convergent.
    s(:, i) = sum(spikes, 1)';
    xrm2(:, i) = mean(xr2, 2);% Calculates the average of the cross correlation, for the Fully independent case.
    s2(:, i) =  sum(spikes2, 1)';
end

figure(1);
subplot(2, 1, 1)
plot(lags, mean(xrm, 2)); 
hold on; 
plot(lags, mean(xrm2, 2)); 
xlim([-20, 20]);
title('Laser on FSI Synchronny (Average)');
legend('Fully Convergent', 'Independent 30% Lat');
xlabel('time(ms)');
ylabel('c');

subplot(2, 1, 2)
plot(lags, xrm(:, 1)); 
hold on; 
plot(lags, xrm2(:, 1)); 
xlim([-20, 20]);
title('Laser on FSI Synchronny (One sample)');
legend('Fully Convergent', 'Independent 30% Lat');
xlabel('time(ms)');
ylabel('c');


figure(2)
histogram(s, 5);
hold on
histogram(s2, 10);
legend('independent', 'connected');

%% Laser off

xrm = zeros(41, 20);
xrm2 = zeros(41, 20);
for i = 1:20
    [xr, lags, spikes] = FSI_modelV2(n, n, 0, p0, 0); %Fully Convergent
    [xr2, lags, spikes2] = FSI_modelV2(n, n^2, pfsi, p0, 0); %Fully independent
      xrm(:, i) = mean(xr, 2); % Calculates the average of the cross correlation, for the convergent.
    s(:, i) = sum(spikes, 1)';
    xrm2(:, i) = mean(xr2, 2);% Calculates the average of the cross correlation, for the Fully independent case.
    s2(:, i) =  sum(spikes2, 1)';
end
figure(3);
subplot(2, 1, 1)
plot(lags, mean(xrm, 2));% Is this kosher.
hold on; 
plot(lags, mean(xrm2, 2)); 
xlim([-20, 20]);
title('Laser off FSI Synchronny (Average)');
legend('Fully Convergent', 'Independent 30% Lat');
xlabel('time(ms)');
ylabel('c');

subplot(2, 1, 2)
plot(lags, xrm(:, 1)); 
hold on; 
plot(lags, xrm2(:, 1)); 
xlim([-20, 20]);
title('Laser off FSI Synchronny (1 sample)');
legend('Fully Convergent', 'Independent 30% Lat');
xlabel('time(ms)');
ylabel('c');


figure(4)
histogram(s, 5);
hold on
histogram(s2, 10);
legend('independent', 'connected');

%% Testing changing the probability of lateral connections. 
pfsi = 0:0.2:1;
n = 10;
p0 = 0.005;
xrm = zeros(41, 20);
xrm2 = zeros(41, 20);
xrm3 = zeros(41, 20);
xrm4 = zeros(41, 20);
xrm5 = zeros(41, 20);
xrm6 = zeros(41, 20);
s = zeros(n, 20);
s2 = zeros(n, 20);
s3 = zeros(n, 20);
s4 = zeros(n, 20);
s5 = zeros(n, 20);
s6 = zeros(n, 20);

for i = 1:20
    [xr, lags, spikes] = FSI_modelV2(n, n^2, pfsi(1), p0, 0); %Fully Convergent
    [xr2, lags, spikes2] = FSI_modelV2(n, n^2, pfsi(2), p0, 0); 
    [xr3, lags, spikes3] = FSI_modelV2(n, n^2, pfsi(3), p0, 0);
    [xr4, lags, spikes4] = FSI_modelV2(n, n^2, pfsi(4), p0, 0);
    [xr5, lags, spikes5] = FSI_modelV2(n, n^2, pfsi(5), p0, 0);
    [xr6, lags, spikes6] = FSI_modelV2(n, n^2, pfsi(6), p0, 0);
    xrm(:, i) = mean(xr, 2); % Calculates the average of the cross correlation, for the convergent.
    xrm2(:, i) = mean(xr2, 2);% Calculates the average of the cross correlation, for the Fully independent case.
    xrm3(:, i) = mean(xr3, 2);
    xrm4(:, i) = mean(xr4, 2);
    xrm5(:, i) = mean(xr5, 2);
    xrm6(:, i) = mean(xr6, 2);
    s(:, i) = sum(spikes, 1)';
    s2(:, i) = sum(spikes2, 1)';
    s3(:, i) = sum(spikes3, 1)';
    s4(:, i) = sum(spikes4, 1)';
    s5(:, i) = sum(spikes5, 1)';
    s6(:, i) = sum(spikes6, 1)';
end

figure(5)
plot(lags, mean(xrm, 2)); 
hold on; 
plot(lags, mean(xrm2, 2)); 
plot(lags,mean(xrm3, 2));
plot(lags,mean(xrm4, 2));
plot(lags,mean(xrm5, 2));
plot(lags,mean(xrm6, 2));
title('Changing lateral connectivity (laser off)');
ylabel('c'); xlabel('time (ms)');
legend('% laterals = 0%', '% laterals = 20%', '% laterals = 40%', '% laterals = 60%', '% laterals = 80%', '% laterals = 100%')
ylim([-.01, 0.1])

figure(6)
histogram(s, 5);
hold on
histogram(s2, 10);
histogram(s3, 10);
histogram(s4, 10);
histogram(s5, 10);
histogram(s6, 10);
legend('% laterals = 0%', '% laterals = 20%', '% laterals = 40%', '% laterals = 60%', '% laterals = 80%', '% laterals = 100%')

%% Changing Convergence. 
pfsi = 0;
n = 10;
p0 = 0.005;
cort = [1, 1.5, 3, 5, 7.5, n].*n;
xrm = zeros(41, 20);
xrm2 = zeros(41, 20);
xrm3 = zeros(41, 20);
xrm4 = zeros(41, 20);
xrm5 = zeros(41, 20);
xrm6 = zeros(41, 20);
clear s; clear s2; clear s3; clear s4; clear s5; clear s6;

s = zeros(n, 20);
s2 = zeros(n, 20);
s3 = zeros(n, 20);
s4 = zeros(n, 20);
s5 = zeros(n, 20);
s6 = zeros(n, 20);

for i = 1:20
    [xr, lags, spikes] = FSI_modelV2(n, cort(1), pfsi, p0, 0); %Fully Convergent
    [xr2, lags, spikes2] = FSI_modelV2(n, cort(2), pfsi, p0, 0); %Fully independent
    [xr3, lags, spikes3] = FSI_modelV2(n, cort(3), pfsi, p0, 0);
    [xr4, lags, spikes4] = FSI_modelV2(n, cort(4), pfsi, p0, 0);
    [xr5, lags, spikes5] = FSI_modelV2(n, cort(5), pfsi, p0, 0);
    [xr6, lags, spikes6] = FSI_modelV2(n, cort(5), pfsi, p0, 0);
    xrm(:, i) = mean(xr, 2); % Calculates the average of the cross correlation, for the convergent.
    xrm2(:, i) = mean(xr2, 2);% Calculates the average of the cross correlation, for the Fully independent case.
    xrm3(:, i) = mean(xr3, 2);
    xrm4(:, i) = mean(xr4, 2);
    xrm5(:, i) = mean(xr5, 2);
    xrm6(:, i) = mean(xr6, 2);
    s(:, i) = sum(spikes, 1)';
    s2(:, i) = sum(spikes2, 1)';
    s3(:, i) = sum(spikes3, 1)';
    s4(:, i) = sum(spikes4, 1)';
    s5(:, i) = sum(spikes5, 1)';
    s6(:, i) = sum(spikes6, 1)';
end

figure(7)
plot(lags, mean(xrm, 2)); 
hold on; 
plot(lags, mean(xrm2, 2)); 
plot(lags,mean(xrm3, 2));
plot(lags,mean(xrm4, 2));
plot(lags,mean(xrm5, 2));
plot(lags,mean(xrm6, 2));
title('Changing percent Convergence (laser off)');
ylabel('c'); xlabel('time (ms)');
legend('100% Convergent', '85% convergent','70% convergent', '50% convergent', '25% convergent', '0% convergent');
ylim([-.01, 0.1])

figure(8)
histogram(s, 5);
hold on
histogram(s2, 10);
histogram(s3, 10);
histogram(s4, 10);
histogram(s5, 10);
histogram(s6, 10);
legend('100% Convergent', '85% convergent','70% convergent', '50% convergent', '25% convergent', '0% convergent');
c = {xrm; xrm2};
filename = sprintf('FSI_%d_data.xlsx', n);
writecell(c,filename, 'Sheet', 2, 'Range', 'A1:T41');