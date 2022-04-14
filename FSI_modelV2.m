function [x_corr, lags, FSI_spikes] = FSI_modelV2(n, cort, pFSI_l, p0, laser)
%% FSI Model part 2
% n = Number of FSI neurons. 
% cort = number of cortical neurons
% pFSI_l = Rate at which FSIs make lateral connections.
% p0 = % Baseline probability of firing. 
% laser = either = 0 or 1, which corresponds to 'off' or 'on' respectively, this will determine whether or not the
% trial includes laser stimulation. 

p_cort = 0.005; %Rate at which cortical neurons fire.
pmax = 0.7; %max firing probability if all cortical neurons fire. 
cc = 0.03; % Correlation Coefficient for lateral connections
conv = n^2 - cort; %Calculates the number of connections the model has to make more connections for. 
conv = round(conv);
conv_n = randsample((1:cort), conv, true); %Randomly draws the convergent connections. 
connections = [1:cort, conv_n];

cort = round(cort);
matrix = zeros(n, cort);
c_mat = zeros(n, n);
% The following lines of code determine the matrix of connectivity. The
% first conditional is to overcome issues of with pulling the connectivity
% randomly. 
if cort == n
    for i = 1:n
        c_mat(i, :) = randperm(n, n);
    end
else
    for i = 1:n
        a = randsample(connections, n); % Randomly draws from the possible cortical neurons. 
        c_mat(i, :) = a;
        [~, w] = unique(c_mat(i,:), 'stable'); % Checks to ensure that the same cortical neurons do not connect to the same FSI. 
        its = setdiff( 1:numel(c_mat(i, :)), w );% Second portion of the checking for overlapping.
        
        if isempty(its) == 0 && i < n-2 % If "its" is not empty the code will go and randomly draw new cortical neurons to connnect to. 
            j = 1;
            while isempty(its) == 0 % This while loop ensures the new cortical neuron drawn is also not repetitive. 
                tt = find(connections == c_mat(i, its(j)));
                h = connections;
                h(tt) = [];
                x = randsample(h,1);
                if ismember(x, c_mat(i, :)) == 0
                    c_mat(i, its(j)) = x;
                    its(j) = [];
                end
            end
        else
            for pp = 1:length(its)
                oo = 1;
                th = ismember(c_mat(i, its(pp)), c_mat(oo, :));
                if th == 0
                    [~, ia] = setdiff(c_mat(i, :), c_mat(oo,:));
                    cx = c_mat(i, its(pp));
                    c_mat(i, its(pp)) = c_mat(oo, ia(1));
                    c_mat(oo, ia) = cx;
                else
                    oo = oo + 1;
                end
            end
        end
        for z = 1:n % This for loop goes through and removes the cortical neurons that were sampled already.
            hi = find(connections == c_mat(i, z));
            if isempty(hi) == 0
                connections(hi(1)) = [];
            end
        end
    end
end

for i = 1:n % Once the connectivity has been checked this then binarizes the connectivity map. 
    matrix(i, c_mat(i, :)) = ones([1, n]);
end
% Now have to make the probability of having an FSI connection, the below
% lines calculate the FSIxFSI matrix of lateral connections. The matrix is
% semetric as electrical coupling is bi-directional. 
lat_c = zeros(n, n);
latc = rand(n);
lat_c(latc < (pFSI_l)) = 1;
lat_c = tril(lat_c, -1) + tril(lat_c, -1)';



%% Building the probability matrix
% The following lines calculate the slope of how the probability of an FSI
% firing changes given the number of inputs an FSI receives from cortex. 
prob_f = zeros(n, 1);
for i = 1:n % Estimate as a line
    prob_f(i) = p0 + (i-1)*((pmax-p0)/n);  
end
% for i = 1:n %Estimate as an asymptote.
%     prob_f(i) = pmax.*(log(i+1)./log(n+1));
% end

prob_f;

%% Setting up cortical input 
time = 1000; %ms; may make this parameter a input parameter for a function.
% Since the cortical neurons are completely independent, we can make our
% life easier and calculate all cortical neuron activity before running
% through our model. This is done by the following lines:
c_in = zeros(time, cort);
c_in = poissrnd(p_cort,[time, cort]); %Randomly draws values from a normal distribution centered @0.5 
c_in(c_in >= (1-p_cort)) = 1; %Compares all probability values if they are larger than the threshold it is one
c_in(c_in < (1-p_cort)) = 0;% If lower than the threshold value, firing is 0. 

if laser == 1
    laser_on = [50, 150, 250, 350, 450, 550, 650, 750, 850, 950]; % Hard code when the "laser" is on.
    c_in(laser_on, :) = 1; %Sets cortical neuron values to 1 for the laser being on.
elseif laser == 0
    
end
% The following lines introduce the refractory period and ensure that no
% cortical neuron is firing <1ms after it had previous fired. 
[aa, bb] = find(c_in == 1);
for i = 1:length(aa)
    if aa(i)+1 >= time
        
    else
        c_in(aa(i)+1, bb(i)) = 0;
    end
end
%% Simulating the network with given input and time
% The following lines are the bread and butter of the model and go through
% the iterative process of calculating an FSIs firing rate through
% iterative sampling given the cortical inputs and the FSIs prior
% experience. 

fsi_sp = zeros(time, n);%Matrix storing the probability of an FSI firing.
FSI_spikes = zeros(time, n); %Matrix of ones and zeros to represent spikes.
t = 0;
% Experimentally derived values for tau:
Tau = [0.0027, 0.0024, 0.0016, 0.0019, 0.0023, 0.002, 0.002, 0.0015, 0.0025, 0.0028, 0.0025, 0.0015, 0.0023, 0.002, 0.0018, 0.0022, 0.0013, 0.0015,0.003, 0.0027, 0.002];
tau = mean(Tau);
while t < time
    if t == 0
%         Sets intial FSI firing probabilities at the start of a trial.
        fsi_s = poissrnd(1,[1, n]);
        fsi_s(fsi_s < 1 - p0) = 0;
        fsi_s(fsi_s >= 1 - p0) = 1;
        FSI_spikes(1, :) = fsi_s;
        t = t+1;
    else
        xx = matrix*c_in(t,:)'; %Calculates the input onto an FSI. 
        if t == 1
            for i = 1:n
                fsi_sp(t, i) = prob_f(xx(i)+1); %Work around to avoid errors if t-1 = 0. (Matlab hates 0).
            end
        else
            for i = 1:n
                p_lat = cc.*((fsi_sp(t-1,:)-p0) * lat_c(:, i)); % The minus p0, ensures that only deflections in the probability are shared. Otherwise the FSIs cause each other to fire. 
                if FSI_spikes(t-1, i) == 1
                    fsi_sp(t, i) = 0;
                elseif xx(i) == 0 && fsi_sp(t-1, i) <= p0
                    fsi_sp(t, i) = p0 + p_lat;
                elseif xx(i) == 0 && fsi_sp(t-1, i) > p0
                    fsi_sp(t, i) = fsi_sp(t-1, i) - (fsi_sp(t-1, i) * tau) + p_lat;
                elseif xx(i) > 0
                    fsi_sp(t, i) = prob_f(xx(i)) + (fsi_sp(t-1, i)-fsi_sp(t-1, i)*tau) + p_lat;
                end
                prob = 1 - rand(1);
                if prob <= fsi_sp(t, i)
                    FSI_spikes(t, i) = 1;
                else
                    FSI_spikes(t, i) = 0;
                end
            end
        end
        t = t + 1;
    end
end

lag = 20; % Determines the lag for cross correlation analyses. 
% This calculates the overall correlation among the neurons. 
[x_corr, lags] = xcorr(FSI_spikes, lag, 'coeff');
sh_spikes = zeros(time, n);

