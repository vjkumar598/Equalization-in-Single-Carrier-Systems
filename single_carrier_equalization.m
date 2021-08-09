
% Equlization in Single Carrier Sytems
% To understand the need for equalization in digital communication systems,
% and to implement Linear MMSE equalizers adaptively.
%%
clear;clc;
close all;
%%
m = 4; % oversampling factor
l = 5; % time at which the transmit filter impulse response truncated on one side
alpha = 0.22; % roll-off or excess bandwith parameter
ts = 1/m; % Assuming the symbol duration T=1
%%
% random bits generate
nbits = 10100; % number of bits in simulation
a = randbit(nbits); % vector of bits
%%
% bits to symbol map
[b,A] = bpskmap(a); % BPSK Constellation Mapping
nsym = length(b); % number of transmitted symbols 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntraining = 100; npayload = 10000; % nsym = ntraining + npayload
% assuming the first 100 symbols are known and are used for training
b_training = b(1:ntraining);
b_payload = b(ntraining+1:nsym); % payload symbols are the actual symbols that carry information bits
t1 = 0:length(b)-1; % symbol times

%%
[gtx,t2] = root_raised_cosine(alpha,m,l); % setting transmit filter impulse response to be RRC(Root Raised Cosine) pulse
grx = fliplr(gtx); % setting receive filter impulse response to be Matched filter to transmit filter
 
%%
% sending upsampled symbols through tx filter 
% transmit filter output is a linearly modulated BPSK signal
b_upsampled = upsample(b,m);
t3 = 0:length(b_upsampled)-1;
%%
x = conv(b_upsampled,gtx); % BPSK signal
ts1 = t3(1) + t2(1);
t4= ts1:ts:ts1+ts*(length(x)-1);
%%
% specifying channel filter to be non ideal
gch = [-0.7,-0.3,0.3,0.5,1,0.9,0.8,-0.7,-0.8,0.7,0.8,0.6,0.3]; % Impulse response of the dispersive channel of interest
t5 = 0:ts:length(gch)-1;

%%
y1 = conv(x,gch); % noiseless input to receive filter
ts2 = t4(1) + t5(1);
t6= ts2:ts:ts2+ts*(length(y1)-1);
%%
z1 = ts*conv(y1,grx); % noiseless output of receive filter
ts3 = t6(1) + t2(1);
t7= ts3:ts:ts3+ts*(length(z1)-1);

%%
% cascade filters
q = conv(gtx,gch);
ts4 = t2(1)+t5(1);
t8 = ts4:ts:ts4+ts*(length(q)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = ts*conv(q,grx); % response of a single symbol b[n] is b[n]*p(t)
% Also we say this as the effective channel seen by a single symbol b[n]
ts5 = t8(1)+ t2(1);
t9 = ts5:ts:ts5+ts*(length(p)-1);

%%
%Code Fragment (Eye diagram) to visualize ISI introduced by dispersive
%channel

%remove edge effects before doing eye diagram
r1 = z1(100:499);
%horizontal display length in number of symbol intervals
K=2;
%break into non-overlapping traces
R1=reshape(r1,K*m,length(r1)/(K*m));
%now enforce continuity across traces
%(append to each trace the first element of the next trace)
row1 = R1(1,:);
L=length(row1);
row_pruned = row1(2:L);
R_pruned = R1(:,1:L-1);
R2 = [R_pruned;row_pruned];
time = (0:K*m)/m; %time as a multiple of symbol interval
figure;
plot(time,R2);
xlabel('t/T');

%%
% Adding noise
% # of samples in noise vector (i.e length of noise vector) = # of samples
% in y1
% Each sample in the noise vector is complex and IID Gussian distributed
% with mean 0 and some non-zero variance sigma_squared
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bn_energy = 1; % Average transmitted energy per symbol for BPSK with symbol alphabet A = [-1,1]  
Es = bn_energy*(q*q'); % Average received energy per symbol
Eb = Es/(log2(length(A))); % Average received bit energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For a given BER or Probability of bit error target. Using Ideal BER vs
% SNR plot we find the required Eb/N0 or SNR 
% Assuming the target BER is 10^-2 or 0.01 for which from Ideal BER cuve we
% find the required SNR which is aprox. to 4.3 dB
% for analysis
ebnodb = 4.3;
ebno = 10^(ebnodb/10);
N0 = Eb/ebno;
sigma = sqrt(N0/2);  % computing the noise standard deviation (or variance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise vector given as input  to receive filter 
    noise_real = sigma*randn(1,length(y1));
    noise = noise_real; % for real valued constellation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % for complex valued constellations
%     noise_imag = noise_real;
%     noise = noise_real + 1j*noise_imag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Assuming noise samples are at rate 1/ts
% Treating noise and receive filter (matched to transmit filter) as
% continuous functions
    z2 = ts * conv(noise,grx); % output of receive filter when noise as input
%%
% Recieve filter input and output equations
    y = y1 + noise;
    z = z1 + z2;
%%
% Adaptive LMMSE Equalizer
m1 = 2; % downsampling factor
r = z(1:m/m1:length(z)); % receive filter output downsampled
%%
h = p(1:m/m1:length(p)); % effective channel downsampled
L = length(h); % equalizer length
offset = 0;
%%
% initializing variables for adaptive implementation
phat = zeros(1,L);
Rhat = zeros(L,L);
for n=1:ntraining
    rn = r(1+m1*(n-1)+offset : L+m1*(n-1)+offset); % nth vector used to decide for symbol b[n]
    phat = phat + b(n)*rn;
    Rhat = Rhat + rn'*rn;
end
cLS = phat/Rhat; % solution to the linear equation xA=b form and the solution is x = b/A
%%
% implementing equalizer as a filter with impulse response g_equalizer 
% matched to cLs
g_equalizer = fliplr(cLS); % matched filter to cLS
zeq = conv(r,g_equalizer); % output of equalizer
%%
delay = length(g_equalizer) + offset; % start time value to sample the output of equalizer
zeq_samples = zeq(delay:m1:delay+(nsym-1)*m1); 
%%
b_payload_estimates = sign(zeq_samples(ntraining+1:nsym)); % Decision rule based on taking the sign of the 
% samples taken
errors = (b_payload_estimates ~= b_payload); % error vector
num_of_errors = sum(errors);
prob_of_error = num_of_errors/npayload; % BER of an equalized system
%%
% comparing with unequalized estimates
[maxval,maxval_loc] = max(h);
sampling_times = maxval_loc:m1:(nsym-1)*m1 + maxval_loc;
z_samples = r(sampling_times); % unequalized samples
%%
b_payload_estimates_unequalized = sign(z_samples(ntraining+1:nsym));
errors_unequalized = (b_payload_estimates_unequalized ~= b_payload);
num_of_errors_unequalized = sum(errors_unequalized ); 
prob_of_error_unequalized = num_of_errors_unequalized/npayload; % BER of an unequalized system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On running the code we see that the BER of an unequalized system , based
% on the decisions based on single samples is significantly worser than an
% equalized system, where more than one sample (here m1=2 samples )are
% taken to decide on a single symbol. 
% Also we see that there is a significant gap between the simulated BER of an
% equalized system with  that of the IDEAL BER computed theoretically at a
% given SNR in dB.






