clear;clc;close all;

%% load real waveform
[x,fs]=audioread('rcap00004_6.wav');

x=x(:,1)';
x=x-mean(x);

fprintf('========== Real Recording ==========\n');
fprintf('Sampling frequency = %.1f kHz\n',fs/1e3);
fprintf('Duration = %.3f ms\n',length(x)/fs*1e3);


%% extract bat event
% main event is around 69-107 ms
t1=65e-3;
t2=112e-3;

idx1=round(t1*fs)+1;
idx2=round(t2*fs);

s=x(idx1:idx2);

s=s/(max(abs(s))+1e-12);

N=length(s);

n=0:N-1;
t=n/fs;


fprintf('Selected segment = %.2f - %.2f ms\n', ...
    t1*1e3,t2*1e3);

fprintf('Segment duration = %.2f ms\n', ...
    N/fs*1e3);


%% plot real waveform
figure;

plot(t*1e3,s);

grid on;

xlabel('Time (ms)');
ylabel('Normalized Amplitude');

title('Real Bat Waveform');


%% source STFT
Nw_show=512;
Noverlap_show=384;
Nfft_show=2048;

win_show=hann(Nw_show,'periodic');

[Sshow,f_show,t_show]=spectrogram( ...
    s, ...
    win_show, ...
    Noverlap_show, ...
    Nfft_show, ...
    fs);

SdB=20*log10( ...
    abs(Sshow)/max(abs(Sshow(:)))+1e-12);


figure;

imagesc( ...
    t_show*1e3, ...
    f_show/1e3, ...
    SdB);

axis xy;

ylim([20 130]);
caxis([-70 0]);

colorbar;

xlabel('Time (ms)');
ylabel('Frequency (kHz)');

title('STFT of Real Bat Waveform');


%% four-microphone model
M=4;

% use approximately the same physical maximum TDOA
% as the previous simulation: 14 us
tau_max=14e-6;

Dmax=tau_max*fs;

Ds=linspace(0,-Dmax,M);

fprintf('\n========== Array Model ==========\n');

fprintf('Maximum TDOA = %.2f us\n',tau_max*1e6);

fprintf('Maximum TDOA = %.4f samples\n',Dmax);

fprintf('Relative delays (samples):\n');
disp(Ds)


%% absolute propagation delays
% common delay added only to keep every delay positive
Lref=20;

Lm=Lref+Ds;

fprintf('Absolute delays (samples):\n');
disp(Lm)


%% generate standard microphone signals
Smic=zeros(M,N);

for m=1:M

    Smic(m,:)=frac_delay( ...
        s,Lm(m));

end


%% conventional time alignment
% Mic 1 is reference
% Ds <= 0, therefore -Ds is a positive compensation delay
Smic_align=zeros(M,N);

for m=1:M

    Smic_align(m,:)=frac_delay( ...
        Smic(m,:), ...
        -Ds(m));

end


%% STFT parameters
% 256 samples at 750 kHz = 0.341 ms
Nw=256;
Noverlap=192;
Nfft=1024;

win_stft=hann(Nw,'periodic');


fprintf('\n========== STFT ==========\n');

fprintf('Window length = %d samples\n',Nw);

fprintf('Window duration = %.3f ms\n', ...
    Nw/fs*1e3);

fprintf('Hop size = %d samples\n', ...
    Nw-Noverlap);

fprintf('Frequency spacing = %.2f Hz\n', ...
    fs/Nfft);


%% standard STFT
for m=1:M

    [Xm,f_stft,t_stft]=spectrogram( ...
        Smic(m,:), ...
        win_stft, ...
        Noverlap, ...
        Nfft, ...
        fs);

    if m==1

        K=size(Xm,1);
        L=size(Xm,2);

        Sstd=zeros(K,L,M);

    end

    Sstd(:,:,m)=Xm;

end


%% time-aligned STFT
for m=1:M

    [Xm,~,~]=spectrogram( ...
        Smic_align(m,:), ...
        win_stft, ...
        Noverlap, ...
        Nfft, ...
        fs);

    if m==1

        Salign=zeros(K,L,M);

    end

    Salign(:,:,m)=Xm;

end


%% steering vector
A=zeros(K,M);

for k=1:K

    for m=1:M

        A(k,m)=exp( ...
            -1j*2*pi*f_stft(k)*Ds(m)/fs);

    end

end


%% phase-only compensation
Sphase=zeros(K,L,M);

for m=1:M

    phase_comp=exp( ...
        1j*2*pi*f_stft*Ds(m)/fs);

    Sphase(:,:,m)= ...
        Sstd(:,:,m).*phase_comp;

end


%% proposed: time alignment + phase restoration
Sprop=zeros(K,L,M);

for m=1:M

    phase_restore=exp( ...
        -1j*2*pi*f_stft*Ds(m)/fs);

    Sprop(:,:,m)= ...
        Salign(:,:,m).*phase_restore;

end


%% active TF mask
mag_ref=abs(Sstd(:,:,1));

% only analyse bat-frequency region
freq_mask=(f_stft>=20e3) & ...
          (f_stft<=130e3);

threshold=max( ...
    mag_ref(freq_mask,:), ...
    [], ...
    'all')*10^(-30/20);

mask=zeros(K,L);

for k=1:K

    if freq_mask(k)==1

        for l=1:L

            if mag_ref(k,l)>threshold

                mask(k,l)=1;

            end

        end

    end

end


fprintf('Active TF cells = %d\n',sum(mask(:)));


%% model error
res_std=0;
res_phase=0;
res_align=0;
res_prop=0;

power_std=0;
power_phase=0;
power_align=0;
power_prop=0;


for k=1:K

    a=A(k,:).';

    a_one=ones(M,1);


    for l=1:L

        if mask(k,l)==1

            x_std=squeeze( ...
                Sstd(k,l,:));

            x_phase=squeeze( ...
                Sphase(k,l,:));

            x_align=squeeze( ...
                Salign(k,l,:));

            x_prop=squeeze( ...
                Sprop(k,l,:));


            %% standard STFT
            S_est=(a'*x_std)/(a'*a);

            e=x_std-a*S_est;

            res_std=res_std+norm(e)^2;

            power_std= ...
                power_std+norm(x_std)^2;


            %% phase-only
            S_est=(a_one'*x_phase)/ ...
                  (a_one'*a_one);

            e=x_phase-a_one*S_est;

            res_phase=res_phase+norm(e)^2;

            power_phase= ...
                power_phase+norm(x_phase)^2;


            %% time alignment
            S_est=(a_one'*x_align)/ ...
                  (a_one'*a_one);

            e=x_align-a_one*S_est;

            res_align=res_align+norm(e)^2;

            power_align= ...
                power_align+norm(x_align)^2;


            %% proposed
            S_est=(a'*x_prop)/(a'*a);

            e=x_prop-a*S_est;

            res_prop=res_prop+norm(e)^2;

            power_prop= ...
                power_prop+norm(x_prop)^2;

        end

    end

end


%% final errors
Err_std=sqrt( ...
    res_std/(power_std+1e-30));

Err_phase=sqrt( ...
    res_phase/(power_phase+1e-30));

Err_align=sqrt( ...
    res_align/(power_align+1e-30));

Err_prop=sqrt( ...
    res_prop/(power_prop+1e-30));


fprintf('\n========== Real Waveform Validation ==========\n');

fprintf('Standard STFT             = %.6e\n', ...
    Err_std);

fprintf('Phase-only compensation   = %.6e\n', ...
    Err_phase);

fprintf('Time alignment            = %.6e\n', ...
    Err_align);

fprintf('Proposed aligned + phase  = %.6e\n', ...
    Err_prop);


fprintf('\n========== Consistency Check ==========\n');

fprintf('|Standard - Phase-only| = %.6e\n', ...
    abs(Err_std-Err_phase));

fprintf('|Time-align - Proposed| = %.6e\n', ...
    abs(Err_align-Err_prop));

fprintf('Standard / Time-align = %.2f\n', ...
    Err_std/(Err_align+1e-30));


%% comparison figure
Err=[ ...
    Err_std ...
    Err_phase ...
    Err_align ...
    Err_prop];

figure;

bar(Err);

set(gca,'YScale','log');

grid on;

xticklabels({ ...
    'Standard', ...
    'Phase-only', ...
    'Time-align', ...
    'Proposed'});

ylabel('Normalized Spatial Model Error');
ylim([1e-13 1e-1]);
title('Real Bat Waveform Validation');


%% fractional delay
function y=frac_delay(x,D)

N=length(x);

Dint=floor(D);
Dfrac=D-Dint;

L=81;
mid=(L-1)/2;

k=0:L-1;

h=sinc(k-mid-Dfrac);

h=h.*hamming(L)';

h=h/sum(h);

z=conv(x,h);

z=z(mid+1:mid+N);

y=zeros(1,N);

if Dint<N

    y(Dint+1:end)= ...
        z(1:N-Dint);

end

end