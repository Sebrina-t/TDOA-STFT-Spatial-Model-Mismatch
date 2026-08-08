clear;clc;close all;

%% basic parameters
fs=192e3;
c=343;

T=0.05;
N=round(T*fs);

n=0:N-1;
t=n/fs;

M=4;
d=0.002;


%% microphone positions
p=zeros(M,2);

for m=1:M
    p(m,:)=[(m-1)*d 0];
end


%% target position
qs=[0.803 0.6];


%% chirp parameters
f1=30e3;
f2=70e3;

Ts=4e-3;

mu=(f2-f1)/Ts;

t_start=20e-3;


%% source-microphone distance
rs=zeros(1,M);

for m=1:M
    rs(m)=norm(qs-p(m,:));
end


%% propagation delay
Ts_m=rs/c;

Ls=Ts_m*fs;

Ds=Ls-Ls(1);

tau=Ts_m-Ts_m(1);


fprintf('========== Propagation ==========\n');

fprintf('Relative delay (samples):\n');
disp(Ds)

fprintf('Relative delay (us):\n');
disp(tau*1e6)

fprintf('Maximum relative delay = %.4f samples\n', ...
    max(Ds)-min(Ds));

fprintf('Maximum relative delay = %.2f us\n', ...
    (max(tau)-min(tau))*1e6);


%% standard microphone signals
% all microphones sampled at the same absolute times t

Smic_std=zeros(M,N);

for m=1:M

    u=t-t_start-Ts_m(m);

    Smic_std(m,:)=source_chirp(u,f1,mu,Ts);

end


%% aligned sampling
% microphone m is sampled at t+tau(m)
% original continuous microphone waveform is unchanged

Smic_align=zeros(M,N);

for m=1:M

    t_sample=t+tau(m);

    u=t_sample-t_start-Ts_m(m);

    Smic_align(m,:)=source_chirp(u,f1,mu,Ts);

end


%% check aligned time-domain signals
align_error=zeros(1,M);

for m=1:M

    align_error(m)=norm( ...
        Smic_align(m,:)-Smic_align(1,:)) / ...
        (norm(Smic_align(1,:))+1e-12);

end


fprintf('\n========== Time Alignment Check ==========\n');

fprintf('Relative error to Mic 1:\n');
disp(align_error)


%% plot standard signals
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Smic_std(m,:));

    grid on;

    xlim([21 28]);

    ylabel('Amplitude');

    title(['Mic ' num2str(m)]);

end

xlabel('Absolute Time (ms)');

sgtitle('Standard Synchronous Sampling');


%% plot aligned signals
figure;

for m=1:M

    plot(t*1e3,Smic_align(m,:));
    hold on;

end

grid on;

xlim([21 28]);

xlabel('Mic 1 Reference Time (ms)');
ylabel('Amplitude');

legend('Mic 1','Mic 2','Mic 3','Mic 4');

title('TDOA-Aligned Sampling');


%% STFT parameters
Nw=64;
Noverlap=48;
Nfft=256;

win_stft=hann(Nw,'periodic');


%% standard STFT
for m=1:M

    [Xm,f_stft,t_stft]=spectrogram( ...
        Smic_std(m,:), ...
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


fprintf('\n========== Standard STFT ==========\n');

fprintf('Frequency bins = %d\n',K);
fprintf('Time frames = %d\n',L);

fprintf('Window length = %d samples\n',Nw);

fprintf('Window duration = %.3f ms\n', ...
    Nw/fs*1e3);

fprintf('Hop size = %d samples\n', ...
    Nw-Noverlap);

fprintf('Frequency spacing = %.2f Hz\n', ...
    fs/Nfft);


%% aligned STFT
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


%% restore propagation phase
Srestore=zeros(K,L,M);

for m=1:M

    phase_restore=exp( ...
        -1j*2*pi*f_stft*Ds(m)/fs);

    Srestore(:,:,m)= ...
        Salign(:,:,m).*phase_restore;

end


%% steering vectors
A=zeros(K,M);

for k=1:K

    for m=1:M

        A(k,m)=exp( ...
            -1j*2*pi*f_stft(k)*Ds(m)/fs);

    end

end


%% spatial model error
Estd=zeros(K,L);
Erestore=zeros(K,L);

for k=1:K

    a=A(k,:).';

    for l=1:L

        x_std=squeeze(Sstd(k,l,:));

        x_restore=squeeze(Srestore(k,l,:));


        % best scalar source coefficient
        S_est_std=(a'*x_std)/(a'*a);

        S_est_restore=(a'*x_restore)/(a'*a);


        % normalized spatial model error
        Estd(k,l)=norm( ...
            x_std-a*S_est_std) / ...
            (norm(x_std)+1e-12);


        Erestore(k,l)=norm( ...
            x_restore-a*S_est_restore) / ...
            (norm(x_restore)+1e-12);

    end

end


%% target-active TF points
mag_ref=abs(Sstd(:,:,1));

threshold=max(mag_ref(:))*10^(-30/20);

mask=mag_ref>threshold;


mean_std=mean(Estd(mask));
mean_restore=mean(Erestore(mask));

median_std=median(Estd(mask));
median_restore=median(Erestore(mask));


fprintf('\n========== Spatial Model Error ==========\n');

fprintf('Standard mean error = %.6e\n', ...
    mean_std);

fprintf('Proposed mean error = %.6e\n', ...
    mean_restore);

fprintf('Standard median error = %.6e\n', ...
    median_std);

fprintf('Proposed median error = %.6e\n', ...
    median_restore);


%% error maps
Estd_dB=20*log10(Estd+1e-12);

Erestore_dB=20*log10(Erestore+1e-12);

Estd_dB(~mask)=NaN;
Erestore_dB(~mask)=NaN;


figure;

subplot(2,1,1);

imagesc( ...
    t_stft*1e3, ...
    f_stft/1e3, ...
    Estd_dB);

axis xy;

ylim([25 75]);

caxis([-100 0]);

colorbar;

xlabel('Time (ms)');
ylabel('Frequency (kHz)');

title('Standard STFT Spatial Model Error (dB)');


subplot(2,1,2);

imagesc( ...
    t_stft*1e3, ...
    f_stft/1e3, ...
    Erestore_dB);

axis xy;

ylim([25 75]);

caxis([-100 0]);

colorbar;

xlabel('Time (ms)');
ylabel('Frequency (kHz)');

title('Aligned STFT + Phase Restoration Error (dB)');


%% compare one frequency bin
[~,k0]=min(abs(f_stft-50e3));

figure;

plot( ...
    t_stft*1e3, ...
    Estd(k0,:), ...
    'LineWidth',1.2);

hold on;

plot( ...
    t_stft*1e3, ...
    Erestore(k0,:), ...
    'LineWidth',1.2);

grid on;

xlim([21 28]);

xlabel('Time (ms)');
ylabel('Normalized Model Error');

legend( ...
    'Standard STFT', ...
    'Aligned + Phase Restoration');

title('Spatial Model Error at 50 kHz');


%% continuous-time source model
function y=source_chirp(u,f1,mu,Ts)

y=zeros(size(u));

idx=(u>=0)&(u<Ts);

u0=u(idx);

w=0.54-0.46*cos(2*pi*u0/Ts);

phase=2*pi*( ...
    f1*u0+0.5*mu*u0.^2);

y(idx)=cos(phase).*w;

end