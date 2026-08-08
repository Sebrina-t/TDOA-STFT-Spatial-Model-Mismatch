clear;clc;close all;

%% basic parameters
fs=192e3;
c=343;
T=0.1;

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


%% source positions
qs=[0.003 1];          % target position
qi=[0.803 0.6];        % interference position


%% target signal s[n]
f1=35e3;
f2=55e3;
Ts=20e-3;

Ns=round(Ts*fs);
ns=0:Ns-1;
ts=ns/fs;

mu=(f2-f1)/Ts;         % chirp rate

s0=cos(2*pi*(f1*ts+0.5*mu*ts.^2));

% target window
win_s=hamming(Ns)';
s0=s0.*win_s;

% target starts at 30 ms
start=round(30e-3*fs)+1;

s=zeros(1,N);
s(start:start+Ns-1)=s0;


%% interference i[n]
rng(1);

i=randn(1,N);

fc=[30e3 60e3]/(fs/2);

b=fir1(200,fc,'bandpass');

i=filter(b,1,i);


%% set SIR
SIR=0;                 % dB

idx=start:start+Ns-1;

Ps=mean(s(idx).^2);
Pi=mean(i(idx).^2);

alpha=sqrt(Ps/(Pi*10^(SIR/10)));

i=alpha*i;


%% source-microphone distance
rs=zeros(1,M);
ri=zeros(1,M);

for m=1:M
    rs(m)=norm(qs-p(m,:));
    ri(m)=norm(qi-p(m,:));
end


%% absolute propagation delay
Ts_m=rs/c;
Ti_m=ri/c;

Ls=Ts_m*fs;            % target delay, samples
Li=Ti_m*fs;            % interference delay, samples


%% relative delay to Mic 1
Ds=Ls-Ls(1);
Di=Li-Li(1);

disp('Target relative delay (samples):')
disp(Ds)

disp('Interference relative delay (samples):')
disp(Di)


%% generate four microphone signals
Smic=zeros(M,N);
Imic=zeros(M,N);
Xmic=zeros(M,N);

for m=1:M

    Smic(m,:)=frac_delay(s,Ls(m));

    Imic(m,:)=frac_delay(i,Li(m));

    Xmic(m,:)=Smic(m,:)+Imic(m,:);

end


%% check SIR
fprintf('\nTarget power = %.4f\n',Ps);
fprintf('Interference power = %.4f\n',mean(i(idx).^2));
fprintf('SIR = %.2f dB\n',10*log10(Ps/mean(i(idx).^2)));


%% plot geometry
figure;

plot(p(:,1),p(:,2),'o','LineWidth',1.5);
hold on;

plot(qs(1),qs(2),'^','LineWidth',1.5);
plot(qi(1),qi(2),'s','LineWidth',1.5);

for m=1:M
    text(p(m,1),p(m,2),['  Mic ' num2str(m)]);
end

grid on;
axis equal;

xlabel('x Position (m)');
ylabel('y Position (m)');

legend('Microphones','Target','Interference');

title('Four-Microphone Array Geometry');


%% target components
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Smic(m,:));

    grid on;

    ylabel('Amplitude');

    title(['Mic ' num2str(m)]);

    xlim([0 T*1e3]);

end

xlabel('Time (ms)');

sgtitle('Target Components s_m[n]');


%% interference components
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Imic(m,:));

    grid on;

    ylabel('Amplitude');

    title(['Mic ' num2str(m)]);

    xlim([0 T*1e3]);

end

xlabel('Time (ms)');

sgtitle('Interference Components i_m[n]');


%% final received signals
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Xmic(m,:));

    grid on;

    ylabel('Amplitude');

    title(['Mic ' num2str(m)]);

    xlim([0 T*1e3]);

end

xlabel('Time (ms)');

sgtitle('Received Microphone Signals x_m[n]');


%% Stage 3: STFT
Nw=256;
Noverlap=192;
Nfft=256;

win_stft=hann(Nw,'periodic');

for m=1:M

    [Xm,f_stft,t_stft]=spectrogram( ...
        Xmic(m,:), ...
        win_stft, ...
        Noverlap, ...
        Nfft, ...
        fs);

    if m==1

        K=size(Xm,1);
        L=size(Xm,2);

        Xstft=zeros(K,L,M);

    end

    Xstft(:,:,m)=Xm;

end


fprintf('\n========== STFT ==========\n');

fprintf('Frequency bins K = %d\n',K);

fprintf('Time frames L = %d\n',L);

fprintf('Frequency spacing = %.2f Hz\n',fs/Nfft);


%% STFT magnitude
XdB=20*log10( ...
    abs(Xstft)/max(abs(Xstft(:)))+1e-12);


figure;

for m=1:M

    subplot(M,1,m);

    imagesc( ...
        t_stft*1e3, ...
        f_stft/1e3, ...
        XdB(:,:,m));

    axis xy;

    xlim([0 T*1e3]);

    ylim([25 65]);

    caxis([-40 0]);

    ylabel('Frequency (kHz)');

    title(['Mic ' num2str(m)]);

    colorbar;

end

xlabel('Time (ms)');

sgtitle('Four-Channel STFT Magnitude');


%% Stage 4: STFT of target and interference
Sstft=zeros(K,L,M);
Istft=zeros(K,L,M);

for m=1:M

    [S_temp,~,~]=spectrogram( ...
        Smic(m,:), ...
        win_stft, ...
        Noverlap, ...
        Nfft, ...
        fs);

    [I_temp,~,~]=spectrogram( ...
        Imic(m,:), ...
        win_stft, ...
        Noverlap, ...
        Nfft, ...
        fs);

    Sstft(:,:,m)=S_temp;

    Istft(:,:,m)=I_temp;

end


%% Stage 4: select one time-frequency point
k=65;

t0=45e-3;

[~,l]=min(abs(t_stft-t0));


fprintf('\n========== Spatial Observation ==========\n');

fprintf('Frequency = %.2f kHz\n', ...
    f_stft(k)/1e3);

fprintf('Time = %.2f ms\n', ...
    t_stft(l)*1e3);


%% spatial observation vectors
skl=zeros(M,1);
ikl=zeros(M,1);
xkl=zeros(M,1);

for m=1:M

    skl(m)=Sstft(k,l,m);

    ikl(m)=Istft(k,l,m);

    xkl(m)=Xstft(k,l,m);

end


disp('Target spatial component:')
disp(skl)

disp('Interference spatial component:')
disp(ikl)

disp('Total observation:')
disp(xkl)


%% check STFT linearity
error=norm(xkl-skl-ikl);

fprintf('x = s + i error = %.6e\n',error);


%% show magnitude and phase
disp('Target magnitude:')
disp(abs(skl))

disp('Target phase (degree):')
disp(angle(skl)*180/pi)


disp('Interference magnitude:')
disp(abs(ikl))

disp('Interference phase (degree):')
disp(angle(ikl)*180/pi)


disp('Observation magnitude:')
disp(abs(xkl))

disp('Observation phase (degree):')
disp(angle(xkl)*180/pi)


