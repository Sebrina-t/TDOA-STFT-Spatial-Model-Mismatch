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
qs=[0.003 1];
qi=[0.803 0.6];


%% target signal s[n]
f1=35e3;
f2=55e3;
Ts=20e-3;

Ns=round(Ts*fs);
ns=0:Ns-1;
ts=ns/fs;

mu=(f2-f1)/Ts;

s0=cos(2*pi*(f1*ts+0.5*mu*ts.^2));

win=hamming(Ns)';
s0=s0.*win;

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
SIR=0;

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
Ls=rs/c*fs;
Li=ri/c*fs;


%% relative delay to Mic 1
Ds=Ls-Ls(1);
Di=Li-Li(1);

disp('Target relative delay (samples):')
disp(Ds)

disp('Interference relative delay (samples):')
disp(Di)


%% generate microphone signals
Smic=zeros(M,N);
Imic=zeros(M,N);
Xmic=zeros(M,N);

for m=1:M

    Smic(m,:)=frac_delay(s,Ls(m));

    Imic(m,:)=frac_delay(i,Li(m));

    Xmic(m,:)=Smic(m,:)+Imic(m,:);

end


%% check SIR
Pi2=mean(i(idx).^2);
SIR_check=10*log10(Ps/Pi2);

fprintf('\nTarget power = %.4f\n',Ps);
fprintf('Interference power = %.4f\n',Pi2);
fprintf('SIR = %.2f dB\n',SIR_check);


%% array geometry
figure;

plot(p(:,1),p(:,2),'o','LineWidth',1.5);
hold on;

plot(qs(1),qs(2),'^','LineWidth',1.5);
plot(qi(1),qi(2),'s','LineWidth',1.5);

grid on;
axis equal;

xlabel('x (m)');
ylabel('y (m)');
legend('Microphones','Target','Interference');

title('Array Geometry');


%% clean source signals
figure;

subplot(2,1,1);
plot(t*1e3,s);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Target Signal s[n]');

subplot(2,1,2);
plot(t*1e3,i);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Interference Signal i[n]');


%% target components
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Smic(m,:));

    ylabel(['Mic ' num2str(m)]);
    grid on;

end

xlabel('Time (ms)');
sgtitle('Target Components');


%% interference components
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Imic(m,:));

    ylabel(['Mic ' num2str(m)]);
    grid on;

end

xlabel('Time (ms)');
sgtitle('Interference Components');


%% received signals
figure;

for m=1:M

    subplot(M,1,m);

    plot(t*1e3,Xmic(m,:));

    ylabel(['Mic ' num2str(m)]);
    grid on;

end

xlabel('Time (ms)');
sgtitle('Received Signals');


%% zoom interference
figure;

index=round(40e-3*fs):round(40.1e-3*fs);

for m=1:M

    plot(t(index)*1e3,Imic(m,index));
    hold on;

end

grid on;

xlabel('Time (ms)');
ylabel('Amplitude');

legend('Mic 1','Mic 2','Mic 3','Mic 4');

title('Interference Signal - Zoomed');


%% fractional delay phase validation
f0=48e3;

test=cos(2*pi*f0*t);

Testmic=zeros(M,N);

for m=1:M
    Testmic(m,:)=frac_delay(test,Li(m));
end


%% measure phase
index=round(20e-3*fs):round(80e-3*fs);

z=zeros(1,M);

for m=1:M

    z(m)=sum(Testmic(m,index).* ...
        exp(-1j*2*pi*f0*t(index)));

end


phase_meas=zeros(1,M);

for m=1:M

    phase_meas(m)=angle(z(m)*conj(z(1)));

end

phase_meas=phase_meas*180/pi;


%% theoretical phase
phase_theory=-2*pi*f0*Di/fs;

phase_theory=angle(exp(1j*phase_theory));

phase_theory=phase_theory*180/pi;


%% phase error
phase_error=phase_meas-phase_theory;

phase_error=angle(exp(1j*phase_error*pi/180));

phase_error=phase_error*180/pi;


%% print phase result
fprintf('\n========== Phase Validation ==========\n');

disp('Theory phase difference (degree):')
disp(phase_theory)

disp('Measured phase difference (degree):')
disp(phase_meas)

disp('Phase error (degree):')
disp(phase_error)