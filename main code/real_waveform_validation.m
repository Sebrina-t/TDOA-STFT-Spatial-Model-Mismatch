clear;clc;close all;

%% =========================
% Real waveform files
% =========================

files={ ...
    'rcap00004_6.wav', ...
    'taegyp00004_10n.wav', ...
    'pcap00001_4.wav'};

% event start time
t1_list=[ ...
    65e-3, ...
    0.3e-3, ...
    7.5e-3];

% event end time
t2_list=[ ...
    112e-3, ...
    3.8e-3, ...
    14.5e-3];

Nfile=length(files);

% columns:
% Standard
% Phase-only
% Time-align
% Proposed
Err_all=zeros(Nfile,4);


%% =========================
% Validation loop
% =========================

for q=1:Nfile

    fprintf('\n\n');
    fprintf('============================================\n');
    fprintf('File %d: %s\n',q,files{q});
    fprintf('============================================\n');


    %% load waveform
    [x,fs]=audioread(files{q});

    x=x(:,1)';
    x=x-mean(x);

    fprintf('Sampling frequency = %.1f kHz\n',fs/1e3);
    fprintf('Original duration = %.3f ms\n', ...
        length(x)/fs*1e3);


    %% extract event
    t1=t1_list(q);
    t2=t2_list(q);

    idx1=round(t1*fs)+1;
    idx2=round(t2*fs);

    idx1=max(idx1,1);
    idx2=min(idx2,length(x));

    s=x(idx1:idx2);

    s=s/(max(abs(s))+1e-12);

    N=length(s);

    fprintf('Selected interval = %.3f - %.3f ms\n', ...
        t1*1e3,t2*1e3);

    fprintf('Selected duration = %.3f ms\n', ...
        N/fs*1e3);


    %% =========================
    % Four-microphone propagation
    % =========================

    M=4;

    % same physical maximum TDOA for all recordings
    tau_max=14e-6;

    Dmax=tau_max*fs;

    % relative TDOA in samples
    Ds=linspace(0,-Dmax,M);

    fprintf('Maximum TDOA = %.2f us\n', ...
        tau_max*1e6);

    fprintf('Maximum TDOA = %.4f samples\n', ...
        Dmax);


    %% common propagation delay
    % make all absolute delays positive
    Lref=20;

    Lm=Lref+Ds;


    %% generate microphone signals
    Smic=zeros(M,N);

    for m=1:M

        Smic(m,:)=frac_delay( ...
            s,Lm(m));

    end


    %% =========================
    % Conventional time alignment
    % =========================

    Smic_align=zeros(M,N);

    for m=1:M

        Smic_align(m,:)=frac_delay( ...
            Smic(m,:), ...
            -Ds(m));

    end


    %% =========================
    % STFT
    % =========================

    Nw=256;
    Noverlap=192;
    Nfft=1024;

    win_stft=hann(Nw,'periodic');


    %% Standard STFT
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


    %% Time-aligned STFT
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


    %% =========================
    % Theoretical steering vector
    % =========================

    A=zeros(K,M);

    for k=1:K

        for m=1:M

            A(k,m)=exp( ...
                -1j*2*pi*f_stft(k)*Ds(m)/fs);

        end

    end


    %% =========================
    % Phase-only compensation
    % =========================

    Sphase=zeros(K,L,M);

    for m=1:M

        phase_comp=exp( ...
            1j*2*pi*f_stft*Ds(m)/fs);

        Sphase(:,:,m)= ...
            Sstd(:,:,m).*phase_comp;

    end


    %% =========================
    % Proposed representation
    % time alignment + phase restoration
    % =========================

    Sprop=zeros(K,L,M);

    for m=1:M

        phase_restore=exp( ...
            -1j*2*pi*f_stft*Ds(m)/fs);

        Sprop(:,:,m)= ...
            Salign(:,:,m).*phase_restore;

    end


    %% =========================
    % Active TF mask
    % =========================

    mag_ref=abs(Sstd(:,:,1));

    % ultrasonic / bat frequency range
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

    fprintf('Active TF cells = %d\n', ...
        sum(mask(:)));


    %% =========================
    % Spatial-model error
    % =========================

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


                %% Standard
                S_est=(a'*x_std)/(a'*a);

                e=x_std-a*S_est;

                res_std=res_std+norm(e)^2;

                power_std= ...
                    power_std+norm(x_std)^2;


                %% Phase-only
                S_est=(a_one'*x_phase)/ ...
                      (a_one'*a_one);

                e=x_phase-a_one*S_est;

                res_phase=res_phase+norm(e)^2;

                power_phase= ...
                    power_phase+norm(x_phase)^2;


                %% Time alignment
                S_est=(a_one'*x_align)/ ...
                      (a_one'*a_one);

                e=x_align-a_one*S_est;

                res_align=res_align+norm(e)^2;

                power_align= ...
                    power_align+norm(x_align)^2;


                %% Proposed
                S_est=(a'*x_prop)/(a'*a);

                e=x_prop-a*S_est;

                res_prop=res_prop+norm(e)^2;

                power_prop= ...
                    power_prop+norm(x_prop)^2;

            end

        end

    end


    %% final error
    Err_std=sqrt( ...
        res_std/(power_std+1e-30));

    Err_phase=sqrt( ...
        res_phase/(power_phase+1e-30));

    Err_align=sqrt( ...
        res_align/(power_align+1e-30));

    Err_prop=sqrt( ...
        res_prop/(power_prop+1e-30));


    Err_all(q,:)=[ ...
        Err_std, ...
        Err_phase, ...
        Err_align, ...
        Err_prop];


    %% display result
    fprintf('\n');
    fprintf('Standard STFT           = %.6e\n',Err_std);
    fprintf('Phase-only              = %.6e\n',Err_phase);
    fprintf('Time alignment          = %.6e\n',Err_align);
    fprintf('Proposed                = %.6e\n',Err_prop);

    fprintf('|Standard - Phase|      = %.6e\n', ...
        abs(Err_std-Err_phase));

    fprintf('|Align - Proposed|      = %.6e\n', ...
        abs(Err_align-Err_prop));

end


%% =========================
% Final result table
% =========================

fprintf('\n\n');
fprintf('============================================\n');
fprintf('FINAL REAL-WAVEFORM RESULTS\n');
fprintf('============================================\n');

ResultTable=array2table( ...
    Err_all, ...
    'VariableNames',{ ...
    'Standard', ...
    'PhaseOnly', ...
    'TimeAlign', ...
    'Proposed'}, ...
    'RowNames',{ ...
    'Recording1', ...
    'Recording2', ...
    'Recording3'});

disp(ResultTable);


%% mean error
MeanError=mean(Err_all,1);

fprintf('\nMean errors:\n');

fprintf('Standard    = %.6e\n',MeanError(1));
fprintf('Phase-only  = %.6e\n',MeanError(2));
fprintf('Time-align  = %.6e\n',MeanError(3));
fprintf('Proposed    = %.6e\n',MeanError(4));


%% =========================
% Grouped bar chart
% =========================

figure;

bar(Err_all);

set(gca,'YScale','log');

ylim([1e-4 1e-1]);

yticks([1e-4 1e-3 1e-2 1e-1]);

grid on;

xticklabels({ ...
    'Recording 1', ...
    'Recording 2', ...
    'Recording 3'});

ylabel('Normalized Spatial Model Error');

legend({ ...
    'Standard', ...
    'Phase-only', ...
    'Time-align', ...
    'Proposed'}, ...
    'Location','best');

title('Real-Waveform Validation');


%% =========================
% Fractional delay function
% =========================

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