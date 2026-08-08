clear;clc;close all;

%% =========================================================
% Basic parameters
% ==========================================================

fs=192e3;

T=0.05;
N=round(T*fs);

n=0:N-1;
t=n/fs;

M=4;

t_start=15e-3;

% Propagation delay to Mic 1
Tref=3e-3;


%% =========================================================
% Source parameters
% ==========================================================

fc=50e3;
Ts=4e-3;

mu=10e6;

f1=fc-mu*Ts/2;
f2=fc+mu*Ts/2;


%% =========================================================
% True TDOA
% ==========================================================

Dmax=2;

Ds=linspace(0,-Dmax,M);

tau=Ds/fs;

Ts_m=Tref+tau;


fprintf('========== Baseline Comparison ==========\n');

fprintf('True relative delays (samples):\n');
disp(Ds);

fprintf('Chirp rate = %.1f MHz/s\n',mu/1e6);

fprintf('Chirp range = %.1f -> %.1f kHz\n', ...
    f1/1e3,f2/1e3);


%% =========================================================
% Generate standard microphone signals
%
% All microphones are evaluated on the same
% absolute-time sampling grid.
% ==========================================================

Smic=zeros(M,N);

for m=1:M

    u=t-t_start-Ts_m(m);

    Smic(m,:)=source_chirp( ...
        u,f1,mu,Ts);

end


%% =========================================================
% Ideal source-relative temporal alignment
% ==========================================================

Smic_align=zeros(M,N);

for m=1:M

    t_sample=t+tau(m);

    u=t_sample-t_start-Ts_m(m);

    Smic_align(m,:)=source_chirp( ...
        u,f1,mu,Ts);

end


%% =========================================================
% Check ideal time alignment
% ==========================================================

u_ref=t-t_start-Tref;

idx_check=(u_ref>1/fs) & ...
          (u_ref<Ts-1/fs);

align_error=zeros(1,M);

for m=1:M

    align_error(m)=norm( ...
        Smic_align(m,idx_check)- ...
        Smic_align(1,idx_check)) / ...
        (norm(Smic_align(1,idx_check))+1e-12);

end

fprintf('\nMaximum time-domain alignment error = %.6e\n', ...
    max(align_error));


%% =========================================================
% STFT parameters
% ==========================================================

Nw=64;

Noverlap=round(0.75*Nw);

Nfft=256;

win_stft=hann(Nw,'periodic');


%% =========================================================
% Standard synchronous STFT
% ==========================================================

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


%% =========================================================
% STFT of temporally aligned signals
% ==========================================================

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


%% =========================================================
% Theoretical steering vectors
% ==========================================================

A=zeros(K,M);

for k=1:K

    for m=1:M

        A(k,m)=exp( ...
            -1j*2*pi*f_stft(k)*Ds(m)/fs);

    end

end


%% =========================================================
% Representation 2:
% Phase-only compensation
%
% The spatial phase is removed, but the
% source-relative STFT time displacement remains.
% ==========================================================

Sphase=zeros(K,L,M);

for m=1:M

    phase_comp=exp( ...
        1j*2*pi*f_stft*Ds(m)/fs);

    Sphase(:,:,m)= ...
        Sstd(:,:,m).*phase_comp;

end


%% =========================================================
% Representation 4:
% Phase-referenced aligned representation
%
% First perform temporal alignment, then map the
% aligned STFT back to the original spatial-phase
% coordinate using the steering-vector phase.
% ==========================================================

Sphase_ref=zeros(K,L,M);

for m=1:M

    phase_restore=exp( ...
        -1j*2*pi*f_stft*Ds(m)/fs);

    Sphase_ref(:,:,m)= ...
        Salign(:,:,m).*phase_restore;

end


%% =========================================================
% Target-active TF mask
% ==========================================================

mag_ref=abs(Sstd(:,:,1));

threshold=max(mag_ref(:))*10^(-30/20);

mask_freq=mag_ref>threshold;


%% =========================================================
% Remove source start/end frames
%
% Same interior-frame criterion used in
% Fig. 1, Fig. 2 and Fig. 3.
% ==========================================================

t_source=t_stft-(t_start+Tref);

margin=Nw/(2*fs);

mask_time=(t_source>=margin) & ...
          (t_source<=Ts-margin);

mask=false(K,L);

for l=1:L

    if mask_time(l)

        mask(:,l)=mask_freq(:,l);

    end

end


%% =========================================================
% Error accumulators
% ==========================================================

res_std=0;
res_phase=0;
res_align=0;
res_phase_ref=0;

power_std=0;
power_phase=0;
power_align=0;
power_phase_ref=0;


%% =========================================================
% Calculate spatial-model mismatch
% ==========================================================

for k=1:K

    % Original steering vector
    a=A(k,:).';

    % Steering vector after phase removal / time alignment
    a_one=ones(M,1);


    for l=1:L

        if mask(k,l)

            x_std=squeeze( ...
                Sstd(k,l,:));

            x_phase=squeeze( ...
                Sphase(k,l,:));

            x_align=squeeze( ...
                Salign(k,l,:));

            x_phase_ref=squeeze( ...
                Sphase_ref(k,l,:));


            %% ================================================
            % 1. Standard synchronous STFT
            % ================================================

            S_est= ...
                (a'*x_std)/(a'*a);

            e= ...
                x_std-a*S_est;

            res_std= ...
                res_std+norm(e)^2;

            power_std= ...
                power_std+norm(x_std)^2;


            %% ================================================
            % 2. Phase-only compensation
            % ================================================

            S_est= ...
                (a_one'*x_phase)/(a_one'*a_one);

            e= ...
                x_phase-a_one*S_est;

            res_phase= ...
                res_phase+norm(e)^2;

            power_phase= ...
                power_phase+norm(x_phase)^2;


            %% ================================================
            % 3. Ideal temporal alignment
            % ================================================

            S_est= ...
                (a_one'*x_align)/(a_one'*a_one);

            e= ...
                x_align-a_one*S_est;

            res_align= ...
                res_align+norm(e)^2;

            power_align= ...
                power_align+norm(x_align)^2;


            %% ================================================
            % 4. Phase-referenced aligned representation
            % ================================================

            S_est= ...
                (a'*x_phase_ref)/(a'*a);

            e= ...
                x_phase_ref-a*S_est;

            res_phase_ref= ...
                res_phase_ref+norm(e)^2;

            power_phase_ref= ...
                power_phase_ref+norm(x_phase_ref)^2;

        end

    end

end


%% =========================================================
% Final normalized errors
% ==========================================================

Err_std=sqrt( ...
    res_std/(power_std+1e-30));

Err_phase=sqrt( ...
    res_phase/(power_phase+1e-30));

Err_align=sqrt( ...
    res_align/(power_align+1e-30));

Err_phase_ref=sqrt( ...
    res_phase_ref/(power_phase_ref+1e-30));


%% =========================================================
% Print numerical results
% ==========================================================

fprintf('\n========== Spatial-Model Errors ==========\n');

fprintf('Standard STFT        = %.6e\n', ...
    Err_std);

fprintf('Phase-only           = %.6e\n', ...
    Err_phase);

fprintf('Time-aligned         = %.6e\n', ...
    Err_align);

fprintf('Phase-referenced     = %.6e\n', ...
    Err_phase_ref);


%% =========================================================
% Theoretical equivalence checks
%
% Expected:
%
% E_standard ~= E_phase-only
%
% E_time-aligned ~= E_phase-referenced
% ==========================================================

fprintf('\n========== Equivalence Checks ==========\n');

fprintf('|Standard - Phase-only| = %.6e\n', ...
    abs(Err_std-Err_phase));

fprintf('|Time-aligned - Phase-referenced| = %.6e\n', ...
    abs(Err_align-Err_phase_ref));


%% =========================================================
% PAPER FIGURE 4
% Comparison of four STFT representations
% ==========================================================

Err=[ ...
    Err_std ...
    Err_phase ...
    Err_align ...
    Err_phase_ref];


% Numerical floor for logarithmic plotting only
Err_plot=max(Err,1e-14);


figure(4);
clf;

bar(Err_plot);


set(gca,'YScale','log');

grid on;
box on;


xticks(1:4);

xticklabels({ ...
    'Standard', ...
    'Phase-only', ...
    'Time-aligned', ...
    'Phase-referenced'});


ylabel('Normalized Spatial Model Error');


ylim([1e-14 1e-1]);


set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 12 8]);


set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',1);


%% =========================================================
% Export Fig. 4
% ==========================================================

exportgraphics( ...
    gcf, ...
    'Fig4_Representation_Comparison.pdf', ...
    'ContentType','vector');


exportgraphics( ...
    gcf, ...
    'Fig4_Representation_Comparison.png', ...
    'Resolution',600);


%% =========================================================
% Continuous chirp model
% ==========================================================

function y=source_chirp(u,f1,mu,Ts)

y=zeros(size(u));

idx=(u>=0)&(u<Ts);

u0=u(idx);

phase=2*pi*( ...
    f1*u0+0.5*mu*u0.^2);

y(idx)=cos(phase);

end