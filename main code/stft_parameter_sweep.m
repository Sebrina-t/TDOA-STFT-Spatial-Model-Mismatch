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
% Chirp parameters
% ==========================================================

fc=50e3;
Ts=4e-3;

mu_list=[0 2 5 10 15 20]*1e6;


%% =========================================================
% STFT parameters
% ==========================================================

Nw_list=[32 64 128 256];

Nfft=256;


%% =========================================================
% Maximum relative TDOA
% ==========================================================

Dmax_list=[0.5 1 2 4 6 8];


%% =========================================================
% Result matrices
% ==========================================================

Nmu=length(mu_list);
Nwin=length(Nw_list);
Nd=length(Dmax_list);

Err_std=zeros(Nmu,Nwin,Nd);


fprintf('========== Parameter Sweep ==========\n');

fprintf('Chirp rates = ');
fprintf('%.1f ',mu_list/1e6);
fprintf('MHz/s\n');

fprintf('Window lengths = ');
fprintf('%d ',Nw_list);
fprintf('samples\n');

fprintf('Maximum TDOA = ');
fprintf('%.1f ',Dmax_list);
fprintf('samples\n');


%% =========================================================
% Parameter sweep
% ==========================================================

for a=1:Nmu

    mu=mu_list(a);

    % Chirp is centred at 50 kHz
    f1=fc-mu*Ts/2;
    f2=fc+mu*Ts/2;

    fprintf('\nChirp rate = %.1f MHz/s',mu/1e6);
    fprintf('   %.1f -> %.1f kHz\n', ...
        f1/1e3,f2/1e3);


    for b=1:Nwin

        Nw=Nw_list(b);

        Noverlap=round(0.75*Nw);

        win_stft=hann(Nw,'periodic');


        for q=1:Nd

            Dmax=Dmax_list(q);


            %% Relative delays
            Ds=linspace(0,-Dmax,M);

            tau=Ds/fs;

            Ts_m=Tref+tau;


            %% ================================================
            % Generate standard microphone signals
            % Same absolute-time sampling grid for all mics
            % =================================================

            Smic_std=zeros(M,N);

            for m=1:M

                u=t-t_start-Ts_m(m);

                Smic_std(m,:)=source_chirp( ...
                    u,f1,mu,Ts);

            end


            %% ================================================
            % Standard synchronous STFT
            % =================================================

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


            %% ================================================
            % Theoretical steering vectors
            % =================================================

            A=zeros(K,M);

            for k=1:K

                for m=1:M

                    A(k,m)=exp( ...
                        -1j*2*pi*f_stft(k)*Ds(m)/fs);

                end

            end


            %% ================================================
            % Active TF mask
            % =================================================

            mag_ref=abs(Sstd(:,:,1));

            threshold=max(mag_ref(:))*10^(-30/20);

            mask_freq=mag_ref>threshold;


            %% ================================================
            % Remove start/end frames
            %
            % IMPORTANT:
            % This keeps the analysis focused on frames whose
            % complete STFT window lies inside the chirp.
            % =================================================

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


            %% ================================================
            % Spatial-model mismatch
            % =================================================

            residual_std=0;
            power_std=0;


            for k=1:K

                a_vec=A(k,:).';

                for l=1:L

                    if mask(k,l)

                        x_std=squeeze( ...
                            Sstd(k,l,:));


                        % Best rank-one source coefficient
                        S_est_std= ...
                            (a_vec'*x_std) / ...
                            (a_vec'*a_vec);


                        % Projection residual
                        e_std= ...
                            x_std-a_vec*S_est_std;


                        residual_std= ...
                            residual_std+norm(e_std)^2;

                        power_std= ...
                            power_std+norm(x_std)^2;

                    end

                end

            end


            Err_std(a,b,q)=sqrt( ...
                residual_std/(power_std+1e-30));

        end

    end

end


%% =========================================================
% Publication style
% ==========================================================

marker_list={'o','s','^','d'};
style_list={'-','--','-.',':'};


%% =========================================================
% FIGURE 1
% Effect of maximum TDOA
%
% Fixed:
% chirp rate = 10 MHz/s
%
% Curves:
% Nw = 32, 64, 128, 256
% ==========================================================

mutest=10e6;

[~,a0]=min(abs(mu_list-mutest));


figure(1);
clf;
hold on;


for b=1:Nwin

    plot( ...
        Dmax_list, ...
        squeeze(Err_std(a0,b,:)), ...
        'LineStyle',style_list{b}, ...
        'Marker',marker_list{b}, ...
        'LineWidth',1.5, ...
        'MarkerSize',6);

end


grid on;
box on;

xlabel('Maximum TDOA (samples)');
ylabel('Normalized Spatial Model Error');

legend( ...
    'N_w = 32', ...
    'N_w = 64', ...
    'N_w = 128', ...
    'N_w = 256', ...
    'Location','northwest');

xlim([0 8]);
ylim([0 0.35]);


set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 12 8]);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',1);


exportgraphics( ...
    gcf, ...
    'Fig1_TDOA_Effect.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    gcf, ...
    'Fig1_TDOA_Effect.png', ...
    'Resolution',600);


%% =========================================================
% FIGURE 2
% Effect of STFT window length
%
% Fixed:
% chirp rate = 10 MHz/s
% maximum TDOA = 2 samples
% ==========================================================

mutest=10e6;
Dtest=2;

[~,a0]=min(abs(mu_list-mutest));
[~,q0]=min(abs(Dmax_list-Dtest));


figure(2);
clf;


plot( ...
    Nw_list, ...
    squeeze(Err_std(a0,:,q0)), ...
    '-o', ...
    'LineWidth',1.5, ...
    'MarkerSize',6);


grid on;
box on;

xlabel('STFT Window Length (samples)');
ylabel('Normalized Spatial Model Error');

xlim([0 300]);
ylim([0 0.09]);


set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 12 8]);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',1);


exportgraphics( ...
    gcf, ...
    'Fig2_Window_Length_Effect.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    gcf, ...
    'Fig2_Window_Length_Effect.png', ...
    'Resolution',600);


%% =========================================================
% Print exact paper values
% ==========================================================

fprintf('\n========== FIGURE 1 DATA ==========\n');

fprintf('Fixed chirp rate = %.1f MHz/s\n', ...
    mu_list(a0)/1e6);

for b=1:Nwin

    fprintf('Nw = %d: ',Nw_list(b));

    fprintf('%.6f ', ...
        squeeze(Err_std(a0,b,:)));

    fprintf('\n');

end


fprintf('\n========== FIGURE 2 DATA ==========\n');

fprintf('Fixed chirp rate = %.1f MHz/s\n', ...
    mu_list(a0)/1e6);

fprintf('Fixed Dmax = %.1f samples\n', ...
    Dmax_list(q0));

for b=1:Nwin

    fprintf('Nw = %d, Error = %.6f\n', ...
        Nw_list(b), ...
        Err_std(a0,b,q0));

end


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