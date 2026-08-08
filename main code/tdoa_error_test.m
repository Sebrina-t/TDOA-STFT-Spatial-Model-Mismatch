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


fprintf('========== True TDOA ==========\n');

fprintf('True relative delay (samples):\n');
disp(Ds);

fprintf('Chirp range: %.1f -> %.1f kHz\n', ...
    f1/1e3,f2/1e3);


%% =========================================================
% TDOA estimation error sweep
%
% err is the signed error parameter.
% The microphone-dependent perturbation is
% linearly distributed from 0 to -err samples.
% ==========================================================

err_list=-0.3:0.025:0.3;

Ne=length(err_list);


%% =========================================================
% Result arrays
% ==========================================================

Err_std_true=zeros(1,Ne);
Err_prop_true=zeros(1,Ne);

Err_std_est=zeros(1,Ne);
Err_prop_est=zeros(1,Ne);


%% =========================================================
% Generate true microphone signals
% ==========================================================

Smic_std=zeros(M,N);

for m=1:M

    u=t-t_start-Ts_m(m);

    Smic_std(m,:)=source_chirp( ...
        u,f1,mu,Ts);

end


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


%% =========================================================
% True steering vectors
% ==========================================================

Atrue=zeros(K,M);

for k=1:K

    for m=1:M

        Atrue(k,m)=exp( ...
            -1j*2*pi*f_stft(k)*Ds(m)/fs);

    end

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
% Same interior-frame criterion used in Fig. 1 and Fig. 2
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
% TDOA estimation-error sweep
% ==========================================================

for e=1:Ne

    err=err_list(e);


    %% -------------------------------------------------------
    % Estimated TDOA
    % --------------------------------------------------------

    Derr=linspace(0,err,M);

    Dhat=Ds+Derr;

    tau_hat=Dhat/fs;


    %% -------------------------------------------------------
    % Temporal alignment using estimated TDOA
    %
    % Physical propagation remains determined by the
    % true delays Ts_m.
    % --------------------------------------------------------

    Smic_align=zeros(M,N);

    for m=1:M

        t_sample=t+tau_hat(m);

        u=t_sample-t_start-Ts_m(m);

        Smic_align(m,:)=source_chirp( ...
            u,f1,mu,Ts);

    end


    %% -------------------------------------------------------
    % STFT of temporally aligned signals
    % --------------------------------------------------------

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


    %% -------------------------------------------------------
    % Phase-referenced aligned representation
    %
    % Restore spatial phase using estimated TDOA
    % --------------------------------------------------------

    Srestore=zeros(K,L,M);

    for m=1:M

        phase_restore=exp( ...
            -1j*2*pi*f_stft*Dhat(m)/fs);

        Srestore(:,:,m)= ...
            Salign(:,:,m).*phase_restore;

    end


    %% -------------------------------------------------------
    % Estimated steering vectors
    % --------------------------------------------------------

    Ahat=zeros(K,M);

    for k=1:K

        for m=1:M

            Ahat(k,m)=exp( ...
                -1j*2*pi*f_stft(k)*Dhat(m)/fs);

        end

    end


    %% -------------------------------------------------------
    % Error accumulators
    % --------------------------------------------------------

    res_std_true=0;
    res_prop_true=0;

    res_std_est=0;
    res_prop_est=0;

    power_std=0;
    power_prop=0;


    %% -------------------------------------------------------
    % Spatial-model errors
    % --------------------------------------------------------

    for k=1:K

        a_true=Atrue(k,:).';
        a_hat=Ahat(k,:).';

        for l=1:L

            if mask(k,l)

                x_std=squeeze(Sstd(k,l,:));

                x_prop=squeeze(Srestore(k,l,:));


                %% ============================================
                % Error relative to TRUE spatial model
                % =============================================

                S_std_true= ...
                    (a_true'*x_std)/(a_true'*a_true);

                S_prop_true= ...
                    (a_true'*x_prop)/(a_true'*a_true);


                e_std_true= ...
                    x_std-a_true*S_std_true;

                e_prop_true= ...
                    x_prop-a_true*S_prop_true;


                %% ============================================
                % Error relative to ESTIMATED spatial model
                % =============================================

                S_std_est= ...
                    (a_hat'*x_std)/(a_hat'*a_hat);

                S_prop_est= ...
                    (a_hat'*x_prop)/(a_hat'*a_hat);


                e_std_est= ...
                    x_std-a_hat*S_std_est;

                e_prop_est= ...
                    x_prop-a_hat*S_prop_est;


                %% ============================================
                % Residual energies
                % =============================================

                res_std_true= ...
                    res_std_true+norm(e_std_true)^2;

                res_prop_true= ...
                    res_prop_true+norm(e_prop_true)^2;


                res_std_est= ...
                    res_std_est+norm(e_std_est)^2;

                res_prop_est= ...
                    res_prop_est+norm(e_prop_est)^2;


                %% Signal energies

                power_std= ...
                    power_std+norm(x_std)^2;

                power_prop= ...
                    power_prop+norm(x_prop)^2;

            end

        end

    end


    %% -------------------------------------------------------
    % Final normalized errors
    % --------------------------------------------------------

    Err_std_true(e)=sqrt( ...
        res_std_true/(power_std+1e-30));

    Err_prop_true(e)=sqrt( ...
        res_prop_true/(power_prop+1e-30));


    Err_std_est(e)=sqrt( ...
        res_std_est/(power_std+1e-30));

    Err_prop_est(e)=sqrt( ...
        res_prop_est/(power_prop+1e-30));

end


%% =========================================================
% Print numerical results
% ==========================================================

fprintf('\n========== Results ==========\n');

fprintf(['Error     Std-True      PhaseRef-True   ' ...
         'Std-Estimated   PhaseRef-Estimated\n']);


for e=1:Ne

    fprintf('%6.3f   %.4e   %.4e   %.4e   %.4e\n', ...
        err_list(e), ...
        Err_std_true(e), ...
        Err_prop_true(e), ...
        Err_std_est(e), ...
        Err_prop_est(e));

end


%% =========================================================
% Locate zero-error condition
% ==========================================================

[~,e0]=min(abs(err_list));

fprintf('\n========== Zero TDOA-Error Condition ==========\n');

fprintf('Standard STFT error = %.6e\n', ...
    Err_std_true(e0));

fprintf('Phase-referenced aligned error = %.6e\n', ...
    Err_prop_true(e0));


%% =========================================================
% PAPER FIGURE 3
% Robustness to TDOA estimation error
%
% Error is evaluated relative to the TRUE spatial model.
% ==========================================================

%% =========================================================
% PAPER FIGURE 3
% Robustness to TDOA estimation error
% Logarithmic y-axis
% ==========================================================

figure(3);
clf;
hold on;


% Avoid zero values on logarithmic axis
% This only affects plotting, not the calculated results
%% =========================================================
% PAPER FIGURE 3
% TDOA estimation robustness - LOG Y AXIS
% ==========================================================

% Numerical floor for log plotting only
Err_std_plot  = max(Err_std_true, 1e-14);
Err_prop_plot = max(Err_prop_true,1e-14);

figure(3);
clf;
hold on;

plot( ...
    err_list, ...
    Err_std_plot, ...
    '-o', ...
    'LineWidth',1.5, ...
    'MarkerSize',6);

plot( ...
    err_list, ...
    Err_prop_plot, ...
    '--s', ...
    'LineWidth',1.5, ...
    'MarkerSize',6);

%% FORCE logarithmic Y axis
set(gca,'YScale','log');

grid on;
box on;

xlabel('TDOA Estimation Error, \epsilon_D (samples)');
ylabel('Normalized Spatial-Model Mismatch');

legend( ...
    'Standard STFT', ...
    'Phase-referenced aligned', ...
    'Location','southwest');

xlim([-0.3 0.3]);

% This is the important part
ylim([1e-14 1e-1]);
yticks([1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e-1]);

xticks(-0.3:0.1:0.3);

set(gcf, ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 12 8]);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'LineWidth',1);

%% Export
exportgraphics( ...
    gcf, ...
    'Fig3_TDOA_Estimation_Robustness.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    gcf, ...
    'Fig3_TDOA_Estimation_Robustness.png', ...
    'Resolution',600);
%% continuous source model
function y=source_chirp(u,f1,mu,Ts)

y=zeros(size(u));

idx=(u>=0)&(u<Ts);

u0=u(idx);

phase=2*pi*( ...
    f1*u0+0.5*mu*u0.^2);

y(idx)=cos(phase);

end
