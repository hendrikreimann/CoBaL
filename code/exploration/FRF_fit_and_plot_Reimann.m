function [Fit,mse] = FRF_fit_and_plot_Reimann ...
  ( ...
    FRF, ...       FRF: vector of FRF values (complex values)
    F, ...         F: vector of FRF frequencies (Hz) - I generally don't include FRF data at frequencies above about 2 Hz
    COH, ...       COH: coherence values
    AvgStim, ...   AvgStim: average stimulus waveform
    AvgResp, ...   AvgResp: average response waveform
    ttl1, ...      ttl1: text title used in plots
    ttl2, ...      ttl2: second text title used in plots
    J, ...         J: Subject moment of inertia about ankle joint
    m, ...         m: Subject mass (not including feet - estimated by subtracting 2.5% of body mass from total mass)
    Hcom, ...      Hcom: CoM height above ankle joint (meters)
    mgh, ...       mgh: Subject mass (excluding feet) times gravity constant times CoM height above ankle joint
    SSstim, ...    SSstim: set to 1 if stimulus was a surface tilt, set to 0 otherwise
    VSstim, ...    VSstim: set to 1 if stimulus was a visual tilt, set to 0 otherwise
    Nfit, ...      Nfit: number of optimization fits to perform using different initial values of parameters
    pp ...         pp: set to 1 to plot results, set to 0 to skip plotting
  )
    %
%FRF_fit_and_plot_Reimann.m - fits simplest possible model including parameters:
%   - Gn (gain factor = sensory weight)
%   - Kp (active stiffness = proportional gain)
%   - Kd (active damping = derivative)
%   - Kf (force/torque feedback gain - torque is mathematically integrated)
%   - Td (time delay)
%
% Inputs: FRF: vector of FRF values (complex values)
%         F: vector of FRF frequencies (Hz) 
%               - I generally don't include FRF data at frequencies above
%               about 2 Hz
%         COH: coherence values
%         AvgStim: average stimulus waveform
%         AvgResp: average response waveform
%         ttl1: text title used in plots
%         ttl2: second text title used in plots
%         J: Subject moment of inertia about ankle joint
%         m: Subject mass (not including feet - estimated by subtracting 2.5% of body mass from total mass)
%         Hcom: CoM height above ankle joint (meters)
%         mgh: Subject mass (excluding feet) times gravity constant times
%                CoM height above ankle joint
%         SSstim: set to 1 if stimulus was a surface tilt, set to 0 otherwise
%         VSstim: set to 1 if stimulus was a visual tilt, set to 0 otherwise
%         Nfit: number of optimization fits to perform using different
%                initial values of parameters
%         pp: set to 1 to plot results, set to 0 to skip plotting
% Outputs: Fit: Structure variable holding best fit parameter values
%          mse: Fit error value of best fit
%

figure_size = [600 800];

samprate=200; % sample rate of time series data - used for plotting
%
%
[Fit,mse]=FRF_fit_Reimann(F,FRF,J,mgh,Nfit); % call to function that performs the optimization fit

gn=Fit.gn;
kp=Fit.kp;
kd=Fit.kd;
kf=Fit.kf;
td=Fit.td;
[gn,kp,kd,kf,td,mse];
%
% Calculate gain and phase curves of fit for display
%        
f=logspace(-2,1,300);
f=f(25:250);
w=2*pi*f;
s=sqrt(-1)*w;
B=(ones(size(w)))./(J*(s.*s)-mgh*ones(size(w)));
NC=kd*s+kp;
FF=(kf*ones(size(w)))./(s); % force/torque feedback includes mathematical integration
TD=(cos(w*td)-sqrt(-1)*sin(w*td));
N=NC.*TD;

frf=(gn*B.*N)./(ones(size(w)) - FF.*N + B.*N);

gfit1=abs(frf);				% fit gain and phase for display
pfit1=180/pi*phase(frf);    % this 'phase' function include phase unwrapping - this may be a legacy Matlab function

if pp==1
    figure('position', [0 0 figure_size])
    subplot(321);loglog(F,abs(FRF),'bo',f,gfit1,'r','markersize',4)
    axis([.01 6 .001 10])
    title(ttl1,'Interpreter','none');
    ylabel('CoM/Stim magitude')
    hold on
    %   error_bars_log(F,abs(TF),ERROR,0.05); % ERROR is the 'r' value used
    %       to calculate confidence limits on gain and phase values. 
%     v1=['J  ',num2str(J),' kg-m^2'];
%     v2=['m  ',num2str(m),' kg'];
%     v3=['h  ',num2str(Hcom),'m'];
%     v4=['mgh  ',num2str(mgh)];
%     text(.002,0.2,v1)
%     text(.002,0.1,v2)
%     text(.002,0.050,v3)
%     text(.002,0.025,v4)
    hold off
    % PhaseError=180/pi*asin((ERROR)./abs(TF)); % calculate error bars on phase
    subplot(323);semilogx(F,180/pi*phase(FRF),'bo',f,pfit1,'r','markersize',4)
%     axis([.01 6 -400 100])
    xlim([.01 6])
    ylabel('CoM/Stim phase')
    hold on
    %   error_bars_log(F,180/pi*phase(TF),PhaseError,0.05);
%     if (SSstim==1)&&(VSstim==0)
%         v1=['Fit: Wp  ',num2str(gn)];      % Fit 1 parameters
%     end
%     if (SSstim==0)&&(VSstim==1)
%         v1=['Fit: Wv  ',num2str(gn)];
%     end
%     if (SSstim==1)&&(VSstim==1)
%         v1=['Fit: Wp+Wv  ',num2str(gn)];
%     end
%     v2=['     Kp  ',num2str(kp),' Nm/rad'];
%     v3=['     Kd  ',num2str(kd),' Nms/rad'];
%     v4=['     Kf  ',num2str(kf),' rad/Nm'];
%     v5=['     Td  ',num2str(td),' s'];
%     v6=['     mse ',num2str(mse)];
%     text(.001,-100,v1)
%     text(.001,-150,v2)
%     text(.001,-200,v3)
%     text(.001,-250,v4)
%     text(.001,-300,v5)
%     text(.001,-350,v6)
    hold off
    subplot(325);semilogx(F,COH,'bo-','markersize',4);
    axis([0.01 6 0 1]);
    ylabel('CoM/Stim Coherence')
    xlabel('Freq (Hz)')
    %   v1=['     mean Coherence  ',num2str(MeanCoh)];
    %   text(.01,0.15,v1)
    tc=(1:length(AvgStim))/samprate;
    subplot(322);plot(tc,AvgStim);
    axis([0 length(AvgStim)/samprate -4 4])
    title(ttl2,'Interpreter','none');
    ylabel('Stim Average (deg)')
    %   v1=['     Stim RMS  ',num2str(StimRMS)];
    %   v2=['     Stim Remnant  ',num2str(StimRemnant)];
    %   text(1,-2,v1)
    %   text(1,-2.5,v2)
    subplot(324);plot(tc,AvgResp-mean(AvgResp),'b');
    axis([0 length(AvgResp)/samprate -4 4])
    ylabel('Resp Average (deg)')
    xlabel('Time (s)')
    %   v1=['     Resp RMS  ',num2str(SwayRMS)];
    %   v2=['     Resp Remnant  ',num2str(SwayRemnant)];
    %   text(1,-2,v1)
    %   text(1,-2.5,v2)
    
    subplot(326)
    xlim([0 2]);
    ylim([0 6]);
        v1=['J  ',num2str(J),' kg-m^2'];
    v2=['m  ',num2str(m),' kg'];
    v3=['h  ',num2str(Hcom),'m'];
    v4=['mgh  ',num2str(mgh)];
    text(0.05,4,v1)
    text(0.05,3,v2)
    text(0.05,2,v3)
    text(0.05,1,v4)

    if (SSstim==1)&&(VSstim==0)
        v1=['Fit: Wp  ',num2str(gn)];      % Fit 1 parameters
    end
    if (SSstim==0)&&(VSstim==1)
        v1=['Fit: Wv  ',num2str(gn)];
    end
    if (SSstim==1)&&(VSstim==1)
        v1=['Fit: Wp+Wv  ',num2str(gn)];
    end
    v2=['     Kp  ',num2str(kp),' Nm/rad'];
    v3=['     Kd  ',num2str(kd),' Nms/rad'];
    v4=['     Kf  ',num2str(kf),' rad/Nm'];
    v5=['     Td  ',num2str(td),' s'];
    v6=['     mse ',num2str(mse)];
    text(1,5.5,v1)
    text(1,4.5,v2)
    text(1,3.5,v3)
    text(1,2.5,v4)
    text(1,1.5,v5)
    text(1,0.5,v6)
    set(gca,'visible','off')
    
    if ~directoryExists('figures')
        mkdir('figures')
    end
    saveas(gcf, ['figures' filesep ttl1], 'pdf')
    
end
end

