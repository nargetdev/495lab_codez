
%% reading the data in
i=1;
for i=1:3
    filename=['Rm_' num2str(i+7) '_good'];
    Data=importdata(filename);
    Data = Data.data;
    Rm_time(:,i)=Data(:,1);
    Rm_voltage(:,i)=Data(:,2);
    Rm_current(:,i)=Data(:,3);
    Rm_distaltach(:,i)=Data(:,4);
    Rm_motortach(:,i)=Data(:,5);
end
i=1;
for i=1:5
    filename=['manual_lefteraser_' num2str(i)];
    Data=importdata(filename);
    Data = Data.data;
    lefteraser_time(:,i)=Data(:,1);
    lefteraser_voltage(:,i)=Data(:,2);
    lefteraser_current(:,i)=Data(:,3);
    lefteraser_distaltach(:,i)=Data(:,4);
    lefteraser_motortach(:,i)=Data(:,5);
end
i=1;

righteraser_time=zeros(10000, 5);
righteraser_voltage=zeros(10000, 5);
righteraser_current=zeros(10000, 5);
righteraser_distaltach=zeros(10000, 5);
righteraser_motortach=zeros(10000, 5);
for i=1:5
    filename=['manual_righteraser_' num2str(i)];
    Data=importdata(filename);
    Data = Data.data;
    righteraser_time_a(:,i)=Data(:,1);
    righteraser_voltage_a(:,i)=Data(:,2);
    righteraser_current_a(:,i)=Data(:,3);
    righteraser_distaltach_a(:,i)=Data(:,4);
    righteraser_motortach_a(:,i)=Data(:,5);
    

    if i==1
        N=1;
    elseif i==2
        N=2840;
    elseif i==3
        N=3700;
    elseif i==4
        N=3389;     
    elseif i==5
        N=1;
    end
  
    righteraser_time(1:(10001-N),i)=righteraser_time_a(N:end,1);
    righteraser_voltage(1:(10001-N),i)=righteraser_voltage_a(N:end,i);
    righteraser_current(1:(10001-N),i)=righteraser_current_a(N:end,i);
    righteraser_distaltach(1:(10001-N),i)=righteraser_distaltach_a(N:end,i);
    righteraser_motortach(1:(10001-N),i)=righteraser_motortach_a(N:end,i);
    
%     righteraser_chopped=righteraser(N:end,:);
%     righteraser_time(:,i)=righteraser_chopped(:,1);
%     righteraser_voltage(:,i)=righteraser_chopped(:,2);
%     righteraser_current(:,i)=righteraser_chopped(:,3);
%     righteraser_distaltach(:,i)=righteraser_chopped(:,4);
%     righteraser_motortach(:,i)=righteraser_chopped(:,5);
end

i=1;
for i=1:4
    if i==1
        filename=['freqsweep_1_day2_bodedata'];
    elseif i==2
        filename=['freqsweep_2_day1_bodedata'];
    elseif i==3
        filename=['freqsweep_2_day2_bodedata'];
    elseif i==4
        filename='freqsweep_zoomed_day2_bodedata';
    end
    Data=importdata(filename);
    Data = Data.data;
    if i<4
        freqsweep_freq(:,i)=Data(:,1);
        freqsweep_current(:,i)=Data(:,2);
        freqsweep_bodemag(:,i)=Data(:,3);
        freqsweep_bodephase(:,i)=Data(:,4);
    elseif i==4
        freqsweepzoom_freq=Data(:,1);
        freqsweepzoom_current=Data(:,2);
        freqsweepzoom_bodemag=Data(:,3);
        freqsweepzoom_bodephase=Data(:,4);       
    end
end

i=1;
for i=1:3
    filename=['onoff_' num2str(i+11) '_good'];
    Data=importdata(filename);
    Data = Data.data;
    onoff_time(:,i)=Data(:,1);
    onoff_voltage(:,i)=Data(:,2);
    onoff_current(:,i)=Data(:,3);
    onoff_distaltach(:,i)=Data(:,4);
    onoff_motortach(:,i)=Data(:,5);
end

scale_distaltach=(pi/30)/(20.8/1000);
scale_motortach=(pi/30)/(7/1000);

%% Constants
% flywheel given
D_f=0.137;
L_f=0.0127;
rho_f=7755;
% shaft given
D_s=3.18E-3;
L_s=0.305;
G_s=7.31E10;
% motor given
L_m=0.002;
J_m=4.23E-5;
R_m_given=1.6;
K_m=0.0097;

J_s=pi*D_s^4/32;
J_f=rho_f*pi*D_f^4*L_f/32;
J_c=J_f+J_m;
K_f=G_s*J_s/L_s;
K=K_f;

%% left eraser running motor

Rm_calculated=abs(Rm_voltage./Rm_current);
Rm=mean(mean(Rm_calculated));

%% motor constant calculation
i=1;
count=1;
for i=[11 9 6 0]
    if i>0
    filename=['onoff_' num2str(i)];
   
    elseif i==0
    filename=['onoff_' num2str(18) '_motorunplug'];    
    end
    
    Data=importdata(filename);
    Data = Data.data;
    Km_time(:,count)=Data(1,1);
    Km_voltage(:,count)=Data(1,2);
    Km_current(:,count)=Data(1,3);
    Km_distaltach(:,count)=Data(1,4);
    %in RPM not rad/s
    %Km_distaltach=Km_distaltach*(30/pi);
    %Km_motortach(:,count)=Data(1,5);
    count=count+1;

end
    X=Km_distaltach*(30/pi);
    Y=Km_voltage-Km_current*Rm;
    Km_fit=polyfit(X,Y,1);
    plot(X,Y, 'o',X,Km_fit(1).*X+Km_fit(2),'r');
    %axis([0 26 0 4])
    xlabel('distaltach')
    ylabel('Vin-Rm*im')
    Km_calculated=Km_fit(1);
    
    
    
%% indictance calculation

%% resonant frequency using frequency sweep 
i=1;
for i=[1:3]
    bla=input('do you want bode plots? 1=yes')
    if bla==1
        figure
        plot(freqsweep_freq(:,i), freqsweep_bodemag(:,i));
    end
    [magMax, magMax_pos] = max(freqsweep_bodemag(:,i));
    resfreq_each(1,i)=freqsweep_freq(magMax_pos,i);
end
w_n=mean(resfreq_each);

%% left eraser
i=1;
a=1;
% adjust for scale factor!
%lefteraser_distaltach=lefteraser_distaltach.*scale_distaltach;

for i=[1:5]
    [pks, pkslocs]=findpeaks(lefteraser_distaltach(:,i),'MinPeakWidth',50,'MinPeakDistance',50);
    bla2=input('do you want lefteraser plots? 1=yes')
    
    if bla2==1
        figure
        plot(lefteraser_time(:,i), lefteraser_distaltach(:,i));
        hold on
        plot(lefteraser_time(pkslocs,i), lefteraser_distaltach(pkslocs,i),'r');
    end 
    
    for a=2:(size(pks)-1)
        r_each(1,i)=0;
        r_each(a,i)=(1/(a+1)).*log(lefteraser_distaltach(pkslocs(1),i)/lefteraser_distaltach(pkslocs(a),i));
    end
    
    r_lefte(1,i)=mean(r_each(:,i));
    tau_lefte(1,i)=lefteraser_time(pkslocs(2),i)-lefteraser_time(pkslocs(1),i);
    
end
    zeta_x_wn_lefte=0.04/0.245;
    %zeta_x_wn_lefte=mean(r_lefte)/mean(tau_lefte);
    w_d_nicholasbosio=3/(lefteraser_time(pkslocs(3))-lefteraser_time(pkslocs(1)));
    
   
%% right eraser
i=1;
a=1;

%righteraser_motortach=righteraser_motortach.*scale_motortach;

for i=[1:5]
    [pks, pkslocs]=findpeaks(righteraser_motortach(:,i),'MinPeakWidth',50,'MinPeakDistance',50);
    bla3=input('do you want righteraser plots? 1=yes')
    if bla3==1
        figure
        plot(righteraser_time(:,i), righteraser_motortach(:,i));
        hold on
        plot(righteraser_time(pkslocs,i), righteraser_motortach(pkslocs,i),'r');
    end
%     for a=2:(size(pks)-1)
%         r_each(1,i)=0;
%         r_each(a,i)=(1/(a+1)).*log(righteraser_distaltach(pkslocs(1),i)/righteraser_distaltach(pkslocs(a),i));
%     end
%     
%     r_righte(1,i)=mean(r_each(:,i));
    r_righte(1,i)=log(righteraser_motortach(pkslocs(1),i)/righteraser_motortach(pkslocs(2),i));
    tau_righte(1,i)=righteraser_time(pkslocs(2),i)-righteraser_time(pkslocs(1),i);
    
end

    zeta_x_wn_righte=mean(r_righte)/mean(tau_righte);
%% applying constant voltage and then unplugging
i=1;

%onoff_distaltach=onoff_distaltach.*scale_distaltach;

for i=[1:3]
    %lm=fitlm(onoff_distaltach(4001:12001,i),onoff_time(4001:12001,i));
    onoff_fit=polyfit(onoff_time(4001:12001,i),onoff_distaltach(4001:12001,i),1);
    %slope=onoff_fit(1);
    bla4=input('do you want onoff plot? 1=yes')
    if bla4==1
        figure
        plot(onoff_time(:,i), onoff_distaltach(:,i),onoff_time(:,i),onoff_fit(1).*onoff_time(:,i)+onoff_fit(2))
        xlabel('time(s)')
        ylabel('distaltach')
        
    end
end



%% damping coefficients
%     B_tr=2*zeta_x_wn_lefte*(J_f+J_s);
%     B_tl=2*zeta_x_wn_righte*(J_m+J_f+J_s);
    B_tr=2*zeta_x_wn_lefte*(J_f);
    B_tl=2*zeta_x_wn_righte*(J_m+J_f);
    % slope=-Bonoff/Jm+2Jf
    J_onoff=J_s+J_m+2*J_f;
    B_onoff=(onoff_fit(1))*J_onoff*-1;
    % B onoff = 2 Bb + Bm
    % B tr = Bs + Bb
    % B tl = Bm + Bs+Bb
        % B tl - Btr = Bm
    B_m=B_tl-B_tr;
    B_b=0.5*(B_onoff-B_m);
    B_s=(B_tr-B_b);

%% MODEL
    R_m=mean(mean(Rm_calculated));
    
    NUM=[K_m*B_s  K_m*K];
    DEN=[(R_m*J_c*J_f)  (R_m*J_c*B_tr+J_f*(K_m^2+B_tl*R_m))  ...
        (R_m*(K*(J_c+J_f)+B_b*B_s)+(K_m^2+R_m*B_m)*B_tr)  ...
        (K*(K_m^2+(B_b+B_m)*R_m))];
        
    SYS=tf(NUM,DEN);
    bode(SYS)
    axis([30 50 -50 0])
    
    
    
%% plots

% bode data vs model
bode(SYS)
axis([30 50 -50 0])
figure
plot(freqsweep_freq*(2*pi),freqsweep_bodemag, 'r')




 
