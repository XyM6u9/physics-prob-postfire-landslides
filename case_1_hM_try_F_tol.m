function [hM]=case_1_hM_try_F_tol(L,L_1,Water_Table,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,beta_degree,t,q,Sink,tol)
%Written by Masood Abdollahi
%q=xlsread('Rainfall','San Gabriel Dam','J367:J2192'); %Infiltration(m/day)
%T=xlsread('Transpiration_85','Sheet3','L367:L2192'); %Transpiration(m/day)
%L=0.6;         %unit m
%L_1=0.15;      %Depth up in which there is no root (m)
%gamaW=9.8;     %unit kn/m3
%EMd=3000;     %Elastic Modulus dry state (kPa)
%EMs=50;      %Elastic Modulus saturate state (kPa)
%mu=0.4;        %Poisson's ratio
%theta_r=0.00;
%theta_s=0.40;
%k_s = .4128;  %saturate hydraulic conductivity (m/day)
%alpha=2.0;   %related to air entery (1/m)
%m=1.26;
%t=length(q); %time (day)
%Sink=T/(L-L_1);
%beta_degree=39; %slope angle in degree
beta=beta_degree*pi/180;
isittrue=0;
%Water_Table=2.00; %Depth of water table below the bottom boundary
L=L+Water_Table;
L_1=L_1+Water_Table;

% while isittrue==0
%     analysmode=input("enter analysmode, either coupled or uncoupled:\n",'s');
%     if analysmode=="coupled" || analysmode=="uncoupled"
%         isittrue=1;
%     end
% end
analysmode="coupled";
%%
MAX_ITER = 10000;  % 
%Transient
N=fix((L-Water_Table)*100+1);
hM=zeros(N,t);
for time=1:1:t
    Summation_1=zeros(N,2);
    Summation_2=zeros(N,2);
    for j=1:1:time
        
        epsnon=1e-5;%0.00001;
        r=1*q(j,1)/k_s*cos(beta);
        % q_eff = min(q(j), k_s * cos(beta));
        % r     = q_eff/(k_s * cos(beta));

        S=Sink(j,1);
        htop1=-epsnon;
        htop2=epsnon-L*cos(beta);
        z=L-epsnon;
        h_l_eps_1=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop1,z,time-j+1,S,beta,analysmode,tol);
        r_1=exp(alpha/2*(htop1+h_l_eps_1))*((htop1-h_l_eps_1)/epsnon+cos(beta));
        h_l_eps_2=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop2,z,time-j+1,S,beta,analysmode,tol);
        r_2=exp(alpha/2*(htop2+h_l_eps_2))*((htop2-h_l_eps_2)/epsnon+cos(beta));

        counter_stop1 = 0;
    while 1
        counter_stop1 = counter_stop1 + 1;
        if counter_stop1 > MAX_ITER
            error('Error: Infinite loop');
        end
        htop=(htop2-htop1)/(r_2-r_1)*(r-r_1)+htop1;
        h_l_eps=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop,z,time-j+1,S,beta,analysmode,tol);
        r_test=exp(alpha/2*(htop+h_l_eps))*((htop-h_l_eps)/epsnon+cos(beta));
        if r_test<(r-epsnon)
            r_2=r_test;
            htop2=htop;
        end
        if r_test>(r+epsnon)
            r_1=r_test;
            htop1=htop;
        end
        if abs(r_test-r)<=epsnon
            break
        end    
    end
    counter_1=1;
    for z=Water_Table:.01:L
        hM_1(counter_1,1)=Heavisidestep(time,j-1)*Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop,z,time-j+1,S,beta,analysmode,tol);
        hM_1(counter_1,2)=z;
        counter_1=counter_1+1;
    end
    Summation_1(:,1)=Summation_1(:,1)+hM_1(:,1);
    Summation_1(:,2)=hM_1(:,2);
    end

    for j=1:1:time-1
            
            epsnon=0.00001;
            r=1*q(j,1)/k_s*cos(beta);
            S=Sink(j,1);
            htop1=-epsnon;
            htop2=epsnon-L*cos(beta);
            z=L-epsnon;
            h_l_eps_1=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop1,z,time-j,S,beta,analysmode,tol);
            r_1=exp(alpha/2*(htop1+h_l_eps_1))*((htop1-h_l_eps_1)/epsnon+cos(beta));
            h_l_eps_2=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop2,z,time-j,S,beta,analysmode,tol);
            r_2=exp(alpha/2*(htop2+h_l_eps_2))*((htop2-h_l_eps_2)/epsnon+cos(beta));
            
            
            counter_stop2 = 0;
            while 1
                counter_stop2 = counter_stop2 + 1;
                if counter_stop2 > MAX_ITER
                    error('Error: Infinite loops');
                end

                htop=(htop2-htop1)/(r_2-r_1)*(r-r_1)+htop1;
                h_l_eps=Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop,z,time-j,S,beta,analysmode,tol);
                r_test=exp(alpha/2*(htop+h_l_eps))*((htop-h_l_eps)/epsnon+cos(beta));
                if r_test<(r-epsnon)
                    r_2=r_test;
                    htop2=htop;
                end
                if r_test>(r+epsnon)
                    r_1=r_test;
                    htop1=htop;
                end
                if abs(r_test-r)<=epsnon
                    break
                end    
            end
            counter_1=1;
            for z=Water_Table:.01:L
                hM_2(counter_1,1)=Heavisidestep(time,j)*Flowhead_tran_unif_tol(L,L_1,gamaW,EMd,EMs,mu,theta_r,theta_s,k_s,alpha,m,htop,z,time-j,S,beta,analysmode,tol);
                hM_2(counter_1,2)=z;
                counter_1=counter_1+1;
            end
            Summation_2(:,1)=Summation_2(:,1)+hM_2(:,1);
            Summation_2(:,2)=hM_2(:,2);
    end
hM(:,time)=Summation_1(:,1)-Summation_2(:,1);
end
hM=[Summation_1(:,2) hM];
%%
%Steady State
%counter_2=1;
%c_2=(r-S*L/k_s+S*L_1/k_s);
%for z=0:0.01:L
%if z>=L_1
%        hss_bar=c_2*(1-exp(-1*alpha*z*cos(beta)))+S*(z-L_1)/k_s/cos(beta)...
%         +exp(-1*alpha*cos(beta)*(z-L_1))*S/alpha/k_s/((cos(beta))^2)-S/alpha/k_s/((cos(beta))^2);
%    else
%        hss_bar=c_2*(1-exp(-1*alpha*z*cos(beta)));
%end
%    hss=1/alpha*log(hss_bar+exp(-1*alpha*z*cos(beta)));
%    hssM(counter_2,1)=hss;
%    hssM(counter_2,2)=z;
%    counter_2=counter_2+1;
%end  
end