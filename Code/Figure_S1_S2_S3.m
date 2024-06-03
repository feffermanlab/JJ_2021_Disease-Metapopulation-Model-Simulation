%SIR discrete model-MULTIPLE PATCH
% =========================================================================================================================================================

% Matlab code to generate Figure S1-S3 in Appendix 1 for manuscript titled "The Eco-Epidemiological Dynamics of Evolutionary Rescue in a Host Metapopulation"

% =========================================================================================================================================================

timesteps = 100000;
N= 40;
N1=15;
N2=30;
N3=40;
Time = timesteps * N;


S_RT = zeros(1,Time);
I_RT = zeros(1,Time);
R_RT = zeros(1,Time);

S_WT = zeros(1,Time);
I_WT = zeros(1,Time);
R_WT = zeros(1,Time);
R0 = zeros(1,Time);

%%define initial values at N patches
S_RT(1)= 200;
S_WT(1)= 200;
I_RT(1)= 5;
I_WT(1)= 5;
R_RT(1)= 0;
R_WT(1)= 0;

for j=1:N-1
    S_RT(1+timesteps*j)=200;
    S_WT(1+timesteps*j)=200;
    I_RT(1+timesteps*j)=0;
    I_WT(1+timesteps*j)=0;
    R_RT(1+timesteps*j)=0;
    R_WT(1+timesteps*j)=0;
end


%%define parameters
 K = 3000;
 r_RT = 0.16;     % natural growth rate in S
 r_WT = 0.2;                                 
 rd_RT = 0.01;   % natural growth rate in I
 rd_WT = 0.01;
 rr_RT = 0.2;    %natural growth rate in R
 rr_WT = 0.2;    
 mu_RT = 0.0005; %natural mortality rate 
 mu_WT = 0.0005;
 alpha_RT = 0.05; %disease-induced mortality rate
 alpha_WT = 0.05;
 beta_RR = 0.000005;%infection rate from RT to RT
 beta_WR = 0.000005; %infection rate from WT to RT
 beta_RW = 0.0001;%infection rate from RT to WT
 beta_WW = 0.0001;%infection rate fro WT to WT
 gamma_RT = 0.05;%Recovery rate 
 gamma_WT = 0.05;
 
%%Migration between patches
V=eye(N);
%%%added row
V_a=V(N,:);
V1=[V_a;V];
V1(N+1,:)=[];
%V1(40,1)=1;
mig_matrix = V1*0.00001;
mig_matrix(35,1) = 0.00001/3;
mig_matrix(36,1) = 0.00001/3;
mig_matrix(37,1) = 0.00001/3;

t1=timesteps;

for t=2:t1
     for j=0:N-1
         
        if(I_RT(t-1+timesteps*j)+I_WT(t-1+timesteps*j)<9e-05)
            I_RT(t-1+timesteps*j)=0;
            I_WT(t-1+timesteps*j)=0;
        end
          
        
        Ni = S_RT(t-1+timesteps*j) + I_RT(t-1+timesteps*j) + R_RT(t-1+timesteps*j) + S_WT(t-1+timesteps*j) + I_WT(t-1+timesteps*j)+ R_WT(t-1+timesteps*j);    %population size at patch i: sum up all Ss and Is
        
        %define newborns from all combinations across genotypes (Wild or
        %robust) and stages (S, I or R)

        LS1S1 = Ni*(S_RT(t-1+timesteps*j)*S_RT(t-1+timesteps*j))/(Ni*Ni)*r_RT;
        %S1S2
        LS1S2 = Ni*(S_RT(t-1+timesteps*j)*S_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(r_RT+r_WT)/2;
        %S1I1
        LS1I1 = Ni*(S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)*2)/(Ni*Ni)*(r_RT+rd_RT)/2;
        %S1I2
        LS1I2 = Ni*(S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(r_RT+rd_WT)/2;
        
        LS1R1 = Ni*(S_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)*2)/(Ni*Ni)*(r_RT+rr_RT)/2;   
        %S1R2
        LS1R2 = Ni*(S_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(r_RT+rr_WT)/2; 
        %S2S2
        LS2S2 = Ni*(S_WT(t-1+timesteps*j)*S_WT(t-1+timesteps*j))/(Ni*Ni)*r_WT;
        %S2I1
        LS2I1 = Ni*(S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(r_WT+rd_RT)/2;
        %S2I2
        LS2I2 = Ni*(S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)*2)/(Ni*Ni)*(r_WT+rd_WT)/2;
        
        LS2R1 = Ni*(S_WT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(r_WT+rr_RT)/2;
        %S2R2
        LS2R2 = Ni*(S_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)*2)/(Ni*Ni)*(r_WT+rr_WT)/2;
        %I1I1
        LI1I1 = Ni*(I_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j))/(Ni*Ni)*rd_RT;
        %I1I2
        LI1I2 = Ni*(I_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(rd_RT+rd_WT)/2;
        
        LI1R1 = Ni*(I_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)*2)/(Ni*Ni)*(rd_RT+rr_RT)/2;   
        %I1R2
        LI1R2 = Ni*(I_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(rd_RT+rr_WT)/2;
        %I2I2
        LI2I2 = Ni*(I_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j))/(Ni*Ni)*rd_WT;
        
        LI2R1 = Ni*(I_WT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(rd_WT+rr_RT)/2;   
        %I2R2
        LI2R2 = Ni*(I_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)*2)/(Ni*Ni)*(rd_WT+rr_WT)/2;  
        %R1R1
        LR1R1 = Ni*(R_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j))/(Ni*Ni)*(rr_RT+rr_RT)/2;   
        %R1R2
        LR1R2 = Ni*(R_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)*2)/(Ni*Ni)/2*(rr_RT+rr_WT)/2;   
        %R2R2
        LR2R2 = Ni*(R_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j))/(Ni*Ni)*(rr_WT+rr_WT)/2;  
        
        
        tmpS1 = dot(mig_matrix(j+1,:), S_RT(t-1:timesteps:end)) - S_RT(t-1+timesteps*j) * sum(mig_matrix(:,j+1)) ;% migration between patches
         
        
        S_RT(t+timesteps*j)= S_RT(t-1+timesteps*j)+ (LS1S1 + LS1S2 + LS1I1 + LS1I2 + LS1R1 + LS1R2 + LS2I1 + LS2R1 + LI1I1 + LI1I2 + LI1R1 + LI1R2 + LI2R1 + LR1R1 + LR1R2)*(1-Ni/K) - beta_RR*S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)- beta_WR*S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)- mu_RT*S_RT(t-1+timesteps*j)+ tmpS1 ; %=dS1/dt
        
        

        tmpS2 = dot(mig_matrix(j+1,:), S_WT(t-1:timesteps:end)) - S_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1)) ;
       
        
        S_WT(t+timesteps*j)= S_WT(t-1+timesteps*j)+ (LS1S2 + LS1I2 + LS1R2 + LS2S2 + LS2I1 + LS2I2 + LS2R1 + LS2R2 + LI1I2 + LI1R2 + LI2I2 + LI2R1 + LI2R2 + LR1R2 + LR2R2)*(1-Ni/K) - beta_RW*S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)- beta_WW*S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)- mu_WT*S_WT(t-1+timesteps*j)+ tmpS2; %=dS2/dt
       
        
        tmpI1 = dot(mig_matrix(j+1,:), I_RT(t-1:timesteps:end)) - I_RT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));

        
        I_RT(t+timesteps*j) = I_RT(t-1+timesteps*j)+ beta_RR*S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)+ beta_WR* S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)-(alpha_RT+mu_RT+gamma_RT)*I_RT(t-1+timesteps*j)+ tmpI1;%=dI1/dt
        
        
        tmpI2 = dot(mig_matrix(j+1,:), I_WT(t-1:timesteps:end)) - I_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1)) ;
         

%         
        I_WT(t+timesteps*j) = I_WT(t-1+timesteps*j)+ beta_RW*S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j) + beta_WW*S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)-(alpha_WT+mu_WT+gamma_WT)*I_WT(t-1+timesteps*j) + tmpI2;%=dI2/dt
        
        
        tmpR1 =dot(mig_matrix(j+1,:), R_RT(t-1:timesteps:end)) - R_RT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));

        
        R_RT(t+timesteps*j) = R_RT(t-1+timesteps*j) + gamma_RT*I_RT(t-1+timesteps*j) - mu_RT*R_RT(t-1+timesteps*j)+tmpR1;

        
        
        tmpR2 = dot(mig_matrix(j+1,:), R_WT(t-1:timesteps:end)) - R_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));

        
        R_WT(t+timesteps*j) = R_WT(t-1+timesteps*j) + gamma_WT*I_WT(t-1+timesteps*j) - mu_WT*R_WT(t-1+timesteps*j)+tmpR2;
        

       % R0(t+timesteps*j) = (I_RT(t+timesteps*j) + I_WT(t+timesteps*j))*(S_RT(t+timesteps*j)*beta_RR + S_WT(t+timesteps*j)*beta_WW)/((gamma_RT + mu_RT + alpha_RT)*(I_RT(t+timesteps*j) + I_WT(t+timesteps*j)) - tmpI1 - tmpI2);
     
     end

end


R0_no = (S_RT * beta_RR + S_WT * beta_WW)/(gamma_RT + mu_RT + alpha_RT);
 % % plot everything of the two genotypes

subplot(4,4,1)
plot(1:timesteps,S_RT(1:timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'RT'

hold on
% subplot(2,4,2) 
plot(1:timesteps,I_RT(1:timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
plot(1:timesteps,R_RT(1:timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
%legend('S_R','I_R','R_R');
%egend boxoff

subplot(4,4,2)
plot(1:timesteps,S_WT(1:timesteps),'--b'); 
xlabel 'Time';
ylabel 'WT'

% subplot(2,4,5) 
hold on
plot(1:timesteps,I_WT(1:timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
plot(1:timesteps,R_WT(1:timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,3) 
plot(1:timesteps,S_RT(1:timesteps)+I_RT(1:timesteps)+R_RT(1:timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
plot(1:timesteps,S_WT(1:timesteps)+I_WT(1:timesteps)+R_WT(1:timesteps),'--black','Linewidth',2);
xlabel 'Time';
ylabel 'Total'

hold on
plot(1:timesteps,S_RT(1:timesteps)+I_RT(1:timesteps)+R_RT(1:timesteps)+ S_WT(1:timesteps)+I_WT(1:timesteps)+R_WT(1:timesteps),'-r','Linewidth',2);

subplot(4,4,4)
plot(1:timesteps,R0_no(1:timesteps),'-m');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'R0'



subplot(4,4,5)
plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'RT'

hold on
% subplot(2,4,2) 
plot(1:timesteps,I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
plot(1:timesteps,R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,6)
plot(1:timesteps,S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--b'); 
xlabel 'Time';
ylabel 'WT'

% subplot(2,4,5) 
hold on
plot(1:timesteps,I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
plot(1:timesteps,R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,7) 
plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
plot(1:timesteps,S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--black','Linewidth',2);
xlabel 'Time';
ylabel 'Total'

hold on
plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+ S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-r','Linewidth',2);

subplot(4,4,8)
plot(1:timesteps,R0_no((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-m');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'R0'



subplot(4,4,9)
plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'RT'

hold on
% subplot(2,4,2) 
plot(1:timesteps,I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
plot(1:timesteps,R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,10)
plot(1:timesteps,S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--b'); 
xlabel 'Time';
ylabel 'WT'

% subplot(2,4,5) 
hold on
plot(1:timesteps,I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
plot(1:timesteps,R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,11) 
plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
plot(1:timesteps,S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--black','Linewidth',2);
xlabel 'Time';
ylabel 'Total'

hold on
plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+ S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-r','Linewidth',2);


subplot(4,4,12)

plot(1:timesteps,R0_no((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-m');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'R0'



subplot(4,4,13)
plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'RT'

hold on
% subplot(2,4,2) 
plot(1:timesteps,I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
plot(1:timesteps,R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,14)
plot(1:timesteps,S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--b'); 
xlabel 'Time';
ylabel 'WT'

% subplot(2,4,5) 
hold on
plot(1:timesteps,I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
plot(1:timesteps,R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,15) 
plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
plot(1:timesteps,S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--black','Linewidth',2);
xlabel 'Time';
ylabel 'Total'

hold on
plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+ S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-r','Linewidth',2);

subplot(4,4,16)
plot(1:timesteps,R0_no((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-m');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
xlabel 'Time';
ylabel 'R0'


