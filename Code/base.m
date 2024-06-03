%SIR discrete model-MULTIPLE PATCH
%updated 8-10-19
%during the discussion with Nina on 8-2-19, we decided to test R0. To make
%the latin hypercubic sampling work, here we use "for" loops to test four
%parameters that determine R0: beta of two genotypes, gamma of two
%genotypes, growth rate that relates to the tradeoff of the two genotypes
%and migration that determines the arriving time of the disease. Therefore,
%JJ here made one change before four "for" loops, here are two parameters
%fixed: number of patches is at 40

%8-19-19 JJ realized that under fixed transmission rates, growth rate and
%migration have complicated influences (not linear, sometimes, increase of growth rate with fixed migration produce cycling, but further increase reduce)
%therefore, instead of using growth ratio, JJ here define growth rates
%separately. This time, fixed tranmission rates, but only change growth
%rates and migration rate.

%%8-20-19 JJ and NF talked to check the combined effects of 1. migration
%%rate vs growth rate; 2. migration vs transmission rate; 3. transmisson
%%rate vs growth rate. JJ thoughts to rerun the results to produce the
%%heatmap by only controling two factors but controling the other two.
%%On 8-25-19, JJ checked how transmission rate vs. migration influences the
%%results. Here, JJ fixed the ratio of transmission rate as 20 here and
%%change the transmission rate of RT [0.000003:0.000002:0.000011],
%%migration rate is [0.000001:0.000009:0.000037]. All other parameters are:
%%growth rate = 0.2 and 0.16.

%%define strays for all variables in the system
%9-4-19 based on the discussion between JJ and NF, we want to check the
%dampen pattern beyond the number of cycling to make sure that is stable
%cycling instead of a gradually dampening pattern.
%Since WT and RT may have different cycling pattern: e.g., one dies out but
%the other can keep cycling by itself, so we also want to check the average
%population size
%since patch number matters for Tran_RT and Growth results, we plan to plot
%the R0 with time goes by in different patches.


%%%based on the discussion between JJ and Nina on 9-17-19, 
timesteps = 100000;
N= 40;
N1=15;
N2=30;
N3=40;
Time = timesteps * N;

trans_ratio=18:0.25:25; %10
% trans_WT=0.00001:0.000045:0.00019;
mig_ratio=0.0000005:0.0000013:0.000037;

for i=1:length(trans_ratio) %transmission rate of RT
    for jj=1:length(mig_ratio)  %ratio of transmission rate that determines the transmission of WT
       
S_RT = zeros(1,Time);
I_RT = zeros(1,Time);
R_RT = zeros(1,Time);

S_WT = zeros(1,Time);
I_WT = zeros(1,Time);
R_WT = zeros(1,Time);

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

% S_RT(1:timesteps:end)
% S_WT(1:timesteps:end)
% I_RT(1:timesteps:end)
% I_WT(1:timesteps:end)
% R_RT(1:timesteps:end)
% R_WT(1:timesteps:end)

%%define parameters

%%take out the parameters that will be used for loops: r_WT, beta_RR and
%%_WR, ratio of growth rate, migration rate (growth ratio is fixed at 1.35)


 K = 3000;
 r_RT = 0.15;     % natural growth rate in S
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
 beta_RW = trans_ratio(i)*beta_RR;%infection rate from RT to WT
 beta_WW = trans_ratio(i)*beta_WR;%infection rate fro WT to WT
 gamma_RT = 0.05;%Recovery rate 
 gamma_WT = 0.05;
 
%%Migration between patches
V=eye(N);
%%%added row
V_a=V(N,:);
V1=[V_a;V];
V1(N+1,:)=[];

mig_matrix = V1*mig_ratio(jj);
% mig_matrix_S = V1*0.000;
% mig_matrix(1,N)=0;
% mig_matrix(3,2)=0;
% mig_matrix(4,3)=0;
% mig_matrix(5,4)=0;
% mig_matrix(6,5)=0;
% mig_matrix(7,6)=0;
% mig_matrix(8,7)=0;
% mig_matrix(9,8)=0;
% mig_matrix(10,9)=0;

% mig_matrix
% mig_matrixaa = mig_matrixa';
% mig_matrixab = mig_matrixa + mig_matrixa';
% mig_matrixb = [0 1 1 0 0; 1 0 1 0 0; 1 1 0 1 1;0 0 1 0 1;0 0 1 1 0]*mig(k);         
% mig_matrixc = [0 1 1 0 0; 1 0 1 1 0; 1 1 0 0 1;0 1 0 0 1;0 0 1 1 0]*mig(k);
% mig_matrixd = [0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
% mig_matrixe = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 0;0 1 1 0 0]*mig(k);
% mig_matrixf = [0 1 0 0 1; 1 0 1 0 1; 0 1 0 1 1;0 0 1 0 1;1 1 1 1 0]*mig(k);
% mig_matrixg = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 0;1 0 1 0 1;0 1 0 1 0]*mig(k);
% mig_matrixh = [0 1 1 0 0; 1 0 1 1 1; 1 1 0 1 1;0 1 1 0 1;0 1 1 1 0]*mig(k);
% mig_matrixi = [0 1 1 1 0; 1 0 1 0 1; 1 1 0 1 1;1 0 1 0 1;0 1 1 1 0]*mig(k);
% mig_matrixj = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 0;1 1 1 0 0]*mig(k);
% mig_matrixk = [0 1 1 1 1; 1 0 1 1 1; 1 1 0 1 1;1 1 1 0 1;1 1 1 1 0]*mig(k);




t1=timesteps;
for t=2:t1
     for j=0:N-1
         
        if(I_RT(t-1+timesteps*j)+I_WT(t-1+timesteps*j)<9e-05)
            I_RT(t-1+timesteps*j)=0;
            I_WT(t-1+timesteps*j)=0;
        end
          
        Ni = S_RT(t-1+timesteps*j) + I_RT(t-1+timesteps*j) + R_RT(t-1+timesteps*j) + S_WT(t-1+timesteps*j) + I_WT(t-1+timesteps*j)+ R_WT(t-1+timesteps*j);    %population size at patch i: sum up all Ss and Is
        
        %define newborns for each genotype susceptibles
        %S1S1
%         LS1S1 = check_zero(S_RT(t-1+timesteps*j)*S_RT(t-1+timesteps*j)/Ni*r_RT*(1-Ni/K));
%         %S1S2
%         LS1S2 = check_zero(S_RT(t-1+timesteps*j)*S_WT(t-1+timesteps*j)/(2*Ni)*(r_RT+r_WT)/2*(1-Ni/K));
%         %S1I1
%         LS1I1 = check_zero(S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)/Ni*(r_RT+rd_RT)/2*(1-Ni/K));
%         %S1I2
%         LS1I2 = check_zero(S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)/(2*Ni)*(r_RT+rd_WT)/2*(1-Ni/K));
%         
%         LS1R1 = check_zero(S_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)/Ni*(r_RT+rr_RT)/2*(1-Ni/K));   
%         %S1R2
%         LS1R2 = check_zero(S_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(2*Ni)*(r_RT+rr_WT)/2*(1-Ni/K)); 
%         %S2S2
%         LS2S2 = check_zero(S_WT(t-1+timesteps*j)*S_WT(t-1+timesteps*j)/Ni*r_WT*(1-Ni/K));
%         %S2I1
%         LS2I1 = check_zero(S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)/(2*Ni)*(r_WT+rd_RT)/2*(1-Ni/K));
%         %S2I2
%         LS2I2 = check_zero(S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)/Ni*(r_WT+rd_WT)/2*(1-Ni/K));
%         
%         LS2R1 = check_zero(S_WT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)/(2*Ni)*(r_WT+rr_RT)/2*(1-Ni/K));
%         %S2R2
%         LS2R2 = check_zero(S_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(1*Ni)*(r_WT+rr_WT)/2*(1-Ni/K));
%         %I1I1
%         LI1I1 = check_zero(I_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)/Ni*rd_RT*(1-Ni/K));
%         %I1I2
%         LI1I2 = check_zero(I_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)/(2*Ni)*(rd_RT+rd_WT)/2*(1-Ni/K));
%         
%         LI1R1 = check_zero(I_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)/(1*Ni)*(rd_RT+rr_RT)/2*(1-Ni/K));   
%         %I1R2
%         LI1R2 = check_zero(I_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(2*Ni)*(rd_RT+rr_WT)/2*(1-Ni/K));
%         %I2I2
%         LI2I2 = check_zero(I_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)/Ni*rd_WT*(1-Ni/K));
%         
%         LI2R1 = check_zero(I_WT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)/(2*Ni)*(rd_WT+rr_RT)/2*(1-Ni/K));   
%         %I2R2
%         LI2R2 = check_zero(I_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(1*Ni)*(rd_WT+rr_WT)/2*(1-Ni/K));  
%         %R1R1
%         LR1R1 = check_zero(R_RT(t-1+timesteps*j)*R_RT(t-1+timesteps*j)/(1*Ni)*(rr_RT+rr_RT)/2*(1-Ni/K));   
%         %R1R2
%         LR1R2 = check_zero(R_RT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(2*Ni)*(rr_RT+rr_WT)/2*(1-Ni/K));   
%         %R2R2
%         LR2R2 = check_zero(R_WT(t-1+timesteps*j)*R_WT(t-1+timesteps*j)/(1*Ni)*(rr_WT+rr_WT)/2*(1-Ni/K));  

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
%          tmpS1=0.0000*t;
        
        S_RT(t+timesteps*j)= S_RT(t-1+timesteps*j)+ (LS1S1 + LS1S2 + LS1I1 + LS1I2 + LS1R1 + LS1R2 + LS2I1 + LS2R1 + LI1I1 + LI1I2 + LI1R1 + LI1R2 + LI2R1 + LR1R1 + LR1R2)*(1-Ni/K) - beta_RR*S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)- beta_WR*S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)- mu_RT*S_RT(t-1+timesteps*j)+ tmpS1 ; %=dS1/dt
        
        

        tmpS2 = dot(mig_matrix(j+1,:), S_WT(t-1:timesteps:end)) - S_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1)) ;
%         tmpS2=0.0000*t;
        
        S_WT(t+timesteps*j)= S_WT(t-1+timesteps*j)+ (LS1S2 + LS1I2 + LS1R2 + LS2S2 + LS2I1 + LS2I2 + LS2R1 + LS2R2 + LI1I2 + LI1R2 + LI2I2 + LI2R1 + LI2R2 + LR1R2 + LR2R2)*(1-Ni/K) - beta_RW*S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)- beta_WW*S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)- mu_WT*S_WT(t-1+timesteps*j)+ tmpS2; %=dS2/dt
       
        
        tmpI1 = dot(mig_matrix(j+1,:), I_RT(t-1:timesteps:end)) - I_RT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));
%         if (t==20)
%             tmpI1 = 0.0;
%        
%         else
%             tmpI1 =0;
%         end
     
        I_RT(t+timesteps*j) = I_RT(t-1+timesteps*j)+ beta_RR*S_RT(t-1+timesteps*j)*I_RT(t-1+timesteps*j)+ beta_WR* S_RT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)-(alpha_RT+mu_RT+gamma_RT)*I_RT(t-1+timesteps*j)+ tmpI1;%=dI1/dt
        
        
        tmpI2 = dot(mig_matrix(j+1,:), I_WT(t-1:timesteps:end)) - I_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1)) ;
         
%        if (t==20)
%             tmpI2 = 0.0;
%        
%         else
%             tmpI2 =0;
%         end
%         
        I_WT(t+timesteps*j) = I_WT(t-1+timesteps*j)+ beta_RW*S_WT(t-1+timesteps*j)*I_RT(t-1+timesteps*j) + beta_WW*S_WT(t-1+timesteps*j)*I_WT(t-1+timesteps*j)-(alpha_WT+mu_WT+gamma_WT)*I_WT(t-1+timesteps*j) + tmpI2;%=dI2/dt
        
        
        tmpR1 =dot(mig_matrix(j+1,:), R_RT(t-1:timesteps:end)) - R_RT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));
%         tmpR1 = 0;
        
        R_RT(t+timesteps*j) = R_RT(t-1+timesteps*j) + gamma_RT*I_RT(t-1+timesteps*j) - mu_RT*R_RT(t-1+timesteps*j)+tmpR1;

        
        
        tmpR2 = dot(mig_matrix(j+1,:), R_WT(t-1:timesteps:end)) - R_WT(t-1+timesteps*j) * sum(mig_matrix(:,j+1));
%         tmpR2 = 0;
        
        R_WT(t+timesteps*j) = R_WT(t-1+timesteps*j) + gamma_WT*I_WT(t-1+timesteps*j) - mu_WT*R_WT(t-1+timesteps*j)+tmpR2;
        
       
        
%         if (I_RT(t+timesteps*j)<I_RT(t-1+timesteps*j) && I_RT(t+timesteps*j)<20)
%             I_RT(t+timesteps*j)=0;
%         
%         end
%         
%         if (I_WT(t+timesteps*j)<I_WT(t-1+timesteps*j) && I_WT(t+timesteps*j)<20)
%             %||(I_WT(t+timesteps*j)<10 && dot(mig_matrix(j+1,:), I_WT(t-1:timesteps:end))==0)
%             I_WT(t+timesteps*j)=0;
%         
%         end
        
%         t
%        I_RT(t+timesteps*j)
%        
%        I_WT(t+timesteps*j) 
%         if(I_RT(t+timesteps*j)<4.89e-10)
%             I_RT(t+timesteps*j)=0;
%         end
% %        
%         if(I_WT(t+timesteps*j)<4.89e-10)
%             I_WT(t+timesteps*j)=0;
%         end 

%         if(I_RT(t+timesteps*j)<4.89e-5)
%             I_RT(t+timesteps*j)=0;
%         end
% %        
%         if(I_WT(t+timesteps*j)<4.89e-5)
%             I_WT(t+timesteps*j)=0;
%         end 

%         if(I_RT(t+timesteps*j)<4.89e-6)
%             I_RT(t+timesteps*j)=0;
%         end
% %        
%         if(I_WT(t+timesteps*j)<4.89e-6)
%             I_WT(t+timesteps*j)=0;
%         end 

% 
   end
% t    
% S_RT(t:timesteps:end)
% S_WT(t:timesteps:end)
% I_RT(t:timesteps:end)
% I_WT(t:timesteps:end)
% R_RT(t:timesteps:end)
% R_WT(t:timesteps:end)
end

%define the cycles of each combinations
%for patch #1, we use S_WT(1:timesteps); for patch #15, we use
%S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps); for patch #30, we use
%S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps); for patch #40, we use
%S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)

%patch 1

%patch 1

burnin1=S_WT(timesteps-50000:timesteps); %only collect the last 50000 steps
burninR1=S_RT(timesteps-50000:timesteps);
R_0_1 = (beta_RR*S_RT(1:timesteps)+ beta_RW*S_WT(1:timesteps))/(gamma_RT+mu_RT+alpha_RT);


[~,locs1] = findpeaks (-burnin1.');
[~,locs1_1] = findpeaks (burnin1.');
% cycles1=0;
if (min(length(locs1),length(locs1_1))==0 | min(length(locs1),length(locs1_1))==1)
%     cycles1=1;
    dampen1=[];
else
    lengt=min(length(locs1),length(locs1_1));
    dampen1 = abs(burnin1(locs1_1(1:lengt))-burnin1(locs1(1:lengt)));
end

if(isempty(dampen1))
    cycles1=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_1=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampen1)>1)
    cycles1=length(dampen1);
    p_1=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampen1);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_1=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    


ave_pop1 = mean(burnin1+I_WT(timesteps-50000:timesteps)+R_WT(timesteps-50000:timesteps));
I_1=mean(I_WT(timesteps-50000:timesteps));

[~,locsR1] = findpeaks (-burninR1.');
[~,locsR1_1] = findpeaks (burninR1.');
% cycles1=0;
if (min(length(locsR1),length(locsR1_1))==0 | min(length(locsR1),length(locsR1_1))==1)
%     cycles1=1;
    dampenR1=[];
else
    lengt=min(length(locsR1),length(locsR1_1));
    dampenR1 = abs(burninR1(locsR1_1(1:lengt))-burninR1(locsR1(1:lengt)));
end

if(isempty(dampenR1))
    cyclesR1=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_R1=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampenR1)>1)
    cyclesR1=length(dampenR1);
    p_R1=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampenR1);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_R1=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    

ave_popR1 = mean(burninR1+I_RT(timesteps-50000:timesteps)+R_RT(timesteps-50000:timesteps));
I_R1=mean(I_RT(timesteps-50000:timesteps));
%patch 15

R_0_15 = (beta_RR*S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+ beta_RW*S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps))/(gamma_RT+mu_RT+alpha_RT);

burninR2=S_RT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps);
[~,locsR2] = findpeaks (-burninR2.');
[~,locsR2_1] = findpeaks (burninR2.');
% cycles1=0;
if (min(length(locsR2),length(locsR2_1))==0 | min(length(locsR2),length(locsR2_1))==1)
%     cycles1=1;
    dampenR2=[];
else
    lengt=min(length(locsR2),length(locsR2_1));
    dampenR2 = abs(burninR2(locsR2_1(1:lengt))-burninR2(locsR2(1:lengt)));
end

if(isempty(dampenR2))
    cyclesR2=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_R2=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampenR2)>1)
    cyclesR2=length(dampenR2);
    p_R2=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampenR2);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_R2=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    

ave_popR2 = mean(burninR2+I_RT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps)+R_RT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps));
I_R2=mean(I_RT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps));

burnin2=S_WT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps);
[~,locs2] = findpeaks (-burnin2.');
[~,locs2_1] = findpeaks (burnin2.');
% cycles1=0;
if (min(length(locs2),length(locs2_1))==0 | min(length(locs2),length(locs2_1))==1)
%     cycles1=1;
    dampen2=[];
else
    lengt=min(length(locs2),length(locs2_1));
    dampen2 = abs(burnin2(locs2_1(1:lengt))-burnin2(locs2(1:lengt)));
end

if(isempty(dampen2))
    cycles2=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_2=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampen2)>1)
    cycles2=length(dampen2);
    p_2=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampen2);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_2=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    
ave_pop2 = mean(burnin2+I_WT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps)+R_WT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps));
I_2=mean(I_WT((N1-1)*timesteps+timesteps-50000:(N1-1)*timesteps+timesteps));
%patch 30

R_0_30 = (beta_RR*S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+ beta_RW*S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps))/(gamma_RT+mu_RT+alpha_RT);

burnin3=S_WT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps);
[~,locs3] = findpeaks (-burnin3.');
[~,locs3_1] = findpeaks (burnin3.');
% cycles1=0;
if (min(length(locs3),length(locs3_1))==0 | min(length(locs3),length(locs3_1))==1)
%     cycles1=1;
    dampen3=[];
else
    lengt=min(length(locs3),length(locs3_1));
    dampen3 = abs(burnin3(locs3_1(1:lengt))-burnin3(locs3(1:lengt)));
end

if(isempty(dampen3))
    cycles3=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen3)==1)
%     cycles3=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_3=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampen3)>1)
    cycles3=length(dampen3);
    p_3=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampen3);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_3=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    
ave_pop3 = mean(burnin3+I_WT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps)+R_WT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps));
I_3=mean(I_WT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps));

burninR3=S_RT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps);
[~,locsR3] = findpeaks (-burninR3.');
[~,locsR3_1] = findpeaks (burninR3.');
% cycles1=0;
if (min(length(locsR3),length(locsR3_1))==0 | min(length(locsR3),length(locsR3_1))==1)
%     cycles1=1;
    dampenR3=[];
else
    lengt=min(length(locsR3),length(locsR3_1));
    dampenR3 = abs(burninR3(locsR3_1(1:lengt))-burninR3(locsR3(1:lengt)));
end

if(isempty(dampenR3))
    cyclesR3=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_R3=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampenR3)>1)
    cyclesR3=length(dampenR3);
    p_R3=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampenR3);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_R3=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    
ave_popR3 = mean(burninR3+I_RT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps)+R_RT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps));
I_R3=mean(I_RT((N2-1)*timesteps+timesteps-50000:(N2-1)*timesteps+timesteps));

%patch 40
R_0_40 = (beta_RR*S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+ beta_RW*S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps))/(gamma_RT+mu_RT+alpha_RT);

burnin4=S_WT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps);
[~,locs4] = findpeaks (-burnin4.');
[~,locs4_1] = findpeaks (burnin4.');
% cycles1=0;
if (min(length(locs4),length(locs4_1))==0 | min(length(locs4),length(locs4_1))==1)
%     cycles1=1;
    dampen4=[];
else
    lengt=min(length(locs4),length(locs4_1));
    dampen4 = abs(burnin4(locs4_1(1:lengt))-burnin4(locs4(1:lengt)));
end

if(isempty(dampen4))
    cycles4=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_4=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampen4)>1)
    cycles4=length(dampen4);
    p_4=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampen4);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_4=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    


ave_pop4 = mean(burnin4+I_WT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps)+R_WT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps));
I_4=mean(I_WT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps));

burninR4=S_RT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps);
[~,locsR4] = findpeaks (-burninR4.');
[~,locsR4_1] = findpeaks (burninR4.');
% cycles1=0;
if (min(length(locsR4),length(locsR4_1))==0 | min(length(locsR4),length(locsR4_1))==1)
%     cycles1=1;
    dampenR4=[];
else
    lengt=min(length(locsR4),length(locsR4_1));
    dampenR4 = abs(burninR4(locsR4_1(1:lengt))-burninR4(locsR4(1:lengt)));
end

if(isempty(dampenR4))
    cyclesR4=0; %picked up area is not enough to include one complete cycle
end

% if(length(dampen1)==1)
%     cycles1=1; %picked up area only contains one cycle, which we still does not call cycling pattern we are interested in
% end

p_R4=0; %by default, cycles are not enough to determine if dampened exists
if (length(dampenR4)>1)
    cyclesR4=length(dampenR4);
    p_R4=1; %by default, once cycles number >=2, those cycles are perminant periodic cycles
    d=diff(dampenR4);
    s=sign(d);
    p = all(s==-1);
    if(p==1)
        p_R4=-1; % pick out dampened (very conservative, so if p_1=-1, it must dampened, in other cases, it may still dampen, but I did not capture)
    end
end    
    



ave_popR4 = mean(burninR4+I_RT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps)+R_RT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps));
I_R4=mean(I_RT((N3-1)*timesteps+timesteps-50000:(N3-1)*timesteps+timesteps));


XXW = [i,jj,cycles1,cycles2,cycles3,cycles4,ave_pop1,ave_pop2,ave_pop3,ave_pop4,I_1,I_2,I_3,I_4];
dlmwrite('SIRmodel_WT_cycles_avepop_aveInfect_Trans_vs_Mig_9_28_19.csv',XXW,'-append');

XXR = [i,jj,cyclesR1,cyclesR2,cyclesR3,cyclesR4,ave_popR1,ave_popR2,ave_popR3,ave_popR4,I_R1,I_R2,I_R3,I_R4];
dlmwrite('SIRmodel_RT_cycles_avepop_aveInfect_Trans_vs_Mig_9_28_19.csv',XXR,'-append');

% YY_P1 = [i,jj,R_0_1];
% dlmwrite('SIRmodel_R_0_1_Trans_vs_Mig_9_9_19.csv',YY_P1,'-append');
% 
% YY_P15 = [i,jj,R_0_15];
% dlmwrite('SIRmodel_R_0_15_Trans_vs_Mig_9_9_19.csv',YY_P15,'-append');
% 
% YY_P30 = [i,jj,R_0_30];
% dlmwrite('SIRmodel_R_0_30_Trans_vs_Mig_9_9_19.csv',YY_P30,'-append');
% 
% YY_P40 = [i,jj,R_0_40];
% dlmwrite('SIRmodel_R_0_40_Trans_vs_Mig_9_9_19.csv',YY_P40,'-append');


% S_RT1=[i,jj,k,h,S_RT(timesteps-35000:timesteps)];
% S_WT1=[i,jj,k,h,S_WT(timesteps-35000:timesteps)];
% I_RT1=[i,jj,k,h,I_RT(timesteps-35000:timesteps)];
% I_WT1=[i,jj,k,h,I_WT(timesteps-35000:timesteps)];
% R_RT1=[i,jj,k,h,R_RT(timesteps-35000:timesteps)];
% R_WT1=[i,jj,k,h,R_WT(timesteps-35000:timesteps)];
% 
% XX1=[S_RT1;S_WT1;I_RT1;I_WT1;R_RT1;R_WT1];
% dlmwrite('SIRmodel_patch1_8_16_19.csv',XX1,'-append');
% 
% 
% S_RT2=[i,jj,k,h,S_RT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% S_WT2=[i,jj,k,h,S_WT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% I_RT2=[i,jj,k,h,I_RT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% I_WT2=[i,jj,k,h,I_WT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% R_RT2=[i,jj,k,h,R_RT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% R_WT2=[i,jj,k,h,R_WT((N1-1)*timesteps+timesteps-35000:(N1-1)*timesteps+timesteps)];
% 
% XX2=[S_RT2;S_WT2;I_RT2;I_WT2;R_RT2;R_WT2];
% dlmwrite('SIRmodel_patchN1_8_16_19.csv',XX2,'-append');
% 
% 
% S_RT3=[i,jj,k,h,S_RT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% S_WT3=[i,jj,k,h,S_WT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% I_RT3=[i,jj,k,h,I_RT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% I_WT3=[i,jj,k,h,I_WT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% R_RT3=[i,jj,k,h,R_RT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% R_WT3=[i,jj,k,h,R_WT((N2-1)*timesteps+timesteps-35000:(N2-1)*timesteps+timesteps)];
% 
% XX3=[S_RT3;S_WT3;I_RT3;I_WT3;R_RT3;R_WT3];
% 
% dlmwrite('SIRmodel_patchN2_8_16_19.csv',XX3,'-append');
% 
% 
% S_RT4=[i,jj,k,h,S_RT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% S_WT4=[i,jj,k,h,S_WT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% I_RT4=[i,jj,k,h,I_RT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% I_WT4=[i,jj,k,h,I_WT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% R_RT4=[i,jj,k,h,R_RT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% R_WT4=[i,jj,k,h,R_WT((N3-1)*timesteps+timesteps-35000:(N3-1)*timesteps+timesteps)];
% 
% XX4=[S_RT4;S_WT4;I_RT4;I_WT4;R_RT4;R_WT4];
% dlmwrite('SIRmodel_patchN3_8_16_19.csv',XX4,'-append');
 % plot everything of the two genotypes
% figure;
subplot(4,4,1)
h1=plot(1:timesteps,S_RT(1:timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
% xlabel 'Time';
% ylabel 'RT'

hold on
% subplot(2,4,2) 
h1=plot(1:timesteps,I_RT(1:timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
h1=plot(1:timesteps,R_RT(1:timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,2)
h1=plot(1:timesteps,S_WT(1:timesteps),'--b'); 
% xlabel 'Time';
% ylabel 'WT'

% subplot(2,4,5) 
hold on
h1=plot(1:timesteps,I_WT(1:timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
h1=plot(1:timesteps,R_WT(1:timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,3) 
h1=plot(1:timesteps,S_RT(1:timesteps)+I_RT(1:timesteps)+R_RT(1:timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'
% 
% subplot(2,4,8) 
hold on
h1=plot(1:timesteps,S_WT(1:timesteps)+I_WT(1:timesteps)+R_WT(1:timesteps),'--black','Linewidth',2);
% xlabel 'Time';
% ylabel 'Total'

hold on
h1=plot(1:timesteps,S_RT(1:timesteps)+I_RT(1:timesteps)+R_RT(1:timesteps)+ S_WT(1:timesteps)+I_WT(1:timesteps)+R_WT(1:timesteps),'-m','Linewidth',2);


subplot(4,4,4)

h1=plot(1:length(R_0_1),R_0_1);
% xlabel 'Time';
% ylabel 'R0'


subplot(4,4,5)
h1=plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
% xlabel 'Time';
% ylabel 'RT'

hold on
% subplot(2,4,2) 
h1=plot(1:timesteps,I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
h1=plot(1:timesteps,R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,6)
h1=plot(1:timesteps,S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--b'); 
% xlabel 'Time';
% ylabel 'WT'

% subplot(2,4,5) 
hold on
h1=plot(1:timesteps,I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
h1=plot(1:timesteps,R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,7) 
h1=plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'
% 
% subplot(2,4,8) 
hold on
h1=plot(1:timesteps,S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'--black','Linewidth',2);
% xlabel 'Time';
% ylabel 'Total'

hold on
h1=plot(1:timesteps,S_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_RT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+ S_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+I_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps)+R_WT((N1-1)*timesteps+1:(N1-1)*timesteps+timesteps),'-m','Linewidth',2);




subplot(4,4,8)
h1=plot(1:length(R_0_15),R_0_15);
% xlabel 'Time';
% ylabel 'R0'


subplot(4,4,9)
h1=plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
% xlabel 'Time';
% ylabel 'RT'

hold on
% subplot(2,4,2) 
h1=plot(1:timesteps,I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
h1=plot(1:timesteps,R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff

subplot(4,4,10)
h1=plot(1:timesteps,S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--b'); 
% xlabel 'Time';
% ylabel 'WT'

% subplot(2,4,5) 
hold on
h1=plot(1:timesteps,I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
h1=plot(1:timesteps,R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,11) 
h1=plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
h1=plot(1:timesteps,S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'--black','Linewidth',2);
% xlabel 'Time';
% ylabel 'Total'

hold on
h1=plot(1:timesteps,S_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_RT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+ S_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+I_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps)+R_WT((N2-1)*timesteps+1:(N2-1)*timesteps+timesteps),'-m','Linewidth',2);

% 
% 
%
subplot(4,4,12)
h1=plot(1:length(R_0_30),R_0_30);
% xlabel 'Time';
% ylabel 'R0'
% 
subplot(4,4,13)
h1=plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-b');
%lgd.FontSize = 8;
%axis([0, 100, 0, 50])
% xlabel 'Time';
% ylabel 'RT'

hold on
% subplot(2,4,2) 
h1=plot(1:timesteps,I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-r');
% xlabel 'Time';
% ylabel 'RT in I'

% subplot(2,4,3) 
hold on
h1=plot(1:timesteps,R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-g');
% xlabel 'Time';
% ylabel 'RT in R'
% legend('S_R','I_R','R_R');
% legend boxoff
% 
subplot(4,4,14)
h1=plot(1:timesteps,S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--b'); 
% xlabel 'Time';
% ylabel 'WT'

% subplot(2,4,5) 
hold on
h1=plot(1:timesteps,I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--r');
% xlabel 'Time';
% ylabel 'WT in I'

% subplot(2,4,6) 
hold on
h1=plot(1:timesteps,R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--g');
% xlabel 'Time';
% ylabel 'WT in R'

subplot(4,4,15) 
h1=plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-black','Linewidth',2);
% xlabel 'Time';
% ylabel 'RT total'

% subplot(2,4,8) 
hold on
h1=plot(1:timesteps,S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'--black','Linewidth',2);
% xlabel 'Time';
% ylabel 'Total'

hold on
h1=plot(1:timesteps,S_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_RT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+ S_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+I_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps)+R_WT((N3-1)*timesteps+1:(N3-1)*timesteps+timesteps),'-m','Linewidth',2);

subplot(4,4,16)
h1=plot(1:length(R_0_40),R_0_40);
% xlabel 'Time';
% ylabel 'R0'

% tightfig;

filename = sprintf('Trans_Ratio %d vs Mig %d_9_28_19.png',i,jj);
% counter=[i,jj,k,h];
% figname = strcat('SIR',num2str(counter),'_8_13_19.png') ;
% saveas(gcf, figname); 
% clf;

% name=[num2str(counter),'_SIR_8_13_19.png'];
saveas(h1,filename);
pause(1);
clf;
pause(1);

            end
        end
    

