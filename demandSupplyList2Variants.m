clc
close all 
clear all 

%% VARIABLES 

% Input (LAUNCH "Run" TO ENTER INPUT INTO WINDOW) 
prompt = {'T : number of periods', ...
    'N : size of the population', ...
    'Average number of contact per person per period at peak infectiousness (Variant 1)', ...
    'probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness (Variant 1)', ...
    'Average number of contact per person per period at peak infectiousness (Variant 2)', ...
    'probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness (Variant 2)', ...
    'number of persons infected at time 0 (Variant 1)',...
    'number of persons infected at time 0 (Variant 2)',...
    'number of person who have recovered at time 0 (i.e: person who are resistant (recovered, vaccinated, ...)',...
    'h0 : ratio of Infected that are treated in hospital if capacity is enough',...
    'ahV1 : ratio of Infected treated in hospital that recovered (Variant 1)' ,...
    'aaV1 : ratio of Infected treated at home that recovered (Variant 1)',...
    'ahV2 : ratio of Infected treated in hospital that recovered (Variant 2)' ,...
    'aaV2 : ratio of Infected treated at home that recovered (Variant 2)',...
    'hospital capacity (number of bed)'};
dlgtitle = 'Inputs';
dims = [1 120];
definput = {'300','47e6','10',...
'0.20','10','0.20','16','16','0','0.04','0.96','0.99','0.96','0.99','1.6e5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);  
 
T = str2num(answer{1}) ;  % number of periods
N = str2num(answer{2}) ; % size of the population
averageContactV1 = str2num(answer{3}) ; % Average number of contact per person per period at peak infectiousness
probaTransmissionV1 = str2num(answer{4}) ; % probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness
averageContactV2 = str2num(answer{5}) ; 
probaTransmissionV2 = str2num(answer{6}) ;
I0V1 = str2num(answer{7}) ; % number of persons infected at time 0
I0V2 = str2num(answer{8}) ; % number of persons infected at time 0
R0 = str2num(answer{9}) ; % number of person who have recovered at time 0 (i.e: person who are resistant (recovered, vaccinated, ...) 
h0 = str2num(answer{10}) ;
ahV1 = str2num(answer{11}) ; % ratio of Infected treated in hospital that recovered
aaV1 = str2num(answer{12}) ; % ratio of Infected treated at home that recovered
ahV2 = str2num(answer{13}) ; % ratio of Infected treated in hospital that recovered
aaV2 = str2num(answer{14}) ; % ratio of Infected treated at home that recovered
capacity0 = str2num(answer{15}) ; % hospital capacity (number of bed) 

% Determination of the others variables
t = 0 : T ; % vector of time (period) 
BETAV1 = averageContactV1 * probaTransmissionV1 ; % number of transmission per person infected per day if all contact are with Susceptible
BETAV2 = averageContactV2 * probaTransmissionV2 ;
S0 = N - I0V1 - R0 ; % number of person susceptible of being infected at time 0

% Distribution of Infectiousness relative to symptom onset = Gamma distribution : (a,b) = (2,1.5)
infectiousnessDistribition = makedist('Gamma','a',2,'b',1.5);
xInfectiousness = 0:1:10; 
yInfectiousnessV1 = pdf(infectiousnessDistribition,xInfectiousness);
yInfectiousnessV1 = (BETAV1./max(yInfectiousnessV1)).*yInfectiousnessV1 ; % choose the maximum value of the distribution as BETA
yInfectiousnessV2 = pdf(infectiousnessDistribition,xInfectiousness);
yInfectiousnessV2 = (BETAV2./max(yInfectiousnessV2)).*yInfectiousnessV2 ;

% Distribution of Incubation period over time = Lognormal distribution : (mu,s) = (1.63,0.5)
incubationDistribution = makedist('Lognormal','mu',1.63,'sigma',0.5) ;
xIncubation = 0:1:19 ; 
yIncubation = pdf(incubationDistribution,xIncubation);

% Distribution of Infectiousness after Infectious contact ( product of the Incubation period and the distribution of infectiousness)
betaV1 = zeros(numel(yIncubation)+numel(yInfectiousnessV1)-3,1) ; 
for i = 1 : numel(yIncubation)    
    mini_beta = min([2 , i-1]) ;    
    betaV1(i-mini_beta:i+8) = betaV1(i-mini_beta:i+8) + (yInfectiousnessV1(3-mini_beta:end).*yIncubation(i))';  
end
betaV2 = zeros(numel(yIncubation)+numel(yInfectiousnessV2)-3,1) ; 
for i = 1 : numel(yIncubation)    
    mini_beta = min([2 , i-1]) ;    
    betaV2(i-mini_beta:i+8) = betaV2(i-mini_beta:i+8) + (yInfectiousnessV2(3-mini_beta:end).*yIncubation(i))';  
end    

% Distribution of Recovered period over time = Distribution of Incubation + 8
xRecover = 0:1:19+8 ;
yRecover(1:8) = zeros(8,1) ; 
yRecover(9:20+8) = yIncubation ; 
yRecover = yRecover' ; 

%% ORDINARY DIFFERENTIAL EQUATIONS 

% Initial conditions
newIV1(1) = I0V1 ; % newly infected with variant 1
newIV2(1) = I0V2 ; % newly infected with variant 2
newI_hospitalV1(1,1) = h0 * newIV1(1,1) ; % number of newly infected admitted at hospital 
newI_homeV1(1,1) = (1-h0) * newIV1(1,1) ; % numner of newly infected cared at home 
newI_hospitalV2(1,1) = h0 * newIV2(1,1) ; % number of newly infected admitted at hospital 
newI_homeV2(1,1) = (1-h0) * newIV2(1,1) ; % numner of newly infected cared at home 
listV1(1,1) = 0 ; 
added_listV1(1,1) = 0 ;
listV2(1,1) = 0 ; 
added_listV2(1,1) = 0 ;
list(1,1) = 0 ; 
waiting(1,1) = 0 ; % waiting time in the queue
ratio(1,1) = 0 ; % ratio of infected with variant 1
ratio(i,2) = 0 ; % ratio of infected with variant 2
added_hospital(1,1) = 0 ; 
added_hospitalV1(1,1) = 0 ; 
added_hospitalV2(1,1) = 0 ; 
old_listV1(1,1) = 0 ; 
old_listV2(1,1) = 0 ; 

dS(1) = 0 ; % variation of Succeptible
dIV1(1) = I0V1 ; % variation of Infected (variant 1)
dIV2(1) = I0V2 ; % variation of Infected (variant 2)
dRV1(1) = 0 ; % variation of Recovered 

S(1) = S0 ; % Succeptible
IV1(1) = I0V1 ; % Infected
IV2(1) = I0V2 ; % Infected
I_hospitalV1(1) = newI_hospitalV1(1,1) ; % Infected at hospital
I_hospitalV2(1) = newI_hospitalV2(1,1) ; % Infected at hospital
I_homeV1(1) = newI_homeV1(1,1) ; % Infected at home
I_homeV2(1) = newI_homeV2(1,1) ; % Infected at home
R(1) = R0 ; % Recovered

h(1) = h0 ; % ratio of Infected that are treated in hospital
capacity(1) = capacity0 ; % number of available beds

% Iterations
for i = 2 : numel(t) 
    
   mini_beta = min(i-1,numel(betaV1)) ; 
   mini_yRecover = min(i-1,numel(yRecover));
   newIflipV1 = flip(newIV1) ; 
   newIflipV2 = flip(newIV2) ; 
   newI_hospitalflipV1 = flip(newI_hospitalV1) ;
   newI_hospitalflipV2 = flip(newI_hospitalV2) ;
   newI_homeflipV1 = flip(newI_homeV1) ;
   newI_homeflipV2 = flip(newI_homeV2) ;
   
   % newly recovered
   newR_hospitalV1(i,1) = sum(yRecover(1:mini_yRecover) .* newI_hospitalflipV1(1:mini_yRecover)) ; % number of newly recover from hospital
   newR_homeV1(i,1) = sum(yRecover(1:mini_yRecover) .* newI_homeflipV1(1:mini_yRecover)) ; % number of newly recover from home
   newRV1(i,1) = newR_hospitalV1(i,1) + newR_homeV1(i,1)  ; % number of newly recover
   newR_hospital_deathV1(i,1) = newR_hospitalV1(i,1) * (1-ahV1) ; % number of newly recover from hospital that died 
   newR_hospital_recoveredV1(i,1) = newR_hospitalV1(i,1) * ahV1 ; % number of newly recover from hospital that recovered
   newR_home_deathV1(i,1) = newR_homeV1(i,1) * (1-aaV1) ; % number of newly recover from home that died
   newR_home_recoveredV1(i,1) = newR_homeV1(i,1) * aaV1 ; % number of newly recover from home that recovered
   
   newR_hospitalV2(i,1) = sum(yRecover(1:mini_yRecover) .* newI_hospitalflipV2(1:mini_yRecover)) ; % number of newly recover from hospital
   newR_homeV2(i,1) = sum(yRecover(1:mini_yRecover) .* newI_homeflipV2(1:mini_yRecover)) ; % number of newly recover from home
   newRV2(i,1) = newR_hospitalV2(i,1) + newR_homeV2(i,1)  ; % number of newly recover
   newR_hospital_deathV2(i,1) = newR_hospitalV2(i,1) * (1-ahV2) ; % number of newly recover from hospital that died 
   newR_hospital_recoveredV2(i,1) = newR_hospitalV2(i,1) * ahV2 ; % number of newly recover from hospital that recovered
   newR_home_deathV2(i,1) = newR_homeV2(i,1) * (1-aaV2) ; % number of newly recover from home that died
   newR_home_recoveredV2(i,1) = newR_homeV2(i,1) * aaV2 ; % number of newly recover from home that recovered

   capacity(i,1) = capacity(i-1,1) - (-newR_hospitalV1(i,1)-newR_hospitalV2(i,1)) ;
   
   % newly infected
   newIV1(i,1) = sum(betaV1(1:mini_beta) .* newIflipV1(1:mini_beta)) .* S(i-1)/N ; % number of newly infected
   newIV2(i,1) = sum(betaV2(1:mini_beta) .* newIflipV2(1:mini_beta)) .* S(i-1)/N ; % number of newly infected
   ratio(i,1) = newIV1(i,1) / ( newIV1(i,1) + newIV2(i,1) ) ;
   ratio(i,2) = newIV2(i,1) / ( newIV1(i,1) + newIV2(i,1) ) ;
   
   % added to hospital from the list 
   added_hospital(i,1) = min(list(i-1,1),capacity(i,1)) ; % people that were on the list and that are accepted in hospital
   
   tol1 = 1e-2 ; 
   if listV1(i-1,1) < tol1 && listV2(i-1,1) < tol1
   added_hospitalV1(i,1) = 0 ; 
   added_hospitalV2(i,1) = 0 ;
   else
   added_hospitalV1(i,1) = listV1(i-1,1)/(listV1(i-1,1)+listV2(i-1,1))*added_hospital(i,1) ; 
   added_hospitalV2(i,1) = listV2(i-1,1)/(listV1(i-1,1)+listV2(i-1,1))*added_hospital(i,1) ;
   end
   
   % modification of lists 
   h(i,1) = min( (capacity(i,1)-added_hospital(i,1))/(newIV1(i,1)+newIV2(i,1)) , h0 ) ; % ratio of Infected that are treated in hospital this period
   new_listV1(i,1) = (h0-h(i,1))*(newIV1(i,1)) ;
   new_listV2(i,1) = (h0-h(i,1))*(newIV2(i,1)) ;
   old_listV1(i,1) = new_listV1(i,1) ; 
   old_listV2(i,1) = new_listV2(i,1) ; 
   
   added_listV1(i,1) = new_listV1(i,1) - added_hospitalV1(i,1) ; % number of newly infected added to the list minus those who left the list (waiting for hospital bed)
   added_listV2(i,1) = new_listV2(i,1) - added_hospitalV2(i,1) ;
   listV1(i,1) = listV1(i-1,1) + added_listV1(i,1) ; % number of infected on the list (waiting for hospital bed)
   listV2(i,1) = listV2(i-1,1) + added_listV2(i,1) ; % number of infected on the list (waiting for hospital bed)
   list(i,1) = listV1(i,1) + listV2(i,1) ; 
   newI_hospitalV1(i,1) = h(i,1) * newIV1(i,1) + added_hospitalV1(i,1) ; % number of newly infected admitted at hospital  
   newI_homeV1(i,1) = (1-h0) * newIV1(i,1) ; % number of newly infected cared at home
   newI_hospitalV2(i,1) = h(i,1) * newIV2(i,1) + added_hospitalV2(i,1) ; % number of newly infected admitted at hospital  
   newI_homeV2(i,1) = (1-h0) * newIV2(i,1) ; % number of newly infected cared at home
  
   capacity(i,1) = capacity(i,1) - (newI_hospitalV1(i,1)+newI_hospitalV2(i,1)) ; 
   
   % Derivatives  
   dS(i,1) = -(newIV1(i,1)+newIV2(i,1)) ; % variation of Susceptible person
   dIV1(i,1) = newIV1(i,1) - newRV1(i,1); % variation of Infected person
   dIV2(i,1) = newIV2(i,1) - newRV2(i,1); % variation of Infected person
   dI_hospitalV1(i,1) = newI_hospitalV1(i,1)-newR_hospitalV1(i,1) ; % variation of Infected person in hospital
   dI_hospitalV2(i,1) = newI_hospitalV2(i,1)-newR_hospitalV2(i,1) ; % variation of Infected person in hospital
   dI_homeV1(i,1) = newI_homeV1(i,1)-newR_homeV1(i,1) ; % variation of Infected person at home
   dI_homeV2(i,1) = newI_homeV2(i,1)-newR_homeV2(i,1) ; % variation of Infected person at home
   dRV1(i,1) = newRV1(i,1) ; % variation of Recovered person including dead
   dRV2(i,1) = newRV2(i,1) ; % variation of Recovered person including dead
   
   % update
   S(i,1) = S(i-1,1) + dS(i,1) ; % Susceptible person at time t(i)
   IV1(i,1) = IV1(i-1,1) + dIV1(i,1) ; % Infected person at time t(i
   IV2(i,1) = IV2(i-1,1) + dIV2(i,1) ; % Infected person at time t(i)
   I_hospitalV1(i,1) = I_hospitalV1(i-1,1) + dI_hospitalV1(i,1) ;
   I_hospitalV2(i,1) = I_hospitalV2(i-1,1) + dI_hospitalV2(i,1) ;
   I_homeV1(i,1) = I_homeV1(i-1,1) + dI_homeV1(i,1) ;
   I_homeV2(i,1) = I_homeV2(i-1,1) + dI_homeV2(i,1) ;
   R(i,1) = R(i-1,1) + dRV1(i,1) +dRV2(i,1) ; % Recovered person at time t(i) including dead 
   
   % waiting time
   hospital_from_list(1) = 0;
   infected_to_list(1) = 0 ; 
   hospital_from_list(i) = sum(added_hospital(:));
   tol2 = 1e-3 ; 
   for j = 2 : i
     infected_to_list(j) = sum(new_listV1(1:j))+sum(new_listV2(1:j)) ;  
     if hospital_from_list(i) >= infected_to_list(j) - tol2
         waiting(i) = 0 ; 
     elseif hospital_from_list(i) < infected_to_list(j) && hospital_from_list(i) > infected_to_list(j-1) 
         waiting(i) = i-(j-1) ;
     end
   end 
   
end
waiting = waiting' ; 
%% OUTPUT 

% Figure 


figure 
plot(-2:1:8,yInfectiousnessV1,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day] (related to the onset of symptoms)','FontSize',20,'fontWeight','bold')
ylabel('Distribution of Infectiousness : \beta','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(xIncubation,yIncubation,'Linewidth',2.5)   
hold on 
plot(xRecover,yRecover,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('Incubation/Recover period distribution','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Incubation','Recover')

figure 
plot(0:27,betaV1,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day] (related to the infectious contact)','FontSize',20,'fontWeight','bold')
ylabel('Distribution of Infectiousness : \beta','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,S,'Linewidth',2.5) 
hold on
plot(t,IV1+IV2,'Linewidth',2.5) 
plot(t,R,'Linewidth',2.5)
legend('Susceptible' , 'Infected' , 'recovered')
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('Number of people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,100*S/N,'Linewidth',2.5) 
hold on
plot(t,100*(IV1+IV2)/N,'Linewidth',2.5) 
plot(t,100*R/N,'Linewidth',2.5)
legend('Susceptible' , 'Infected' , 'recovered')
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('Number of people [% of the total population]','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,ratio(:,1),'Linewidth',2.5)
hold on 
plot(t,ratio(:,2),'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('proportion of Variant in the new infections','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Variant 1' , 'Variant 2')

figure 
plot(t,newIV1,'Linewidth',2.5)
hold on 
plot(t,newIV2,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('newly infected people infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Variant 1' , 'Variant 2')

figure 
plot(t,IV1,'Linewidth',2.5)
hold on 
plot(t,IV2,'Linewidth',2.5)
plot(t,IV1+IV2,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Variant 1' , 'Variant 2' , 'total')



figure 
plot(t,I_hospitalV1,'Linewidth',2.5)
hold on 
plot(t,I_homeV1,'Linewidth',2.5)
plot(t,listV1,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people (Variant 1)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital' , 'Home' , 'List')

figure 
plot(t,I_hospitalV2,'Linewidth',2.5)
hold on 
plot(t,I_homeV2,'Linewidth',2.5)
plot(t,listV2,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people (Variant 2)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital' , 'Home' , 'List')

figure 
plot(t,I_hospitalV2+I_hospitalV1,'Linewidth',2.5)
hold on 
plot(t,I_homeV2+I_homeV1,'Linewidth',2.5)
plot(t,listV2+listV1,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people (total)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital' , 'Home' , 'List')

figure 
plot(t,newR_hospital_recoveredV1,'Linewidth',2.5) 
hold on 
plot(t,newR_hospital_deathV1,'Linewidth',2.5) 
plot(t,newR_home_recoveredV1,'Linewidth',2.5) 
plot(t,newR_home_deathV1,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly recovered people (Variant 1)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,newR_hospital_recoveredV2,'Linewidth',2.5) 
hold on 
plot(t,newR_hospital_deathV2,'Linewidth',2.5) 
plot(t,newR_home_recoveredV2,'Linewidth',2.5) 
plot(t,newR_home_deathV2,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly recovered people (Variant 2)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,cumsum(newI_hospitalV1+newI_hospitalV2),'Linewidth',2.5) 
hold on 
plot(t,cumsum(newI_homeV1+newI_homeV2),'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('cumulative number of infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home')


figure 
plot(t,cumsum(newR_hospital_recoveredV1),'Linewidth',2.5) 
hold on 
plot(t,cumsum(newR_hospital_deathV1),'Linewidth',2.5) 
plot(t,cumsum(newR_home_recoveredV1),'Linewidth',2.5) 
plot(t,cumsum(newR_home_deathV1),'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('cumulative number of recovered people (Variant 1)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,cumsum(newR_hospital_recoveredV2),'Linewidth',2.5) 
hold on 
plot(t,cumsum(newR_hospital_deathV2),'Linewidth',2.5) 
plot(t,cumsum(newR_home_recoveredV2),'Linewidth',2.5) 
plot(t,cumsum(newR_home_deathV2),'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('cumulative number of recovered people (Variant 2)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,newI_hospitalV1,'Linewidth',2.5) 
hold on 
plot(t,newI_homeV1,'Linewidth',2.5)  
plot(t,newIV1,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly infected people (Variant 1)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home','Total')

figure 
plot(t,newI_hospitalV2,'Linewidth',2.5) 
hold on 
plot(t,newI_homeV2,'Linewidth',2.5)  
plot(t,newIV2,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly infected people (Variant 2)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home','Total')

figure 
plot(t,h,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('ratio of new Infected that are treated in hospital : h','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,waiting,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('the time that people entering the hospital from the list have waited [day]','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

t_adj = t-(waiting') ; 
figure 
plot(t_adj,waiting,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('waiting time before entering the hospital [day]','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')


figure 
plot(t,new_listV1,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('new people on the list (Variant 1)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,new_listV2,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('new people on the list (Variant 2)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,added_hospitalV1 + added_hospitalV2 ,'Linewidth',2.5) 
hold on 
plot(t,added_hospitalV1,'Linewidth',2.5)
plot(t,added_hospitalV2,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('new people accepted at hospital from list','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('total' , 'from list : Variant 1' , 'from list : Variant 2')
%{
figure 
plot(t,newI,'Linewidth',2.5) 
hold on
plot(t,[36 45 0 37 141 100 0 400 622 0 0 2955 1159 2144 1806 2162 0 2447 4964 0 6368 4749 9630 0 7939 7516],'Linewidth',2.5)
%plot(t,[16 21 0 67 83 117 0 201 270 0 405 556 674 0 1399 2444 0 5619 6516 0 10432 10433 14634 0 17330 21066],'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Model','Spain')

%plot(t,[36 45 0 37 141 100 0 400 622 0 0 2955 1159 2144 1806 2162 0 2447 4964 0 6368 4749 9630 0 7939 7516],'Linewidth',2.5)
%}