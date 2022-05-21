clc
close all 
clear all 

%% VARIABLES 

% Input (LAUNCH "Run" TO ENTER INPUT INTO WINDOW) 
prompt = {'T : number of periods', ...
    'N : size of the population', ...
    'Average number of contact per person per period at peak infectiousness', ...
    'probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness', ...
    'number of persons infected at time 0',...
    'number of person who have recovered at time 0 (i.e: person who are resistant (recovered, vaccinated, ...)',...
    'h0 : ratio of Infected that are treated in hospital if capacity is enough',...
    'ah : ratio of Infected treated in hospital that recovered' ,...
    'aa : ratio of Infected treated at home that recovered',...
    'hospital capacity (number of bed)'};
dlgtitle = 'Inputs';
dims = [1 120];
definput = {'300','47e6','10',...
'0.20','16','0','0.04','0.96','0.99','1.6e5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);  
 
T = str2num(answer{1}) ;  % number of periods
N = str2num(answer{2}) ; % size of the population
averageContact = str2num(answer{3}) ; % Average number of contact per person per period at peak infectiousness
probaTransmission = str2num(answer{4}) ; % probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness
I0 = str2num(answer{5}) ; % number of persons infected at time 0
R0 = str2num(answer{6}) ; % number of person who have recovered at time 0 (i.e: person who are resistant (recovered, vaccinated, ...) 
h0 = str2num(answer{7}) ; % ratio of Infected that are treated in hospital if capacity is enough
ah = str2num(answer{8}) ; % ratio of Infected treated in hospital that recovered
aa = str2num(answer{9}) ; % ratio of Infected treated at home that recovered
capacity0 = str2num(answer{10}) ; % hospital capacity (number of bed) 

% Determination of the others variables
t = 0 : T ; % vector of time (period) 
BETA = averageContact * probaTransmission ; % number of transmission per person infected per day if all contact are with Susceptible
S0 = N - I0 - R0 ; % number of person susceptible of being infected at time 0

% Distribution of Infectiousness relative to symptom onset = Gamma distribution : (a,b) = (2,1.5)
infectiousnessDistribition = makedist('Gamma','a',2,'b',1.5);
xInfectiousness = 0:1:10; 
yInfectiousness = pdf(infectiousnessDistribition,xInfectiousness);
yInfectiousness = (BETA./max(yInfectiousness)).*yInfectiousness ; % choose the maximum value of the distribution as BETA

% Distribution of Incubation period over time = Lognormal distribution : (mu,s) = (1.63,0.5)
incubationDistribution = makedist('Lognormal','mu',1.63,'sigma',0.5) ;
xIncubation = 0:1:19 ; 
yIncubation = pdf(incubationDistribution,xIncubation);

% Distribution of Infectiousness after Infectious contact ( product of the Incubation period and the distribution of infectiousness)
beta = zeros(numel(yIncubation)+numel(yInfectiousness)-3,1) ; 
for i = 1 : numel(yIncubation)    
    mini_beta = min([2 , i-1]) ;    
    beta(i-mini_beta:i+8) = beta(i-mini_beta:i+8) + (yInfectiousness(3-mini_beta:end).*yIncubation(i))';  
end
    
% Distribution of Recovered period over time = Distribution of Incubation + 8
xRecover = 0:1:19+8 ;
yRecover(1:8) = zeros(8,1) ; 
yRecover(9:20+8) = yIncubation ; 
yRecover = yRecover' ; 

%% ORDINARY DIFFERENTIAL EQUATIONS 

% Initial conditions
newI(1) = I0 ; % newly infected
newI_hospital(1,1) = h0 * newI(1,1) ; % number of newly infected admitted at hospital 
newI_home(1,1) = (1-h0) * newI(1,1) ; % numner of newly infected cared at home 
list(1,1) = 0 ; 
added_list(1,1) = 0 ; 
waiting(1,1) = 0 ; % waiting time in the queue

dS(1) = 0 ; % variation of Succeptible
dI(1) = I0 ; % variation of Infected
dR(1) = 0 ; % variation of Recovered

S(1) = S0 ; % Succeptible
I(1) = I0 ; % Infected
I_hospital(1) = newI_hospital(1,1) ; % Infected at hospital
I_home(1) = newI_home(1,1) ; % Infected at home 
R(1) = R0 ; % Recovered

h(1) = h0 ; % ratio of Infected that are treated in hospital
capacity(1) = capacity0 ; % number of available beds

% Iterations
for i = 2 : numel(t) 
    
   mini_beta = min(i-1,numel(beta)) ; 
   mini_yRecover = min(i-1,numel(yRecover));
   newIflip = flip(newI) ; 
   newI_hospitalflip = flip(newI_hospital) ;
   newI_homeflip = flip(newI_home) ;
   
   % newly recovered
   newR_hospital(i,1) = sum(yRecover(1:mini_yRecover) .* newI_hospitalflip(1:mini_yRecover)) ; % number of newly recover from hospital
   newR_home(i,1) = sum(yRecover(1:mini_yRecover) .* newI_homeflip(1:mini_yRecover)) ; % number of newly recover from home
   newR(i,1) = newR_hospital(i,1) + newR_home(i,1)  ; % number of newly recover
   newR_hospital_death(i,1) = newR_hospital(i,1) * (1-ah) ; % number of newly recover from hospital that died 
   newR_hospital_recovered(i,1) = newR_hospital(i,1) * ah ; % number of newly recover from hospital that recovered
   newR_home_death(i,1) = newR_home(i,1) * (1-aa) ; % number of newly recover from home that died
   newR_home_recovered(i,1) = newR_home(i,1) * aa ; % number of newly recover from home that recovered
   
   capacity(i,1) = capacity(i-1,1) - (-newR_hospital(i,1)) ;
   
   % newly infected
   newI(i,1) = sum(beta(1:mini_beta) .* newIflip(1:mini_beta)) .* S(i-1)/N ; % number of newly infected
   
   added_hospital(i,1) = min(list(i-1,1),capacity(i,1)) ; % people that were on the list and that are accepted in hospital
   
   h(i,1) = min( (capacity(i,1)-added_hospital(i,1))/newI(i,1) , h0 ) ; % ratio of Infected that are treated in hospital this period
   new_list(i,1) = (h0-h(i,1))*newI(i,1) ; 
   added_list(i,1) = new_list(i,1) - added_hospital(i,1) ; % number of newly infected added to the list minus those who left the list (waiting for hospital bed)
   list(i,1) = list(i-1,1) + added_list(i,1) ; % number of infected on the list (waiting for hospital bed)
   newI_hospital(i,1) = h(i,1) * newI(i,1) + added_hospital(i,1) ; % number of newly infected admitted at hospital  
   newI_home(i,1) = (1-h0) * newI(i,1) ; % number of newly infected cared at home 
  
   capacity(i,1) = capacity(i,1) - (newI_hospital(i,1)) ; 
   
   % Derivatives  
   dS(i,1) = -newI(i,1) ; % variation of Susceptible person
   dI(i,1) = newI(i,1) - newR(i,1); % variation of Infected person
   dI_hospital(i,1) = newI_hospital(i,1)-newR_hospital(i,1) ; % variation of Infected person in hospital
   dI_home(i,1) = newI_home(i,1)-newR_home(i,1) ; % variation of Infected person at home
   dR(i,1) = newR(i,1) ; % variation of Recovered person including dead
   
   % update
   S(i,1) = S(i-1,1) + dS(i,1) ; % Susceptible person at time t(i)
   I(i,1) = I(i-1,1) + dI(i,1) ; % Infected person at time t(i)
   I_hospital(i,1) = I_hospital(i-1,1) + dI_hospital(i,1) ; 
   I_home(i,1) = I_home(i-1,1) + dI_home(i,1) ; 
   R(i,1) = R(i-1,1) + dR(i,1) ; % Recovered person at time t(i) including dead 
   
   % waiting time
   hospital_from_list(i,1) = sum(added_hospital(:));
   tol = 1e-4 ; 
   for j = 2 : i
     infected_to_list(j,1) = sum(new_list(1:j)) ; 
     if hospital_from_list(i,1) >= infected_to_list(j,1)-tol
         waiting(i,1) = 0 ; 
     elseif hospital_from_list(i,1) < infected_to_list(j,1) && hospital_from_list(i,1) > infected_to_list(j-1,1) 
         waiting(i,1) = i-(j-1) ; 
     end
   end 
   
end

%% OUTPUT 

% Figure 
figure 
plot(t,S,'Linewidth',2.5) 
hold on
plot(t,I,'Linewidth',2.5) 
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
plot(t,100*I/N,'Linewidth',2.5) 
plot(t,100*R/N,'Linewidth',2.5)
legend('Susceptible' , 'Infected' , 'recovered')
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('Number of people [% of the total population]','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(t,I,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(-2:1:8,yInfectiousness,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day] (related to the onset of symptoms)','FontSize',20,'fontWeight','bold')
ylabel('Distribution of Infectiousness : \beta','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(xIncubation,yIncubation,'Linewidth',2.5)   
hold on 
plot(xRecover,yRecover,'Linewidth',2.5)
ylabel('Incubation/Recover period distribution','FontSize',20,'fontWeight','bold')
yyaxis right 
plot(0:27,beta,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day] (related to the infectious contact)','FontSize',20,'fontWeight','bold')
ylabel('Infectiousness period distribution : \beta ','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Incubation period distribution','Recover period distribution','Infectiousness period distribution')

figure 
plot(xRecover,yRecover,'Linewidth',2.5)
ylabel('Recover period distribution : \gamma','FontSize',20,'fontWeight','bold')
grid on 
grid minor
xlabel('time [day] (related to the infectious contact)','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')

figure 
plot(0:27,beta,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day] (related to the infectious contact)','FontSize',20,'fontWeight','bold')
ylabel('Infectiousness period distribution : \beta ','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')


figure 
plot(t,newR_hospital_recovered,'Linewidth',2.5) 
hold on 
plot(t,newR_hospital_death,'Linewidth',2.5) 
plot(t,newR_home_recovered,'Linewidth',2.5) 
plot(t,newR_home_death,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly recovered people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,cumsum(newI_hospital),'Linewidth',2.5) 
hold on 
plot(t,cumsum(newI_home),'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('cumulative number of infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home')


figure 
plot(t,cumsum(newR_hospital_recovered),'Linewidth',2.5) 
hold on 
plot(t,cumsum(newR_hospital_death),'Linewidth',2.5) 
plot(t,cumsum(newR_home_recovered),'Linewidth',2.5) 
plot(t,cumsum(newR_home_death),'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('cumulative number of recovered people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital recovered', 'Hospital death', 'Home recovered', 'Home death')

figure 
plot(t,newI_hospital,'Linewidth',2.5) 
hold on 
plot(t,newI_home,'Linewidth',2.5)  
plot(t,newI,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of newly infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home','Total')

figure 
plot(t,I_hospital,'Linewidth',2.5) 
hold on 
plot(t,I_home,'Linewidth',2.5)
plot(t,list,'Linewidth',2.5)
yline(capacity0,'Linewidth',2.5)
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('number of infected people','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
legend('Hospital','Home','waiting list','maximum capacity at hospital')

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
plot(t,new_list,'Linewidth',2.5) 
grid on 
grid minor
xlabel('time [day]','FontSize',20,'fontWeight','bold')
ylabel('new people on the list','FontSize',20,'fontWeight','bold')
set(gca,'FontSize',20,'fontWeight','bold')
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