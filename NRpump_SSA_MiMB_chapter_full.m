%Daniel Charlebois - MATLAB R2016b - June 2017
%"Skeleton" code for the stochastic simulation of the 
%NRpump gene circuit for the Charlebois and Diao et al. book chapter  
%in Methods in Molecular Biology. 
%*Note that this "full" version of the code explicitly models each reaction 
%in the NRpump system individually (the other code is an approximation, 
%where the "source" and "sink" reactions for each variable are grouped 
%together, which simplified the code/reduced the number of reactions for 
%the MiMB book chapter and still captured the steady-state behavior 
%of the systems considered in the Diao et al., ACS Synth. Biol., 2016 
%study reasonably well). 

%rates
ax = 43;      %maximal TetR production rate 
az = 50;      %maximal PDR5::GFP production rate 
lz = 10;      %leaky PDR5::GFP production rate 
bprime = 10;  %doxycycline-TetR binding rate 
d = 0.12;     %TetR and PDR5::GFP dilution rate   
f = 1.2;      %doxycycline diffusion rate 
n = 4;        %Hill coefficient for repression function
theta = 0.44; %TetR-promoter dissociation constant   
k = 200;      %rate of PDR5 mediated efflux 
K = 50;       %half-maximal pump activation parameter   
h = 3.5;      %Hill coefficient for pump term
c = 10;       %doxycycline influx rate 
dox = 0:15;   %extracellular doxycyline concentration 
t_end = 200;  %simulation end time

dox_cnt = 0;
for C = dox*c
    dox_cnt = dox_cnt+1;

    %species
    x = 0; %unbound/active TetR protein
    y = 0; %intracellular doxycycline
    z = 0; %efflux pump-reporter protein 

    %Gilesspie SSA loop
    t = 0; %initialize simulation time
    while t < t_end
        %calculate reaction propensities
        a = [ax bprime*x*y x*d C ... 
             y*f (k*y^h*z)/(K^h+y^h) ...
             lz az*(theta^n/(theta^n+x^n)) z*d];
        a0 = sum(a); %combined reaction hazard

        %calculate time to next event
        r1 = rand;
        while r1 == 0
            r1 = rand;
        end
        t_next = ((1/a0)*(log(1/r1)));
        %update time
        t = t + t_next;
    
        %determine next reaction
        i = 1; mu = 0; amu = 0; r2 = rand;
        while amu < r2*a0
            mu = mu + 1;
            amu = amu + a(i); 
            i = i + 1;
        end
    
        %reactions
        if mu == 1 
            x = x + 1;
        elseif mu == 2 
            x = x - 1; y = y - 1;
        elseif mu == 3
            x = x - 1;
        elseif mu == 4 
            y = y + 1;
        elseif mu == 5
            y = y - 1;
        elseif mu == 6
            y = y - 1;
        elseif mu == 7
            z = z + 1;
        elseif mu == 8
            z = z + 1;
        elseif mu == 9
            z = z - 1;
        end
        
    end
    
    %save SSA data
    SSA_data(dox_cnt,1:3) = [x y z];
    
end

%figure
plot(dox,SSA_data(:,3),'k.--')
xlabel('doxycycline concentration (ug/ml)'); ylabel('reporter-pump protein number')