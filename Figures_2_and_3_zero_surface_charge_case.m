close all 
clear all

%note at 20C
% only works for 1, -1 symetric electrolytes. 

if(true) %1mM case in paper
    pH = 6;
    pOH = 8;
    pROH2 = 10;
    pRO = 8.83;
    pKa = 16.84;
    pKb = 15.99;
    pure_liquids = false;
    
    %graph limits
    xlimits = [-0.05 0.05]; %um
    ylimits_phi = [-160 160]; %mV
    ylimits_E = [-2e7 2e7]; %V/m
    ylimits_conc = [1e-12 1e-1]; %M
end

if(false) %pure case in paper
    pH = 7;
    pOH = 7;
    pROH2 = 9.415;
    pRO = 9.415;
    pKa = 16.84;
    pKb = 15.99;
    pure_liquids = true;
    
    %graph limits
    xlimits = [-30 30]; %um
    ylimits_phi = [-40 40]; %mV
    ylimits_E = [-5e3 5e3]; %V/m
    ylimits_conc = [1e-10 1e-6]; %M
end

if ~pure_liquids
    water_ionic_strength  = 1e-3;   %in M
    octanol_ionic_strength = 1e-3; %in M
else
    water_ionic_strength  = (10^-pH + 10^-pOH)/2;
    octanol_ionic_strength = (10^-pROH2 + 10^-pRO)/2;
end


%check inputs
if ~(pH + pOH == 14)
    disp("warning: bulk water not at equilibrium Kw \n aborting run")
    return
elseif abs(pROH2+pRO - pKa - pKb + 14) > 1e-12
    disp("warning: bulk alcohol not at equilibrium with estimated K_ROH \n aborting run")
    return
elseif (water_ionic_strength < (10^-pH + 10^-pOH)/2)
    disp("warning: water ionic strength less than autoionzation at given pH \n aborting run")
    return  
elseif (octanol_ionic_strength < (10^-pROH2 + 10^-pRO)/2)
    disp("warning: alcohol ionic strength less than autoionzation at given pROH2 \n aborting run")
    return
end

%constants
N_A = 6.0221409 *10^23;
Vt = 0.025261709807677;                                                    %thermal voltage at 20C
In_e = 1.60217662 * 10^-19;                                                %charge on an electron in C 
epsilon0 = 8.854*10^-12;                                                   %in F/m   
epsilonWater =  80.2;                                                      %dielectric for water at 20C 
epsilonOctanol = 10.3;                                                     %dielectric for 1-octanol at 20C

%Calculate simple things

H_bulk = 10^-pH;
OH_bulk = 10^-pOH;
ROH2_bulk = 10^-pROH2;
RO_bulk = 10^-pRO;

water_sum_n = 2 * water_ionic_strength * N_A * 1000;    %in ions per m^3
ROH_sum_n = 2 * octanol_ionic_strength * N_A * 1000;  %in ions per m^3    


K_H2O = sqrt(In_e * water_sum_n / (epsilon0*epsilonWater * Vt));          %Hunter eq 2.3.9
K_ROH = sqrt(In_e * ROH_sum_n  / (epsilon0*epsilonOctanol * Vt));          %Hunter eq 2.3.9


%Calc phi
phi_ROH = - Vt/2 * log(10) * (pKa + pOH + pROH2 - pKb - pH - pRO); % see derivation in paper
x = linspace(0,max(abs(xlimits))*1e-6,1000);

%calc phi surface with no charge
syms phiS real
dphi_H2O = 4 * (water_ionic_strength * N_A * 1000) * In_e / (K_H2O * epsilonWater * epsilon0) * sinh(phiS/(2*Vt));
dphi_ROH = 4 * (octanol_ionic_strength * N_A * 1000) * In_e / (K_ROH * epsilonOctanol * epsilon0) * sinh((phiS - phi_ROH)/(2*Vt));

eq1 = dphi_H2O == -dphi_ROH;

phi_surface = double(solve(eq1,phiS));
    
phi_H2O = calc_phi(x,K_H2O,phi_surface,0,Vt);
phi_alcohol = calc_phi(x,K_ROH,phi_surface,phi_ROH,Vt);

ion_H = H_bulk * exp(-phi_H2O/Vt); 
ion_OH = OH_bulk * exp(phi_H2O/Vt);
ion_ROH2 = ROH2_bulk * exp((-phi_alcohol + phi_ROH)/Vt); 
ion_RO = RO_bulk * exp((phi_alcohol - phi_ROH)/Vt); 



if (true)
    figure(1)

    if(~pure_liquids)
        ion_all_pos_H2O = water_ionic_strength * exp(-phi_H2O/Vt); 
        ion_all_neg_H2O = water_ionic_strength * exp(phi_H2O/Vt);
        ion_all_pos_ROH = octanol_ionic_strength * exp((-phi_alcohol + phi_ROH)/Vt); 
        ion_all_neg_ROH = octanol_ionic_strength * exp((phi_alcohol - phi_ROH)/Vt); 

        ion_all_pos = [ion_all_pos_H2O, NaN, ion_all_pos_ROH]; 
        ion_all_neg = [ion_all_neg_H2O, NaN, ion_all_neg_ROH]; 
        graphx_conc = [-x, 0, x];
        
        semilogy(-graphx_conc*1e6,ion_all_pos,'LineWidth',2)
        hold on
        semilogy(-graphx_conc*1e6,ion_all_neg,'LineWidth',2)
    end
    
    semilogy(-x*1e6,ion_H,'LineWidth',2)
    hold on
    semilogy(-x*1e6,ion_OH,'LineWidth',2)
    semilogy(x*1e6,ion_ROH2,'LineWidth',2)
    semilogy(x*1e6,ion_RO,'LineWidth',2)  
    
    
    
end

figure(100)
graphphi = [phi_H2O(end:-1:2) phi_alcohol]*1000; 
graphx = [-x(end:-1:2) x];
graphE = -diff(graphphi)/(graphx(2)-graphx(1))/1000;
graphxE = graphx(1:end-1) + (graphx(2)- graphx(1))/2;



yyaxis left
plot(graphx*1e6,graphphi,'LineWidth',2,'DisplayName',"\phi(x) (mV)")
yyaxis right
plot(graphxE*1e6,graphE,'LineWidth',2,'DisplayName',"E(x) (V/m)")
hold on

yyaxis right
xline(0,'k--', 'LineWidth',2,'DisplayName','');
legend(["\phi(x)","E(x)",],'Location','northeast')
yyaxis left

xlim(xlimits);
xlabel("Distance (\mum)");

yyaxis left
ylabel("\phi(x) (mV)");
ylim(ylimits_phi)

yyaxis right
ylabel("E(x) (V/m)");
ylim(ylimits_E)




figure(1)
hold on
xline(0,'k--', 'LineWidth',2,'DisplayName','');

if ~pure_liquids
legend(["sum all +1 ions","sum all -1 ions", "H_3O^+","OH^-","CH_3(CH_2)_7OH_2^+","CH_3(CH_2)_7O^-"],'Location','northeast')
else
legend(["H_3O^+","OH^-","CH_3(CH_2)_7OH_2^+","CH_3(CH_2)_7O^-"],'Location','northeast')
end

xlim(xlimits);
ylim(ylimits_conc);
xlabel("Distance (\mum)");
ylabel("Concentration (M)");


function phi = calc_phi(x,K,phi_surface,phi_bulk,Vt)
     phi0 = phi_surface - phi_bulk;
     phi = 4 * Vt * atanh(tanh(phi0./(4*Vt)).*exp(-K*x)); %from Hunter textbook, valid only for +1, -1 electrolytes 
     phi = phi + phi_bulk;
end