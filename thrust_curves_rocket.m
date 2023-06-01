function [eng_w,liquid_consumpt_rocket,oxidizer_consumpt_rocket,thrust] = thrust_curves_rocket(eng_name,p)
if eng_name == 0 % Nerv
    atm_curve = [0 1; 60 13.88];
    eng_w = 3000; % en kg
    eng_cons = 1.53*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons;
    oxidizer_consumpt_rocket = 0;
elseif eng_name == 1 % Poodle
    atm_curve = [0 1; 250 64.29];
    eng_w = 1750; % en kg
    eng_cons = 14.568*5; % en kg/s
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 2 % Dart
    atm_curve = [0 1; 180 153.53];
    eng_w = 1000; % en kg
    eng_cons = 10.797*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 3 % Swivel
    atm_curve = [0 1; 215 167.97];
    eng_w = 1500; % en kg
    eng_cons = 13.703*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 4 % Reliant
    atm_curve = [0 1; 240 205.16]; % En kN
    eng_w = 1250; % en kg
    eng_cons = 15.789*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 5 % Vector
    atm_curve = [0 1; 1000 936.51];
    eng_w = 4000; % en kg
    eng_cons = 64.745*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 6 % Skipper
    atm_curve = [0 1; 650 568.75];
    eng_w = 3000; % en kg
    eng_cons = 41.426*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
elseif eng_name == 7 % Mainsail
    atm_curve = [0 1; 1500 1379.03];
    eng_w = 6000; % en kg
    eng_cons = 98.683*5; % en kg/s 
    liquid_consumpt_rocket = eng_cons*(0.9/2);
    oxidizer_consumpt_rocket = eng_cons*(1.1/2);
end

xx = min(atm_curve(1,:)):0.01:max(atm_curve(1,:));
atm_curve_F(2,:) = ((xx-atm_curve(1,1))/(atm_curve(1,2)-atm_curve(1,1)))*(atm_curve(2,2)-atm_curve(2,1))+atm_curve(2,1);
atm_curve_F(1,:) = xx;

[~, index] = min(abs(atm_curve_F(1,:) - round(p,2)));
thrust = atm_curve_F(2,index);

end