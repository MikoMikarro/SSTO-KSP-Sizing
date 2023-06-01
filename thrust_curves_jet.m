function [eng_w,eng_cons,thrust_0,vel_mult,atm_mult] = thrust_curves_jet(eng_name,M,p)

if eng_name == 0 % Panther
    vel_curve = [0 0.35 1 1.75 2 2.5; 1 0.932 1.13 1.5 1.38 0]; %velocidad en mach
    atm_curve = [0 0.072 0.17 0.34 1; 0 0.08 0.21 0.39 0.969697]; %presión en atm
    eng_w = 1200; % en kg
    eng_cons = 0.66*5; % en kg/s 
    thrust_0 = 85e3; % en N
elseif eng_name == 1 % C.R.A.P.I.E.R.
    vel_curve = [0 0.2 0.7 1.4 3.75 4.5 5.5 6; 1 0.98 1.8 4 8.5 7.3 3 0];
    atm_curve = [0 0.018 0.08 0.35 1; 0 0.09 0.3 0.5 1.055097];
    eng_w = 2000; % en kg
    eng_cons = 0.67*5; % en kg/s
    thrust_0 = 105e3;
elseif eng_name == 2 % Whiplash
    vel_curve = [0 0.2 0.72 1.36 2.15 3 4.5 5.5; 1 0.98 1.716 3.2 4.9 5.8 3 0];
    atm_curve = [0 0.045 0.16 0.5 1; 0 0.166 0.5 0.6 1.013946];
    eng_w = 1800; % en kg
    eng_cons = 0.66*5; % en kg/s
    thrust_0 = 130e3;
end

xx = min(vel_curve(1,:)):0.01:max(vel_curve(1,:));
vel_curve_F(2,:) = spline(vel_curve(1,:),vel_curve(2,:),xx);
vel_curve_F(1,:) = xx;

xx = min(atm_curve(1,:)):0.01:max(atm_curve(1,:));
atm_curve_F(2,:) = spline(atm_curve(1,:),atm_curve(2,:),xx);
atm_curve_F(1,:) = xx;

[~, index] = min(abs(atm_curve_F(1,:) - round(p,2)));
atm_mult = atm_curve_F(2,index);

[~, index] = min(abs(vel_curve_F(1,:) - round(M,2)));
vel_mult = vel_curve_F(2,index);
    
end