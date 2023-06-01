clear all
close all
clc


jet_eng_name = 1;
rocket_eng_name = 2;
n_jet_engine = 2;
n_rocket_engine = 1;

str_mass = 787.294679999998;                        % kg
payload_mass = 7.5e3;                               % kg
jet_fuel = 750.399999999994;                        % kg
rocket_liquid_fuel = 930.91734;                     % kg
liquid_fuel_mass = jet_fuel + rocket_liquid_fuel;
oxidizer_fuel_mass = 1137.78786;                    % kg

S = 2.5;                                            % m^2

G = 6.674e-11;
M_Kerbin = 5.2915158e22;
R_Kerbin = 6e5;
mu = 3.5316000e12;

cd0 = 0.06;
k = .01;
cl0 = 0.0;
cl_alpha = 2*pi;


[eng_mass_jet,fuel_consumpt_air,air_thrust,vel_mult,atm_mult] = thrust_curves_jet(jet_eng_name,0,1);
[eng_mass_rocket,liquid_consumpt_rocket,oxidizer_consumpt_rocket,rocket_thrust] = thrust_curves_rocket(rocket_eng_name,1);


total_mass = str_mass + n_jet_engine*eng_mass_jet + n_rocket_engine*eng_mass_rocket + payload_mass + liquid_fuel_mass + oxidizer_fuel_mass;
vehicle_dry_mass = total_mass - liquid_fuel_mass - oxidizer_fuel_mass;    % kg

g = 9.81; % [m/s^2]
h_speed = 0; % [m/s]
v_speed = 0; % [m/s]

dt = .01; % [s]

t = 0;
takeoff = 1;
alpha = 0; % [rad]
vehicle_angle = 0; % [rad]
height = 0;     % [m]
advanced = 0;   % [m]
x_points = [0];
y_points = [0];
target_vehicle_angle = 30*pi/180;

archived_terminal_vel = 1;
angle_change_velocity= 0.05; % rad/s
got_inclined = 0;
figure();
hold on
index = 1;
%%
while true
    [h, t_e, p_e, rho_e, c] = Earth_to_Kerbin(height/1000);
    p_e = p_e / 101325;
    speed = sqrt(v_speed^2 + h_speed^2);
    M = speed/c;
    [eng_mass,fuel_consumpt_air,air_thrust,vel_mult,atm_mult] = thrust_curves_jet(jet_eng_name,M,p_e);
    real_th = n_jet_engine*air_thrust*vel_mult*atm_mult;
    
    if (archived_terminal_vel == 0)
        if vehicle_angle < target_vehicle_angle
            vehicle_angle = vehicle_angle + angle_change_velocity*dt;
        else
            vehicle_angle = target_vehicle_angle;
            got_inclined = 1;
        end
    end

    if speed > 0
        speed_angle = atan(v_speed/h_speed);
    else
        speed_angle = 0;
    end
    if takeoff
        perceived_alpha = 10*pi/180;
    else
        perceived_alpha = vehicle_angle - speed_angle + 5*pi/180;
    end

    % Considering alpha 0 
    Cl = cl0 + cl_alpha * perceived_alpha;
    Cd = cd0; %+ k*Cl^2;

    L = Cl * rho_e * speed^2 * S / 2;
    D = Cd * rho_e * speed^2 * S / 2;
    Tm = vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass;
    W = Tm * g;
    
    if (L >= W)
        if (takeoff == 1)
            disp("Takeoff!");
            fprintf("Takeoff speed: %f m/s\n", speed);
            fprintf("Takeoff distance: %f m\n", advanced);
            disp("Put your vehicle in 0º of pitch");
            xline(advanced)
            disp("_____");
            takeoff = 0;
        end
    end

    if liquid_fuel_mass <= 0
        xline(advanced)
        disp("Out of fuel!");
        fprintf("Max speed %f m/s\n", speed);
        break;
    end
    total_front_force = real_th*cos(vehicle_angle) - D*cos(vehicle_angle) - L*sin(vehicle_angle);
    total_vertical_force = real_th*sin(vehicle_angle) + L*cos(vehicle_angle) -D*sin(vehicle_angle) - W;
    dhV = total_front_force / (vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass);
    dvV = total_vertical_force / (vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass);

    if (dvV < 0) && (height == 0) % Aún no ha arrancad
        dvV = 0;
    end

    new_hspeed = h_speed + dhV * dt;
    new_vspeed = v_speed + dvV * dt;

    new_speed = sqrt(new_vspeed^2 + new_hspeed^2);


    if(archived_terminal_vel)
        if dhV < 0.1
            fprintf("Reached terminal speed: %f m/s\n", speed);
            fprintf("Remaining fuel: %f kg\n", liquid_fuel_mass);
            fprintf("Current height: %f m\n", height);
            archived_terminal_vel = 0;
            xline(advanced);
            disp("Put your vehicle in 30º of pitch");
            disp("_____");
        end
    end
    
    if speed_angle < target_vehicle_angle/2
        if (got_inclined && archived_terminal_vel == 0)
            fprintf("Reached terminal height: %f m/s\n", speed);
            fprintf("Remaining fuel (Liquid): %f kg\n", liquid_fuel_mass);
            fprintf("Current height: %f m\n", height);
            xline(advanced);
            disp("Put your vehicle in 30º of pitch and enable rocket boosters!");
            disp("_____");
            break
        end
    end
    h_speed = new_hspeed;
    v_speed = new_vspeed;
    advanced = advanced + h_speed*dt;
    height = height + v_speed*dt;
    liquid_fuel_mass = liquid_fuel_mass - n_jet_engine * fuel_consumpt_air * dt;

    x_points(index) = advanced;
    y_points(index) = height;
    index = index +1;

end

plot(x_points, y_points)
ylabel("Height [m]")
xlabel("Distance [m]")
title("Aerodynamic phase")
daspect([1 1 1])
hold off;

%% Parte aerodinámica DONE, empieza la orbital

target_vehicle_angle = 30 * pi/180;
h_speed = h_speed + 174; %% Adding sideral rotational velocity

while true
    [h, t_e, p_e, rho_e, c] = Earth_to_Kerbin(height/1000);
    p_e = p_e / 101325;
    
    [eng_mass_rocket,liquid_consumpt_rocket,oxidizer_consumpt_rocket,rocket_thrust] = thrust_curves_rocket(rocket_eng_name,p_e);

    
    real_th = n_rocket_engine * rocket_thrust * 1000;
    
    if (archived_terminal_vel == 0)
        if vehicle_angle < target_vehicle_angle
            vehicle_angle = vehicle_angle + angle_change_velocity*dt;
        else
            vehicle_angle = target_vehicle_angle;
        end
    end

    speed = sqrt(v_speed^2 + h_speed^2);

    if speed > 0
        speed_angle = atan(v_speed/h_speed);
    end

    perceived_alpha = vehicle_angle - speed_angle;

    Cl = cl0 + cl_alpha * perceived_alpha;
    Cd = cd0 ;
    
    air_speed =  sqrt(v_speed^2 + (h_speed-174)^2);
    L = Cl * rho_e * air_speed^2 * S / 2;
    D = Cd * rho_e * air_speed^2 * S / 2;

    W = (vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass) * g;

    r0 = height+R_Kerbin;

    total_front_force = real_th*cos(vehicle_angle) - D*cos(vehicle_angle) - L*sin(vehicle_angle);
    total_vertical_force = real_th*sin(vehicle_angle) + L*cos(vehicle_angle) -D*sin(vehicle_angle) - W;
    dhV = total_front_force / (vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass);
    dvV = total_vertical_force / (vehicle_dry_mass + liquid_fuel_mass + oxidizer_fuel_mass);
    dvV = dvV + h_speed^2/r0;


    new_hspeed = h_speed + dhV * dt;
    new_vspeed = v_speed + dvV * dt;

    new_speed = sqrt(new_vspeed^2 + new_hspeed^2);

    h_speed = new_hspeed;
    v_speed = new_vspeed;

    height = height + v_speed*dt;
    
    liquid_fuel_mass = liquid_fuel_mass - n_rocket_engine * liquid_consumpt_rocket * dt;
    oxidizer_fuel_mass = oxidizer_fuel_mass - n_rocket_engine * oxidizer_consumpt_rocket * dt;

    

    energia_cinetica = 1/2 * speed^2 - (mu/(r0));

    if speed > 0
        speed_angle = atan(v_speed/h_speed);
    end

    semieje_mayor = -mu/(2*energia_cinetica);
    momento_angular = r0*speed*cos(speed_angle);
    excentricidad = sqrt(1-((momento_angular^2)/(mu*semieje_mayor)));
    apoapsis = semieje_mayor*(1-excentricidad^2)/(1 - excentricidad)- R_Kerbin;
    anomalia_verdadera = real(acos(((semieje_mayor*(1-excentricidad^2)/r0)-1)/excentricidad));
    
    if apoapsis > 100000
        fprintf("Final speed %f m/s\n", speed);
        fprintf("Final height %f m/s\n", height);
        fprintf("Final Fuel (L) %f Kg\n", liquid_fuel_mass);
        fprintf("Final Fuel (O) %f Kg\n", oxidizer_fuel_mass);
        fprintf("Apoapsis: %f Km \n", apoapsis*1e-3);
        disp("________");
        break
    end

    if liquid_fuel_mass <= 0 || oxidizer_fuel_mass <= 0
        xline(advanced)
        disp("Out of fuel!");
        fprintf("Max speed %f m/s\n", speed);
        fprintf("Final Fuel (L) %f Kg\n", liquid_fuel_mass);
        fprintf("Final Fuel (O) %f Kg\n", oxidizer_fuel_mass);
        fprintf("Current height: %f m\n", height);
        disp("_____");
        break;
    end

    
end


N = 1000;
plot_orbit(N, semieje_mayor, excentricidad, anomalia_verdadera, R_Kerbin);

plot_orbit(N, semieje_mayor, excentricidad, pi, R_Kerbin);

