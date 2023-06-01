clear all;
clc;

%1.5
%5
%10
target_payload      = 12.5;% [Metric Tons] 3mf

max_jet_count       = 10;
max_rocket_count    = 10;
min_jet_count       = 1;
min_rocket_count    = 1;

iteration_relaxation = 0.8;
strut_dens_liq_dens = 0.3;

G           = 6.674e-11;
M_Kerbin    = 5.2915158e22;
R_Kerbin    = 6e5;
mu          = 3.5316000e12;
g           = 9.81; % [m/s^2]
runaway_length = 2000 ;% [m]
cd0 = 0.06;
cl0 = 0.0;
k = 0.01;
cl_alpha = 2*pi;

iteration_added_jet_fuel_mass       = 0; % [kg]
iteration_added_liquid_fuel_mass    = 0; % [kg]
iteration_added_oxidizer_fuel_mass  = 0; % [kg]

h_speed = 0; % [m/s]
v_speed = 0; % [m/s]
height  = 0; % [m]
advanced= 0; % [m]

dt = .1; % [s]
t = 0;

takeoff = 1; % Used to know if have taken off
archived_terminal_vel = 1;
got_inclined = 0;

alpha = 0; % [rad]
vehicle_angle = 0; % [rad]
target_vehicle_angle = 30*pi/180;
angle_change_velocity= 0.05; % rad/s

current_jet_eng_name = 0;
current_rocket_eng_name = 0;

current_n_jet_engine = min_jet_count;
current_n_rocket_engine = min_rocket_count;

current_jet_fuel_mass               = 0; % [kg]
curent_liquid_fuel_mass             = 0; % [kg]
current_oxidizer_fuel_mass          = 0; % [kg]

target_payload = target_payload * 1000;

current_configuration = 1;
configurations = [0 0 0 0 0 0 0 0 0];
while true

current_iteration = 1;
    while true
        
        disp("###############");
        fprintf("Current iteration: %d\n", current_iteration);
        disp("###############");
        fprintf("Jet count %d\n", current_n_jet_engine);
        fprintf("Rocket count %d \n", current_n_rocket_engine);
        
        fprintf("Jet name %d\n", current_jet_eng_name);
        fprintf("Rocket name %d \n", current_rocket_eng_name);
        
    
        h_speed = 0; % [m/s]
        v_speed = 0; % [m/s]
        height  = 0; % [m]
        advanced= 0; % [m]
        
        dt = .1; % [s]
        t = 0;
    
        [h, t_e, p_e, rho_e, c] = Earth_to_Kerbin(height/1000);
        p_e = p_e / 101325;
    
        [eng_mass_jet,fuel_consumpt_air,air_thrust,vel_mult,atm_mult] = thrust_curves_jet(current_jet_eng_name,0,p_e);
        [eng_mass_rocket,liquid_consumpt_rocket,oxidizer_consumpt_rocket,rocket_thrust] = thrust_curves_rocket(current_rocket_eng_name,p_e);
        
        current_struct_mass = (current_jet_fuel_mass + curent_liquid_fuel_mass + current_oxidizer_fuel_mass) * strut_dens_liq_dens;
        total_dry_weight = current_struct_mass + ...
            eng_mass_jet * current_n_jet_engine + eng_mass_rocket * current_n_rocket_engine + ...
            target_payload;
    
        total_fuel_weight = current_jet_fuel_mass + curent_liquid_fuel_mass + current_oxidizer_fuel_mass;
    
        total_weight = total_dry_weight + total_fuel_weight;
    
        air_acceleration = air_thrust/total_dry_weight; % [m/s^2]
    
        perceived_alpha = 10*pi/180; % [rads]
    
        Cl = cl0 + cl_alpha * perceived_alpha;
    
        time_to_accelerate = sqrt(2*runaway_length/air_acceleration); % [s]
        end_of_runaway_speed = air_acceleration * time_to_accelerate;
        S = total_weight * g * 2 / (end_of_runaway_speed^2 * Cl * rho_e);
    
        %% Starts aerodynamic simulation:
        iteration_added_jet_fuel_mass       = 0; % [kg]
    
        available_jet_fuel_mass         = current_jet_fuel_mass;
        available_liquid_fuel_mass      = curent_liquid_fuel_mass;
        available_oxidizer_fuel_mass    = current_oxidizer_fuel_mass;
        iteration_added_strut_mass      = 0;
        aerodynamic_iterations = 0;
    
        takeoff = 1; % Used to know if have taken off
        archived_terminal_vel = 1;
        got_inclined = 0;
    
        status = 0; % Check if iteration gone wrong
    
        while true
            [h, t_e, p_e, rho_e, c] = Earth_to_Kerbin(height/1000);
            p_e = p_e / 101325;
            speed = sqrt(v_speed^2 + h_speed^2);
            M = speed/c;
            [eng_mass_jet,fuel_consumpt_air,air_thrust,vel_mult,atm_mult] = thrust_curves_jet(current_jet_eng_name,M,p_e);
            [eng_mass_rocket,liquid_consumpt_rocket,oxidizer_consumpt_rocket,rocket_thrust] = thrust_curves_rocket(current_rocket_eng_name,p_e);
        
            real_th = current_n_jet_engine*air_thrust*vel_mult*atm_mult;
            
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
            Cd = cd0;%+ k*Cl^2;
        
            L = Cl * rho_e * speed^2 * S / 2;
            D = Cd * rho_e * speed^2 * S / 2;
            Tm = total_dry_weight + ...
                available_oxidizer_fuel_mass + ...
                available_jet_fuel_mass + ...
                available_liquid_fuel_mass + ...
                iteration_added_strut_mass + ...
                eng_mass_jet * current_n_jet_engine + ...
                eng_mass_rocket * current_n_rocket_engine;
    
            W = Tm * g;
            
            if (L >= W)
                if (takeoff == 1)
                    disp("Takeoff!");
                    fprintf("Takeoff speed: %f m/s\n", speed);
                    fprintf("Takeoff distance: %f m\n", advanced);
                    %disp("_____");
                    takeoff = 0;
                end
            end
        
            if available_jet_fuel_mass<= 0
                iteration_added_jet_fuel_mass  = iteration_added_jet_fuel_mass + ...
                    fuel_consumpt_air * current_n_jet_engine * dt;
                iteration_added_strut_mass = iteration_added_strut_mass + fuel_consumpt_air * current_n_jet_engine * dt * strut_dens_liq_dens;
            end
            
            total_front_force = real_th*cos(vehicle_angle) - D*cos(vehicle_angle) - L*sin(vehicle_angle);
            total_vertical_force = real_th*sin(vehicle_angle) + L*cos(vehicle_angle) -D*sin(vehicle_angle) - W;
            dhV = total_front_force / Tm;
            dvV = total_vertical_force / Tm;
        
            if (dvV < 0) && (height == 0) % Aún no ha arrancado
                dvV = 0;
            end
        
            new_hspeed = h_speed + dhV * dt;
            new_vspeed = v_speed + dvV * dt;
        
            new_speed = sqrt(new_vspeed^2 + new_hspeed^2);
        
        
            if(archived_terminal_vel)
                if dhV < 0.1
                    fprintf("Reached terminal speed: %f m/s\n", speed);
                    fprintf("Remaining fuel: %f kg\n", available_jet_fuel_mass);
                    fprintf("Current height: %f m\n", height);
                    archived_terminal_vel = 0;
                    if takeoff == 1
                        status = 1;
                        break
                    end
                    disp("_____");
                end
            end
            
            if speed_angle < target_vehicle_angle/2
                if (got_inclined && archived_terminal_vel == 0)
                    fprintf("Reached terminal height: %f m/s\n", speed);
                    fprintf("Remaining fuel (Liquid): %f kg\n", available_jet_fuel_mass);
                    fprintf("Current height: %f m\n", height);
                    disp("_____");
                    break
                end
            end
            h_speed = new_hspeed;
            v_speed = new_vspeed;
            if height < -100
                disp(height)
                status = 1;
                break
            end
            advanced = advanced + h_speed*dt;
            height = height + v_speed*dt;
            if available_jet_fuel_mass > 0
                available_jet_fuel_mass = available_jet_fuel_mass - (current_n_jet_engine * fuel_consumpt_air * dt);
            end
        end
    
        if status == 1
            %d
            disp("Found inviable configuration");
            if current_n_jet_engine == max_jet_count
                if current_jet_eng_name == 2
                    break
                else
                    current_jet_eng_name = current_jet_eng_name + 1;
                    current_n_jet_engine = min_jet_count;
                end
                current_rocket_eng_name = 0;
                current_n_rocket_engine = min_rocket_count;
            else
                current_n_jet_engine = current_n_jet_engine + 1; 
            end
        
            current_oxidizer_fuel_mass = 0;
            curent_liquid_fuel_mass = 0;
            current_jet_fuel_mass = 0;
            continue
        end
    
    
        
        current_jet_fuel_mass = current_jet_fuel_mass + iteration_added_jet_fuel_mass*iteration_relaxation;
        if available_jet_fuel_mass > 0
            current_jet_fuel_mass = current_jet_fuel_mass - available_jet_fuel_mass*(1-iteration_relaxation);
        end
        iteration_added_jet_fuel_mass = 0;
    
        %% Starts orbital iterations
        %% Parte aerodinámica DONE, empieza la orbital
    
        target_vehicle_angle = 30 * pi/180;
        h_speed = h_speed + 174; %% Adding sideral rotational velocity
        
        while true
            [h, t_e, p_e, rho_e, c] = Earth_to_Kerbin(height/1000);
            p_e = p_e / 101325;
            
            [eng_mass_rocket,liquid_consumpt_rocket,oxidizer_consumpt_rocket,rocket_thrust] = thrust_curves_rocket(current_rocket_eng_name,p_e);
            
            real_th = current_n_rocket_engine * rocket_thrust * 1e3;
            
            if vehicle_angle < target_vehicle_angle
                vehicle_angle = vehicle_angle + angle_change_velocity*dt;
            else
                vehicle_angle = target_vehicle_angle;
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
            Tm = total_dry_weight + ...
                available_oxidizer_fuel_mass + ...
                available_jet_fuel_mass + ...
                available_liquid_fuel_mass + ...
                iteration_added_strut_mass + ...
                eng_mass_jet * current_n_jet_engine + ...
                eng_mass_rocket * current_n_rocket_engine;
    
        
            W = Tm * g;
        
            r0 = height+R_Kerbin;
        
            total_front_force = real_th*cos(vehicle_angle) - D*cos(vehicle_angle) - L*sin(vehicle_angle);
            total_vertical_force = real_th*sin(vehicle_angle) + L*cos(vehicle_angle) -D*sin(vehicle_angle) - W;
            dhV = total_front_force / Tm;
            dvV = total_vertical_force / Tm;
            dvV = dvV + h_speed^2/r0;
        
        
            new_hspeed = h_speed + dhV * dt;
            new_vspeed = v_speed + dvV * dt;
        
            new_speed = sqrt(new_vspeed^2 + new_hspeed^2);
        
            h_speed = new_hspeed;
            v_speed = new_vspeed;
        
            height = height + v_speed*dt;
            
            if available_liquid_fuel_mass > 0
                available_liquid_fuel_mass = available_liquid_fuel_mass - current_n_rocket_engine * liquid_consumpt_rocket * dt;
            end
    
            if available_oxidizer_fuel_mass > 0
                available_oxidizer_fuel_mass  = available_oxidizer_fuel_mass - current_n_rocket_engine * oxidizer_consumpt_rocket * dt;
            end
        
            if speed > 0
                speed_angle = atan(v_speed/h_speed);
            end
    
            if v_speed < 0
                status = 1;
                break
            end
    
            energia_cinetica = 1/2 * speed^2 - (mu/(r0));
        
            semieje_mayor = -mu/(2*energia_cinetica);
            momento_angular = r0*speed*cos(speed_angle);
            excentricidad = sqrt(1-((momento_angular^2)/(mu*semieje_mayor)));
            apoapsis = semieje_mayor*(1-excentricidad^2)/(1 - excentricidad)- R_Kerbin;
            anomalia_verdadera = real(acos(((semieje_mayor*(1-excentricidad^2)/r0)-1)/excentricidad));
            
            if apoapsis > 100000
                fprintf("Final speed %f m/s\n", speed);
                fprintf("Final height %f m/s\n", height);
                fprintf("Final Fuel (L) %f Kg\n", available_liquid_fuel_mass);
                fprintf("Final Fuel (O) %f Kg\n", available_oxidizer_fuel_mass);
                fprintf("Apoapsis: %f Km", apoapsis*1e-3);
                disp("________");
                break
            end
    
            if available_oxidizer_fuel_mass <= 0
                iteration_added_oxidizer_fuel_mass  = iteration_added_oxidizer_fuel_mass + ...
                    oxidizer_consumpt_rocket * current_n_rocket_engine * dt;
                iteration_added_strut_mass = iteration_added_strut_mass + oxidizer_consumpt_rocket * current_n_rocket_engine * dt * strut_dens_liq_dens;
            end
            if available_liquid_fuel_mass <= 0
                iteration_added_liquid_fuel_mass  = iteration_added_liquid_fuel_mass + ...
                    liquid_consumpt_rocket * current_n_rocket_engine * dt;
                iteration_added_strut_mass = iteration_added_strut_mass + liquid_consumpt_rocket * current_n_rocket_engine * dt * strut_dens_liq_dens;
            end        
        end
        
        if status == 1
            disp("Unviable configuration during orbital phase")
            if current_n_rocket_engine == max_rocket_count
                if current_rocket_eng_name == 7
                    if current_jet_eng_name == 2
                        break
                    else
                        current_jet_eng_name = current_jet_eng_name +1;
                        current_n_jet_engine = min_jet_count;
                        current_n_rocket_engine = min_rocket_count;
                        current_rocket_eng_name = 0;
                    end
                else
                    current_rocket_eng_name = current_rocket_eng_name + 1;
                    current_n_rocket_engine = min_rocket_count;
                    current_n_jet_engine = min_jet_count;
                end
            else
                current_n_rocket_engine = current_n_rocket_engine + 1;
            end
    
            current_oxidizer_fuel_mass = 0;
            curent_liquid_fuel_mass = 0;
            current_jet_fuel_mass = 0;
            continue
        end
    
        curent_liquid_fuel_mass = curent_liquid_fuel_mass + iteration_added_liquid_fuel_mass*iteration_relaxation;
        current_oxidizer_fuel_mass = current_oxidizer_fuel_mass + iteration_added_oxidizer_fuel_mass*iteration_relaxation;
        if available_liquid_fuel_mass > 0
            curent_liquid_fuel_mass = curent_liquid_fuel_mass - available_liquid_fuel_mass*(1-iteration_relaxation);
        end
    
        if available_oxidizer_fuel_mass > 0
            current_oxidizer_fuel_mass = current_oxidizer_fuel_mass - available_oxidizer_fuel_mass*(1-iteration_relaxation);
        end
        iteration_added_liquid_fuel_mass = 0;
        iteration_added_oxidizer_fuel_mass = 0;
    
        if iteration_added_strut_mass < 100
            break
        end
        current_iteration = current_iteration + 1;
        iteration_added_strut_mass = 0;
    end
    
    disp("final configuration");
    fprintf("Jet fuel mass %f [kg]\n", current_jet_fuel_mass);
    fprintf("Liquid fuel mass %f [kg]\n", curent_liquid_fuel_mass);
    fprintf("Oxidier fuel mass %f [kg]\n", current_oxidizer_fuel_mass);
    fprintf("Strut mass %f [kg]\n", current_struct_mass);
    
    fprintf("Jet count %d\n", current_n_jet_engine);
    fprintf("Rocket count %d \n", current_n_rocket_engine);
    
    fprintf("Jet name %d\n", current_jet_eng_name);
    fprintf("Rocket name %d \n", current_rocket_eng_name);
    
    fprintf("Aicraft wing %f [m^2]\n", S);

    configurations(current_configuration,1) = current_jet_fuel_mass;
    configurations(current_configuration,2) = curent_liquid_fuel_mass;
    configurations(current_configuration,3) = current_oxidizer_fuel_mass;
    configurations(current_configuration,4) = current_struct_mass;
    configurations(current_configuration,5) = current_n_jet_engine;
    configurations(current_configuration,6) = current_n_rocket_engine;
    configurations(current_configuration,7) = current_jet_eng_name;
    configurations(current_configuration,8) = current_rocket_eng_name;
    configurations(current_configuration,9) = S;

    current_configuration = current_configuration + 1;
    
    if current_rocket_eng_name == 7
        if current_jet_eng_name == 2
            break
        else
            current_jet_eng_name = current_jet_eng_name + 1;
            current_rocket_eng_name = 0;
            
        end
    else
        current_rocket_eng_name  = current_rocket_eng_name  + 1;
    end

    current_oxidizer_fuel_mass = 0;
    curent_liquid_fuel_mass = 0;
    current_jet_fuel_mass = 0;

    current_n_rocket_engine = min_rocket_count;
    current_n_jet_engine    = min_jet_count;
end

JetFuel         = configurations(:,1);
LiquidFuel      = configurations(:,2);
OxidizerFuel    = configurations(:,3);
StructMass      = configurations(:,4);
N_Jets          = configurations(:,5);
N_Rockets       = configurations(:,6);
Jet_ID          = configurations(:,7);
Rocket_ID       = configurations(:,8);
Surface         = configurations(:,9);

T = table(JetFuel, LiquidFuel, OxidizerFuel, StructMass, ...
    N_Jets, N_Rockets, Jet_ID, Rocket_ID, Surface);
output_file_base = "payload_%.2f_t.csv";
output_file = sprintf(output_file_base,target_payload/1000);
writetable(T,output_file)  