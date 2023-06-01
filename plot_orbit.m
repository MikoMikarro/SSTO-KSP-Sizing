function [] = plot_orbit(N,semieje_mayor, excentricidad, anomalia_verdadera, R_planet)
%PLOT_ORBI Summary of this function goes here
%   Detailed explanation goes here

figure();
angles = 1:N;
angles = angles .* (2*pi/N);
X_Kerbin = zeros(1,N);
Y_Kerbin = zeros(1,N);
X_My_Orbit = zeros(1,N);
Y_My_Orbit = zeros(1,N);
X_radius = zeros(1,2);
Y_radius = zeros(1,2);

R_radius = semieje_mayor*(1-excentricidad^2)/(1+ excentricidad * cos(anomalia_verdadera));
X_radius(2) = R_radius * cos(anomalia_verdadera);
Y_radius(2) = R_radius * sin(anomalia_verdadera);
for i=1:N
    X_Kerbin(i) = R_planet * cos(angles(i));
    Y_Kerbin(i) = R_planet * sin(angles(i));

    %r_0 = a(1-e^2)/(1+e cos(thet))
    R_My_Orbit = semieje_mayor*(1-excentricidad^2)/(1+ excentricidad * cos(angles(i)));
    X_My_Orbit(i) = R_My_Orbit * cos(angles(i));
    Y_My_Orbit(i) = R_My_Orbit * sin(angles(i));
end

hold on;
plot(X_Kerbin, Y_Kerbin);
plot(X_My_Orbit, Y_My_Orbit);
plot(X_radius, Y_radius);
daspect([1 1 1])
hold off;

end

