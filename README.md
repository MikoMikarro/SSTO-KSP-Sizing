# SSTO-KSP-Sizing

Simple matlab code that computes the necessary fuel, engines and surface area of an SSTO to get to Orbit (and circularization if added some fuel ; D )

Code tested on KSP 1.12.4

### How to use it

Use the `optimizer.m` code to generate the different viable designs by using the following parameters as inputs

* `target_payload` : The aumount of Metric Tons you want to throw at orbit
* `max_jet_count, min_jet_count, min_rocket_count, max_rocket_count` : Recommended to have `min_jet_count` at 2 for symetric vessels
* `iteration_relaxation` :Reduces de speed of the algorithm but makes it more stable
* `strut_dens_liq_dens` : fraction of kilograms of structural mass expected to be added from any kilogram of fuel, recomended values are:
  * 0.3 - For small vessles 		 (< 5 t)
  * 0.2 - For medium vessels 	 (    5 T < 15 T)
  * 0.1 - For large vessels		 ( > 15 T)

This program will return a table with all the possible designs that archive this feat. Use the one you like more! Use the values it gives as design guidance for your SSTO.

The name of the Rockets and Jets can be found in the `thrust_curves_jet` and `thrust_curves_rocket.m` files.

### Flight assistant

Once you've designed your aircract use the real values of your aircraft and input them tothe `plotter.m`

The algorithm will give you 3 plots and a flight guide.

The plots are visuals of the expected trayectory of the flight if the guide is followed. For ex:

```
Takeoff!
Takeoff speed: 305.844348 m/s
Takeoff distance: 2248.705431 m
Put your vehicle in 0º of pitch
_____
Reached terminal speed: 1794.998978 m/s
Remaining fuel: 1195.031340 kg
Current height: 5879.712005 m
Put your vehicle in 30º of pitch
_____
Reached terminal height: 1522.313406 m/s
Remaining fuel (Liquid): 703.586340 kg
Current height: 52733.050643 m
Put your vehicle in 30º of pitch and enable rocket boosters!
_____
Final speed 1968.315461 m/s
Final height 63162.301676 m/s
Final Fuel (L) 79.735680 Kg
Final Fuel (O) 375.303720 Kg
Apoapsis: 100.001099 Km
```

# Have fun! 🚀
