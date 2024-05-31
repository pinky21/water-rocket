# Water Rocket
Simple project solving the differtial equations for a water rocket using a Cash-Karp Runge-Kutta solver  

## Theory

### A Simplistic Model of a Water Rocket
The equation of motion for the rocket, assuming all forces are directed along the z-axis:
```math
F_{rocket} = F_{trust} + F_{gravitation} + F_{drag},
```
where `F_{rocket}=m_{rocket}\frac{dv_{rocket}}{dt}`
