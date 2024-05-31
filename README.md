# Water Rocket
Simple project solving the differtial equations for a water rocket using a Cash-Karp Runge-Kutta solver  

## Theory

### A Simplistic Model of a Water Rocket
The equation of motion for the rocket, assuming all forces are directed along the z-axis:
```math
F_{rocket} = F_{trust} + F_{gravitation} + F_{drag},
```
where $F_{rocket}=m_{rocket}\frac{dv_{rocket}}{dt}$, $F_{thrust}=v_{exaust}\frac{dm_{rocket}}{dt}$, $F_{gravitation}=m_{rocket}g$, $F_{drag}=0.5 \rho_{ambient} v_{rocket}^2 C_D A_{rocket}$.
