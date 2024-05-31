# Water Rocket
Simple project solving the differtial equations for a water rocket using a Cash-Karp Runge-Kutta solver  

## Theory

### A Simplistic Model of a Water Rocket
The equation of motion for the rocket, assuming all forces are directed along the z-axis:
```math
F_{rocket} = F_{trust} + F_{gravitation} + F_{drag},
```
where $F_{rocket}=m_{rocket}\frac{dv_{rocket}}{dt}$, $F_{thrust}=v_{exaust}\frac{dm_{rocket}}{dt}$, $F_{gravitation}=m_{rocket}g$, and $F_{drag}=\frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}$, assuming quadratic drag for objects at high Reynolds numbers. Inserting these relations obove, we obtain:
```math
m_{rocket}\frac{dv_{rocket}}{dt} = v_{exaust}\frac{dm_{rocket}}{dt} + m_{rocket}g + \frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}.
```
The mass of the rocket consists of several contributions, $m_{rocket} = m_{hull}+m_{water}+m_{air}\approx m_{rocket} = m_{hull}+m_{water}$, where we neglect for simplicity the contribution from the pressurized air, which should be in the order of a few grams for a 1l bottle and pressures up to five bars.  

