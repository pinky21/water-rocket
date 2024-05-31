# Water Rocket
Simple project solving the differtial equations for a water rocket using a Cash-Karp Runge-Kutta solver  

## Theory

### A Simplistic Model of a Water Rocket
The equation of motion for the rocket, assuming all forces are directed along the z-axis:
```math
F_{rocket} = F_{trust} + F_{gravitation} + F_{drag},
```
where 
- $F_{rocket}=m_{rocket}\frac{dv_{rocket}}{dt}$,
- $F_{thrust}=v_{exaust}\frac{dm_{rocket}}{dt}$,
- $F_{gravitation}=-m_{rocket}g$,
- $F_{drag}=-\frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}\frac{v_{rocket}}{|v_{rocket}|}$, assuming quadratic drag for objects at high Reynolds numbers.

Since for our simple setup, we do not expect the rocket to fly higher than a few meters, we can safely assume the gravitational acceleration to be constant. Inserting these relations obove, we obtain:
```math
m_{rocket}\frac{dv_{rocket}}{dt} = v_{exaust}\frac{dm_{rocket}}{dt} - m_{rocket}g - \frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}\frac{v_{rocket}}{|v_{rocket}|}.
```
The mass of the rocket consists of several contributions, $m_{rocket} = m_{hull}+m_{water}+m_{air}\approx m_{rocket} = m_{hull}+m_{water}$, where we neglect for simplicity the contribution from the pressurized air, which should be in the order of a few grams for a 1l bottle and pressures up to five bars.  

While the mass of the rocket hull is constant over time, the mass of the water is a function of time. Further, also the exaust velocity is a function of time, since the pressure in the bottle will decrase as more water or air is expelled from the bottle. To relate the water pressure in the bottle to the pressure at the nozze, the Bernoulli equation is applied:
```math
p_{water}+\frac{1}{2}\rho_{water} v_{sink}^2=p_{ambient} + \frac{1}{2}\rho_{water} v_{exaust}^2
```
where $v_{sink}$ is the sink velocity of the water in the bottle and v_{exaust} the exaust velocity of the water at the nozzle. Here, we assume water to be incompressible and neglect the hydrostatic pressure of the system.


