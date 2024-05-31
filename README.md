# Water Rocket
Simple project solving the differtial equations for a water rocket using a Cash-Karp Runge-Kutta solver  

## Theory

### A Simplistic Model of a Water Rocket

In the following, we assume a cylindrical rocket with the cross-section area $A_{rocket}$, the nozzle opening area $A_{exaust}$. The rocket has the hight $H_{rocket}$. The rocket can be filled with water and the filling level is counted from the bottom of the rocket and is a function of time $h_{water}(t)$. Hence, the air volume in the rocket can simply be calculated as 
```math
V_{air}(t)=A_{rocket}(H_{rocket}-h_{water}(t))
```

The equation of motion for the rocket, assuming all forces are directed along the z-axis:
```math
F_{rocket} = F_{trust} + F_{gravitation} + F_{drag},
```
where 
- $F_{rocket}=m_{rocket}\frac{dv_{rocket}}{dt}$,
- $F_{thrust}=v_{exaust}\frac{dm_{rocket}}{dt}$,
- $F_{gravitation}=-m_{rocket}g$,
- $F_{drag}=-\frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}\frac{v_{rocket}}{|v_{rocket}|}$

Above, we assume quadratic drag for objects at high Reynolds numbers and since for our simple setup, we do not expect the rocket to fly higher than a few meters, we can safely assume the gravitational acceleration to be constant. Inserting these relations obove, we obtain:
```math
m_{rocket}\frac{dv_{rocket}}{dt} = v_{exaust}\frac{dm_{rocket}}{dt} - m_{rocket}g - \frac{1}{2} \rho_{ambient} v_{rocket}^2 C_D A_{rocket}\frac{v_{rocket}}{|v_{rocket}|}.
```
The mass of the rocket consists of several contributions, $m_{rocket} = m_{hull}+m_{water}+m_{air}\approx m_{rocket} = m_{hull}+m_{water}$, where we neglect for simplicity the contribution from the pressurized air, which should be in the order of a few grams for a 1l bottle and pressures up to five bars.  

While the mass of the rocket hull is constant over time, the mass of the water is a function of time. Further, also the exaust velocity is a function of time, since the pressure in the bottle will decrease as more water or air is expelled from the bottle. 

In summary, we have the the unknown variables $m_{water}$, $v_{rocket}$, $v_{exaust}$, and several known variables such as the cross sectional area of the bottle $A_{rocket}$, the mass of the bottle $m_{hull}$, or the ambient air density $\rho_{ambient}$.

To relate the water pressure in the bottle to the pressure at the nozze, the Bernoulli equation is applied:
```math
p_{water}+\frac{1}{2}\rho_{water} v_{sink}^2=p_{ambient} + \frac{1}{2}\rho_{water} v_{exaust}^2
```
where $v_{sink}$ is the sinking velocity of the water in the bottle and $v_{exaust}$ the exaust velocity of the water at the nozzle. Here, we assume water to be incompressible ($\rho_{water} = const$) and neglect the hydrostatic pressure of the system. Further, we assume that the pressue at the nozze is equal to the ambient air pressue. Here, we introduced two further unknown, $v_{sink}$ and $p_{water}$. To find a closing relation for $v_{sink}$, the mass conservation can be applied:
```math
v_{sink}A_{rocket}=v_{exaust}A_{exaust},
```
where $A_{exaust}$ is the cross sectional area of the nozzle opening. To determine the water pressure $p_{water}$, we assume that the water pressure is in thermodynamic equilibrium with the pressurized air, so that $p_{water}=p_{air}$ holds.

However, with sinking water level in the bottle, also the air pressure is lowering. To relate the air pressure at time $t$ to the initial air pressure at time $t=0$, we assume a polytropic process with the polytropic index $\gamma$:
```math
p_{air}(t)V_{air}^{\gamma}(t) = p_{air}(0)V_{air}^{\gamma}(0)
```





