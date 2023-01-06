# Microscale Models


- describe in three dimensions the charge and mass balance within solid particles and within the electrolyte separately
- can be used to simulate small volumes inside electrodes, comprising small clusters of particles
- give information (idea) of how particle sizes and mixtures of geometries etc. interact
- **not** feasable to simulate an entire cell using microscale models

5 important equations:

1. charge conservation in homogenous solid:
$$\nabla \cdot {\bf i_s}= \nabla \cdot (- \sigma \nabla \phi_s) = 0$$

2. mass conservation in homogenous solid:
$$\frac{\partial c_s}{\partial t} = \nabla \cdot (D_s \nabla \phi_s)$$

3. mass conservation in the homogenous electrolyte:
$$\frac{\partial c_e}{\partial t} = \nabla \cdot (D_e \nabla \phi_e) - \frac{{\bf i_e} \cdot \nabla t^0_+}{z_+ \nu_+ F} - \nabla \cdot (c_e {\bf v_0})$$

4. charge conservation in the homogenous electrolyte:
$$\nabla \cdot {\bf i_e} = \nabla \cdot \left( \kappa \nabla \phi_e - \frac{\nu \kappa RT}{F} \left( \frac{s_+}{n \nu_+} + \frac{s_0 c}{n c_0} + \frac{t^0_+}{\nu_+ z_+} \right) \left(1 + \frac{\partial ln f_\pm}{\partial ln c_e} \right) \nabla ln c_e \right) = 0$$

5. lithium movement between the solid and electrolyte phase
$$j = \frac{i_0}{F} \left( \exp{ \left( \frac{(1 - \alpha) n F}{RT} \eta \right) } -  \exp{ \left( -\frac{\alpha n F}{RT} \eta \right) } \right)$$

with:
- ${\bf i_s}$
- $n$
- $s_+$, $s_0$
- $c, c_0$
- $z_+$
- 