[← Back to the README](../README.md)   

# FFT filters: Grav/Mag, Frequency and Gradient   

## Line-Noise Filter Update   

The **wedge** (half-width, degrees) of the directional line-noise filter is now a setting on the FFT Filters tab of both versions (default 45), and the ArcGIS Pro toolbox has a **Remove Line Noise** tool with it. **Fix**: the filter previously treated only one of the two opposite lobes of the noise spectrum, so about half of the noise was left in; both are now filtered. Because of that, a Scale of 2 or more now over-subtracts, so **Scale** now defaults to 1.   

## Grav/Mag Filters   
   
**Reduction to the Pole**    
$`H_{RTP}(k_x, k_y) = \frac{k \cos I \cos D + i k_y \cos I \sin D + k_x \sin I}{k}`$   
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic pole.
Where     
- k<sub>x</sub> and k<sub>y</sub> : The wavenumber components in the x and y directions.
- k = The total wavenumber magnitude = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>}   
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit.


**Reduction to the Equator**    
$`H_{RTE}(k_x, k_y) = \frac{k \cos I \cos D + i k_y \cos I \sin D + k_x \sin I}{k \cos I \cos D - i k_y \cos I \sin D + k_x \sin I}`$     
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic equator.
Where   
- k<sub>x</sub> and k<sub>y</sub> : The wavenumber components in the x and y directions.
- k = The total wavenumber magnitude = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>}
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit. 
   
**Variable (Differential) Reduction to the Pole**   
Accounts for the spatial variation of inclination and declination across a survey area by applying a Taylor-series expansion of the standard RTP operator (Cooper & Cowan, 2005). The IGRF field is computed at the four grid corners and the centre, then bilinearly interpolated to build spatially varying inclination and declination grids. The expansion is evaluated at 13 perturbed inc/dec pairs and combined to give a single corrected output.   
Outputs: `_DRTP` (corrected field), `_DRTP_inc` and `_DRTP_dec` (the interpolated inclination and declination fields used).   
Recommended for surveys covering several degrees of latitude/longitude where a single inc/dec value would introduce systematic phase or amplitude errors.   
Based on: Cooper, G.R.J. & Cowan, D.R. (2005), *Computers & Geosciences* 31, 989–999. https://doi.org/10.1016/j.cageo.2005.02.005   

**Continuation**    
$`H(k) = e^{-k h}`$   
Where   
h > 0 for upward continuation.   
h < 0  for downward continuation.   
   
**Vertical Integration**   
$`H(k_x, k_y) = \frac{1}{k}`$  
When applied to an RTE or RTP image provides the so called Pseudogravity result    
Where    
k = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>} .   
   

## Frequency Filters   

**What the cutoff values mean**   
The values entered for the High pass, Low pass and Band pass filters are **wavelengths**, in the length units of the grid (metres, or degrees for a geographic grid), not frequencies. A wavelength is the distance over which a feature repeats, so it is the inverse of frequency: *short* wavelength = *high* frequency (small, sharp, usually shallow features) and *long* wavelength = *low* frequency (broad, smooth features). The equations below use the wavenumber $k = 2\pi/\lambda$.   

- **High pass** keeps wavelengths *shorter* than the cutoff and removes longer ones (regionals, deeper sources).   
- **Low pass** keeps wavelengths *longer* than the cutoff and removes shorter ones (noise, shallow sources).   
- **Band pass** keeps the wavelengths between two values. **Low cut wavelength** is the *longest* wavelength kept (longer wavelengths, i.e. lower frequencies, are cut) and **High cut wavelength** is the *shortest* wavelength kept (shorter wavelengths, i.e. higher frequencies, are cut), so the Low cut value is the larger number (for example 50000 and 5000).   
- **Transition width** (same units) is how gradually the filter switches from pass to reject: 0 is a sharp cut, and larger values reduce ringing.   

*Relationship to source depth (approximate).* The signal from a source at depth $z$ falls off with wavelength $\lambda$ as $e^{-2\pi z/\lambda}$, the same factor that governs upward continuation, so a deep source cannot produce short wavelengths. A source at depth $z$ typically has its main response at wavelengths of about 3 to 6 times $z$. As a rough guide, a cutoff of 6000 m separates sources shallower than about 1000 to 2000 m (kept by a high pass) from deeper ones (kept by a low pass), and a band pass between $\lambda_{high}$ and $\lambda_{low}$ favours sources from about $\lambda_{high}/6$ to $\lambda_{low}/3$. This is only a guide: a shallow, broad source also has long wavelengths, and potential-field data are inherently ambiguous, so filtering by wavelength is not a true depth slice.   

**High Pass Filter**

$$H(k) = \begin{cases} 
0 & \text{if } k \leq k_{low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{low}}{k_{high} - k_{low}}\right)\right) & \text{if } k_{low} < k < k_{high} \\
1 & \text{if } k \geq k_{high}
\end{cases}$$

The high-pass filter removes low-frequency components (long wavelengths) while retaining high-frequency components (short wavelengths) with a smooth transition to reduce ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{low} = \frac{2\pi}{\lambda_c + w/2}$ : Lower transition boundary
- $k_{high} = \frac{2\pi}{\lambda_c - w/2}$ : Upper transition boundary  
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (in same units as wavelength)

**Low Pass Filter**

$$H(k) = \begin{cases} 
1 & \text{if } k \leq k_{inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{inner}}{k_{outer} - k_{inner}}\right)\right) & \text{if } k_{inner} < k < k_{outer} \\
0 & \text{if } k \geq k_{outer}
\end{cases}$$

The low-pass filter removes high-frequency components (short wavelengths) while retaining low-frequency components (long wavelengths). Optional smooth transition reduces potential ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{inner} = \frac{2\pi}{\lambda_c + w/2}$ : Inner transition boundary
- $k_{outer} = \frac{2\pi}{\lambda_c - w/2}$ : Outer transition boundary
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (optional, in same units as wavelength)

**Band Pass Filter**

$$H(k) = H_{high}(k) \times H_{low}(k)$$

Where:

***High-pass component:***
$$H_{high}(k) = \begin{cases} 
0 & \text{if } k \leq k_{h,low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{h,low}}{k_{h,high} - k_{h,low}}\right)\right) & \text{if } k_{h,low} < k < k_{h,high} \\
1 & \text{if } k \geq k_{h,high}
\end{cases}$$

***Low-pass component:***
$$H_{low}(k) = \begin{cases} 
1 & \text{if } k \leq k_{l,inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{l,inner}}{k_{l,outer} - k_{l,inner}}\right)\right) & \text{if } k_{l,inner} < k < k_{l,outer} \\
0 & \text{if } k \geq k_{l,outer}
\end{cases}$$

The band-pass filter isolates features within a specific wavelength range by combining high-pass and low-pass components with smooth transitions.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{h,low} = \frac{2\pi}{\lambda_{low} + w_h/2}$ : High-pass lower transition boundary
- $k_{h,high} = \frac{2\pi}{\lambda_{low} - w_h/2}$ : High-pass upper transition boundary
- $k_{l,inner} = \frac{2\pi}{\lambda_{high} + w_l/2}$ : Low-pass inner transition boundary
- $k_{l,outer} = \frac{2\pi}{\lambda_{high} - w_l/2}$ : Low-pass outer transition boundary
- $\lambda_{low}$ : Low cutoff wavelength (removes longer wavelengths)
- $\lambda_{high}$ : High cutoff wavelength (removes shorter wavelengths)
- $w_h$ : High-pass transition width
- $w_l$ : Low-pass transition width

**Directional Band Pass**   
Removes combined directional and high pass filtered data from original data, with scaling function modify extent of feature suppression.   
   
***Butterworth High-Pass Filter***
$`H(k) = \frac{1}{1 + \left(\frac{k_c}{k}\right)^{2n}}`$    
The Butterworth filter attenuates frequencies below the cutoff k<sub>c</sub> while preserving higher frequencies.    
H(k) : Filter response as a function of wavenumber k.    
k : Wavenumber (spatial frequency).    
k<sub>c</sub> : Cutoff wavenumber, related to the cutoff wavelength by k<sub>c</sub> = \frac{1}{\text{cutoff wavelength}}.    
n : Filter order, determining the sharpness of the transition. Higher \( n \) makes the filter more selective.   
    
***Directional Cosine Filter***    
$`H(k_x, k_y) = \left| \cos(\theta - \theta_c) \right|^p`$   
The Directional Cosine Filter emphasizes or suppresses frequency components along a specific direction.   
H(k<sub>x</sub>, k<sub>y</sub>): Filter response as a function of wavenumber components k<sub>x</sub> and k<sub>y</sub>.   
theta = \arctan\left(\frac{k_y}{k_x}\right) : Angle of the frequency component.   
theta<sub>c</sub> : Center direction (in radians), representing the direction to emphasize.   
p : Degree of the cosine function. Higher \( p \) sharpens the directional emphasis.   

**Remove Regional**   
Remove a 1st order (dipping plane) or 2nd order (parabolic plane) regional from data. 
   
**Automatic Gain Control**    
$`AGC(x, y) = \frac{f(x, y)}{\text{RMS}(f(x, y), w)}`$   
Where    
RMS(f, w)  is the root mean square of the data over a window w.   
   
**Radially averaged power spectrum (but needs testing!)**    
$`P(k) = \frac{1}{N_k} \sum_{(k_x, k_y) \in k} |\text{FFT}(f)|^2`$   
Where    
P(k) is the radially averaged power spectrum, and N<sub>k</sub> is the number of samples in the radial bin.   
   

## Gradient Filters   
   
**Derivative**    
$`\frac{\partial f}{\partial u} = \frac{\partial f}{\partial x} \cos\theta + \frac{\partial f}{\partial y} \sin\theta`$   
Where   
theta is the angle defining the direction of the derivative (x,y or z).   
The vertical derivative is calculated from the horizontal gradients using Laplace's equation, $\frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} + \frac{\partial^2 f}{\partial z^2} = 0$, which gravity and magnetic fields satisfy outside their sources (they are Green's functions). In the Fourier domain this multiplies the spectrum by $k$, so the vertical derivative amplifies short wavelengths and often amplifies noise as well.   
   
**Total Horizontal Gradient**   
$`THG(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}`$   
   
**Analytic Signal**    
$`A(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2 + \left(\frac{\partial f}{\partial z}\right)^2}`$   
Computes the total amplitude of the gradients, independent of field inclination or declination.
Useful for locating edges of potential field sources (e.g., faults or contacts).   
       
**Tilt Angle**    
$`T = \tan^{-1}\left(\frac{\frac{\partial f}{\partial z}}{\sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}}\right)`$   
Enhances the contrast of geological features by highlighting gradients relative to the vertical component. The result is saved in degrees (-90 to 90).
Where   
df/dz : Vertical derivative of the field.
df/dx , df/dy : Horizontal derivatives of the field.   
   
