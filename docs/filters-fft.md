[← Back to the README](../README.md)   

# FFT filters: Grav/Mag, Frequency and Gradient   

Each function below gives its primary reference (*Reference*, listed in full at the [end of the page](#references)) and what SGTool's code does with it (*In SGTool*): whether the code is original, translated or adapted from someone else's, and where it differs from the published method.   

## Line-Noise Filter Update   

The **wedge** (half-width, degrees) of the directional line-noise filter is now a setting on the FFT Filters tab of both versions (default 45), and the ArcGIS Pro toolbox has a **Remove Line Noise** tool with it. **Fix**: the filter previously treated only one of the two opposite lobes of the noise spectrum, so about half of the noise was left in; both are now filtered. Because of that, a Scale of 2 or more now over-subtracts, so **Scale** now defaults to 1.   

## How the FFT filters are applied   

*In SGTool* (original code): every FFT filter subtracts the grid's median, fills no-data cells with the value of the nearest valid cell, pads the grid to a size whose only factors are 2, 3 and 5 with a linear ramp to zero, tapers the padded grid towards its edges with a Tukey (cosine) window, multiplies the spectrum by the filter, transforms back, removes the padding, restores the no-data cells and adds the median back (except for derivatives, the line-noise estimate and apparent petrophysics, which are zero-mean). The padding on each side is the grid's smaller dimension, up to the Max buffer setting. RTE and the low-latitude RTP pad by mirroring the edges instead (see below). In the equations, $`k = \sqrt{k_x^2 + k_y^2}`$ is the wavenumber in radians per grid unit, and $`\theta`$ the azimuth of the wavenumber, clockwise from north, $`\theta = \arctan(k_x / k_y)`$.   

## Grav/Mag Filters   
   
**Reduction to the Pole**    
$`H_{RTP}(\theta) = \frac{1}{[\sin I + i \cos I \cos(D-\theta)]^2} = \frac{k^2}{[k \sin I + i \cos I \,(k_x \sin D + k_y \cos D)]^2}`$   
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic pole, assuming the magnetisation is induced (parallel to the field).
Where     
- k<sub>x</sub> and k<sub>y</sub> : The wavenumber components in the x and y directions.
- k = The total wavenumber magnitude = $\sqrt{k_x^2 + k_y^2}$   
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit.

*Reference*: Baranov (1957); Baranov & Naudy (1964); Fourier-domain form from Gunn (1975).   
*In SGTool*: original code. Filter values larger than 10<sup>6</sup> are set to 0; there is no other low-latitude stabilisation, so close to the magnetic equator use the Low-latitude RTP or RTE. (The equation shown here before 2026 was not the one the code uses.)   

**Reduction to the Equator**    
$`H_{RTE}(\theta) = \frac{-\cos^2(D-\theta)}{[\sin I + i \cos I \cos(D-\theta)]^2}`$     
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic equator (same declination); at I = 0 it leaves the data unchanged.
Where   
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit. 

*Reference*: Kis (1990); also Eq. 5 of Hao et al. (2018).   
*In SGTool*: original code, evaluated as $`-\cos^2(D-\theta)\,\bar f^2/|f|^4`$ with $`f = \sin I + i\cos I\cos(D-\theta)`$, so the denominator is real. Corrected in 2026: the earlier code did not square the denominator, so at I = 0 it multiplied the data by $`\cos^2(D-\theta)`$ instead of leaving them unchanged, and at other inclinations it gave anomalies that were too weak. The grid is padded by mirroring its edges, because the operator changes sharply near the direction perpendicular to the declination and would otherwise smear the streaks that ramp padding makes along the declination into the grid.   
   
**Low-latitude Reduction to the Pole (Taylor series iterative)**   
Stable RTP for low magnetic latitudes. The data are first reduced to the equator   
$`C(\theta) = \frac{-\cos^2(D-\theta)}{[i \cos I \cos(D-\theta) + \sin I]^2}`$   
and then to the pole with m iterations of a truncated (two-term Taylor series) equator-to-pole operator $`H_{ea} = -[1 + (D-\theta)^2/2]^2`$, which in closed form is   
$`\tilde T_p = \{1 - [1 - C_{pe} H_{ea}]^{m+1}\} H_e \tilde T_e`$, with $`H_e = -\sec^2(D-\theta)`$ and $`C_{pe} = -\cos^2(D-\theta)`$.   
More iterations get closer to the exact (unstable) RTP; the gain perpendicular to the declination is limited to about 5(m+1). The number of iterations is chosen automatically (0 to 200) as the one whose result correlates best with the analytic signal amplitude of the input grid. Assumes induced magnetisation only.   
*Iterations*: blank (0) chooses them automatically; fewer give less striping along the declination but a less complete reduction.   
*Stripe suppression (°)*: close to the equator the data hold almost no signal in the directions perpendicular to the declination, so what the operator amplifies there is mostly noise, edge effects and levelling errors, which show as stripes parallel to the declination. Within about this angle $`\sigma`$ of that direction the iterations fade out, $`m_{eff} = m\,[1 - e^{-(\phi/\sigma)^2}]`$, $`\phi`$ the angle of the wavenumber from the perpendicular. 0 turns it off (the method as published); 10° (the default) helps at any low inclination, 20° suits data within a few degrees of the equator, and wider values start to blur real features at higher inclinations. When the iterations are chosen automatically, they are chosen before the suppression is applied.   
The grid is padded by mirroring its edges, tapered to zero within the padding: padding that extends each column straight out from the north and south edges (as the other FFT filters do) makes streaks along the declination, which this operator turns into stripes across the grid.   
Output: `_RTP_LL`.   
*Reference*: Hao et al. (2018).   
*In SGTool*: original code (2026) implementing the paper's Eqs. 4, 5, 13 and 14. The closed form (Eq. 13) is used rather than iterating, evaluated so that it stays finite perpendicular to the declination. The automatic choice tries a fixed set of 21 values of m from 0 to 200 rather than every value. The stripe suppression and the mirror padding are SGTool additions, not part of the published method.   

**Variable (Differential) Reduction to the Pole**   
Accounts for the spatial variation of inclination and declination across a survey area by applying a Taylor-series expansion of the standard RTP operator. The IGRF field is computed at the four grid corners and the centre, then bilinearly interpolated to build spatially varying inclination and declination grids. The expansion is evaluated at 13 perturbed inc/dec pairs and combined to give a single corrected output.   
Outputs: `_DRTP` (corrected field), `_DRTP_inc` and `_DRTP_dec` (the interpolated inclination and declination fields used).   
Recommended for surveys covering several degrees of latitude/longitude where a single inc/dec value would introduce systematic phase or amplitude errors.   
*Reference*: Cooper & Cowan (2005).   
*In SGTool*: a line-by-line Python translation (Mark Jessell, 2026) of Gordon Cooper's MATLAB `rtpvariable.m`, supplied by its author. SGTool adds the IGRF calculation of the corner and centre values and the filling and restoring of no-data cells.   

**Continuation**    
$`H(k) = e^{-k h}`$ (upward), $`H(k) = e^{k h}`$ (downward)   
Where   
h > 0 is the continuation distance, in the grid's length units.   
*Reference*: Henderson & Zietz (1949); Fourier-domain form from Bhattacharyya (1965) and Gunn (1975).   
*In SGTool*: original code. Downward continuation is not stabilised, so it amplifies short-wavelength noise very quickly; keep h small compared with the depth to the sources. For a geographic grid the height is converted roughly to degrees.   
   
**Vertical Integration**   
$`H(k_x, k_y) = \frac{1}{k}`$  
When applied to an RTE or RTP image provides the so called Pseudogravity result    
Where    
k = $\sqrt{k_x^2 + k_y^2}$ .   
*Reference*: pseudogravity, Baranov (1957); the integration operator, Gunn (1975).   
*In SGTool*: original code. k is not allowed below 10<sup>-4</sup> rad/m, so the longest wavelengths stay finite. The constants of Poisson's relation (G and the ratio of magnetisation to density) are not applied, so the result is proportional to the pseudogravity, not equal to it. Needs a grid in a metre-based projection.   
   

## Frequency Filters   

**What the cutoff values mean**   
The values entered for the High pass, Low pass and Band pass filters are **wavelengths**, in the length units of the grid (metres, or degrees for a geographic grid), not frequencies. A wavelength is the distance over which a feature repeats, so it is the inverse of frequency: *short* wavelength = *high* frequency (small, sharp, usually shallow features) and *long* wavelength = *low* frequency (broad, smooth features). The equations below use the wavenumber $k = 2\pi/\lambda$.   

- **High pass** keeps wavelengths *shorter* than the cutoff and removes longer ones (regionals, deeper sources).   
- **Low pass** keeps wavelengths *longer* than the cutoff and removes shorter ones (noise, shallow sources).   
- **Band pass** keeps the wavelengths between two values. **Low cut wavelength** is the *longest* wavelength kept (longer wavelengths, i.e. lower frequencies, are cut) and **High cut wavelength** is the *shortest* wavelength kept (shorter wavelengths, i.e. higher frequencies, are cut), so the Low cut value is the larger number (for example 50000 and 5000).   
- **Transition width** (same units) is how gradually the filter switches from pass to reject: 0 is a sharp cut, and larger values reduce ringing.   

*Relationship to source depth (approximate).* The signal from a source at depth $z$ falls off with wavelength $\lambda$ as $e^{-2\pi z/\lambda}$, the same factor that governs upward continuation, so a deep source cannot produce short wavelengths. A source at depth $z$ typically has its main response at wavelengths of about 3 to 6 times $z$. As a rough guide, a cutoff of 6000 m separates sources shallower than about 1000 to 2000 m (kept by a high pass) from deeper ones (kept by a low pass), and a band pass between $\lambda_{high}$ and $\lambda_{low}$ favours sources from about $\lambda_{high}/6$ to $\lambda_{low}/3$. This is only a guide: a shallow, broad source also has long wavelengths, and potential-field data are inherently ambiguous, so filtering by wavelength is not a true depth slice.   

**High Pass Filter**

```math
H(k) = \begin{cases}
0 & \text{if } k \leq k_{low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{low}}{k_{high} - k_{low}}\right)\right) & \text{if } k_{low} < k < k_{high} \\
1 & \text{if } k \geq k_{high}
\end{cases}
```

The high-pass filter removes low-frequency components (long wavelengths) while retaining high-frequency components (short wavelengths) with a smooth transition to reduce ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{low} = \frac{2\pi}{\lambda_c + w/2}$ : Lower transition boundary
- $k_{high} = \frac{2\pi}{\lambda_c - w/2}$ : Upper transition boundary  
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (in same units as wavelength)

**Low Pass Filter**

```math
H(k) = \begin{cases}
1 & \text{if } k \leq k_{inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{inner}}{k_{outer} - k_{inner}}\right)\right) & \text{if } k_{inner} < k < k_{outer} \\
0 & \text{if } k \geq k_{outer}
\end{cases}
```

The low-pass filter removes high-frequency components (short wavelengths) while retaining low-frequency components (long wavelengths). Optional smooth transition reduces potential ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{inner} = \frac{2\pi}{\lambda_c + w/2}$ : Inner transition boundary
- $k_{outer} = \frac{2\pi}{\lambda_c - w/2}$ : Outer transition boundary
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (optional, in same units as wavelength)

**Band Pass Filter**

```math
H(k) = H_{high}(k) \times H_{low}(k)
```

Where:

***High-pass component:***

```math
H_{high}(k) = \begin{cases}
0 & \text{if } k \leq k_{h,low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{h,low}}{k_{h,high} - k_{h,low}}\right)\right) & \text{if } k_{h,low} < k < k_{h,high} \\
1 & \text{if } k \geq k_{h,high}
\end{cases}
```

***Low-pass component:***

```math
H_{low}(k) = \begin{cases}
1 & \text{if } k \leq k_{l,inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{l,inner}}{k_{l,outer} - k_{l,inner}}\right)\right) & \text{if } k_{l,inner} < k < k_{l,outer} \\
0 & \text{if } k \geq k_{l,outer}
\end{cases}
```

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

*Reference* (high, low and band pass): ideal filters with a raised-cosine (Hann) transition band, a standard design with no single primary reference.   
*In SGTool*: original code. With a transition width of 0 the low pass is a sharp cut; when $\lambda_c - w/2 \le 0$ the outer edge of the transition is set to 1.5 times the cutoff wavenumber instead.   

**Directional line-noise removal (Directional Band Pass)**   
Removes noise parallel to the flight lines. A zero-mean noise estimate *N* is made by passing the grid through a Butterworth band pass and a directional wedge, and is subtracted after scaling: $`T' = T - s\,N`$.   
$`H(k, \theta) = \frac{1}{1 + (k_h/k)^{2n}} \cdot \frac{1}{1 + (k/k_l)^{2n}} \cdot W(\theta)`$   
$`W(\theta) = \cos^2\!\left(\frac{\pi\,\Delta}{2w}\right)`$ for $`\Delta \le w`$, 0 beyond   
Where   
- $k_h = 2\pi/(10 \times \text{largest line spacing})$ and $k_l = 2\pi/(2 \times \text{smallest line spacing})$: the band kept, so that the regional, which also varies across the lines, is not taken as noise.   
- n = 4 : Butterworth order.   
- $\Delta$ : angle between the wavenumber and the azimuth, measured modulo 180°, so that both lobes of the spectrum are treated.   
- w : wedge half-width (default 45°).   
- s : Scale (default 1).   

*Reference*: Butterworth (1930) for the filter response; removing line noise with a directional filter follows the micro-levelling idea of Minty (1991).   
*In SGTool*: original code. The band limits (2 to 10 times the line spacing), the cos² wedge and the scale are SGTool's choices, and unlike Minty (1991) the noise estimate is not amplitude-limited before it is subtracted. (The Butterworth high-pass and |cos|<sup>p</sup> directional filter described here before 2026 are not the filters the dialog uses.)   

**Remove Regional**   
Remove a 1st order (dipping plane) or 2nd order (parabolic plane) regional from data.   
*Reference*: least-squares trend-surface fitting, a standard method with no single primary reference.   
*In SGTool*: original code: the plane or quadratic surface is fitted by least squares to the valid cells and subtracted.   
   
**Automatic Gain Control**    
$`AGC(x, y) = \frac{f(x, y)}{G_w * |f|(x, y)}`$   
Where    
$`G_w * |f|`$ is the absolute value of the data smoothed with a Gaussian of standard deviation w (the window size, in cells).   
*Reference*: Rajagopalan & Milligan (1994).   
*In SGTool*: original code. The local amplitude is the Gaussian-smoothed absolute value, not a windowed RMS (the RMS form shown here before 2026 is not what the code does), the edges are reflected, and no other scaling is applied.   
   
**Radially averaged power spectrum (but needs testing!)**    
$`P(k) = \frac{1}{N_k} \sum_{(k_x, k_y) \in k} |\text{FFT}(f)|^2`$   
Where    
P(k) is the radially averaged power spectrum, and N<sub>k</sub> is the number of samples in the radial bin.   
*Reference*: Spector & Grant (1970), for using the radially averaged spectrum to estimate source depths.   
*In SGTool*: original code. No-data cells are set to 0 before the transform, ln P is plotted against wavenumber, and no depths are estimated.   
   

## Gradient Filters   
   
**Derivative**    
$`H_x = (i k_x)^n, \quad H_y = (i k_y)^n, \quad H_z = k^n`$   
Where   
n is the power of the derivative (it can be fractional).   
The vertical derivative is calculated from the horizontal gradients using Laplace's equation, $\frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} + \frac{\partial^2 f}{\partial z^2} = 0$, which gravity and magnetic fields satisfy outside their sources (they are Green's functions). In the Fourier domain this multiplies the spectrum by $k$, so the vertical derivative amplifies short wavelengths and often amplifies noise as well.   
*Reference*: Bhattacharyya (1965); Gunn (1975).   
*In SGTool*: original code; the result is zero-mean. (The rotated-gradient formula shown here before 2026 is not what the code does: the derivative is taken along x, y or z only.)   
   
**Total Horizontal Gradient**   
$`THG(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}`$   
*Reference*: Cordell & Grauch (1985).   
*In SGTool*: original code, with the horizontal derivatives calculated in the Fourier domain.   
   
**Analytic Signal**    
$`A(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2 + \left(\frac{\partial f}{\partial z}\right)^2}`$   
Computes the total amplitude of the gradients, independent of field inclination or declination (exactly so only for 2D sources; approximately in 3D).
Useful for locating edges of potential field sources (e.g., faults or contacts).   
*Reference*: Nabighian (1972) for 2D sources; Roest, Verhoef & Pilkington (1992) for the 3D amplitude.   
*In SGTool*: original code, with the derivatives calculated in the Fourier domain.   
       
**Tilt Angle**    
$`T = \tan^{-1}\left(\frac{\frac{\partial f}{\partial z}}{\sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}}\right)`$   
Enhances the contrast of geological features by highlighting gradients relative to the vertical component. The result is saved in degrees (-90 to 90).
Where   
df/dz : Vertical derivative of the field.
df/dx , df/dy : Horizontal derivatives of the field.   
*Reference*: Miller & Singh (1994).   
*In SGTool*: original code, with the derivatives calculated in the Fourier domain.   
   
## References   

- Baranov, V., 1957. A new method for interpretation of aeromagnetic maps: pseudo-gravimetric anomalies. *Geophysics* 22, 359–382. https://doi.org/10.1190/1.1438369   
- Baranov, V., Naudy, H., 1964. Numerical calculation of the formula of reduction to the magnetic pole. *Geophysics* 29, 67–79. https://doi.org/10.1190/1.1439334   
- Bhattacharyya, B.K., 1965. Two-dimensional harmonic analysis as a tool for magnetic interpretation. *Geophysics* 30, 829–857. https://doi.org/10.1190/1.1439658   
- Butterworth, S., 1930. On the theory of filter amplifiers. *Wireless Engineer* 7, 536–541.   
- Cooper, G.R.J., Cowan, D.R., 2005. Differential reduction to the pole. *Computers & Geosciences* 31, 989–999. https://doi.org/10.1016/j.cageo.2005.02.005   
- Cordell, L., Grauch, V.J.S., 1985. Mapping basement magnetization zones from aeromagnetic data in the San Juan basin, New Mexico. In: Hinze, W.J. (Ed.), *The Utility of Regional Gravity and Magnetic Anomaly Maps*. Society of Exploration Geophysicists, 181–197.   
- Gunn, P.J., 1975. Linear transformations of gravity and magnetic fields. *Geophysical Prospecting* 23, 300–312. https://doi.org/10.1111/j.1365-2478.1975.tb01530.x   
- Hao, M., Zhang, F., Tai, Z., Du, W., Ren, L., Li, Y., 2018. Reduction to the pole at low latitudes by using the Taylor series iterative method. *Journal of Applied Geophysics* 159, 127–134. https://doi.org/10.1016/j.jappgeo.2018.08.004   
- Henderson, R.G., Zietz, I., 1949. The upward continuation of anomalies in total magnetic intensity fields. *Geophysics* 14, 517–534. https://doi.org/10.1190/1.1437560   
- Kis, K.I., 1990. Transfer properties of the reduction of magnetic anomalies to the pole and to the equator. *Geophysics* 55, 1141–1147. https://doi.org/10.1190/1.1442930   
- Miller, H.G., Singh, V., 1994. Potential field tilt — a new concept for location of potential field sources. *Journal of Applied Geophysics* 32, 213–217. https://doi.org/10.1016/0926-9851(94)90022-1   
- Minty, B.R.S., 1991. Simple micro-levelling for aeromagnetic data. *Exploration Geophysics* 22, 591–592. https://doi.org/10.1071/EG991591   
- Nabighian, M.N., 1972. The analytic signal of two-dimensional magnetic bodies with polygonal cross-section: its properties and use for automated anomaly interpretation. *Geophysics* 37, 507–517. https://doi.org/10.1190/1.1440276   
- Rajagopalan, S., Milligan, P., 1994. Image enhancement of aeromagnetic data using automatic gain control. *Exploration Geophysics* 25, 173–178. https://doi.org/10.1071/EG994173   
- Roest, W.R., Verhoef, J., Pilkington, M., 1992. Magnetic interpretation using the 3-D analytic signal. *Geophysics* 57, 116–125. https://doi.org/10.1190/1.1443174   
- Spector, A., Grant, F.S., 1970. Statistical models for interpreting aeromagnetic data. *Geophysics* 35, 293–302. https://doi.org/10.1190/1.1440092   
