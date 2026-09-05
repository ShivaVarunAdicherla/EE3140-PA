#set page(paper: "a4", margin: (top: 1.5cm, left: 1.2cm, right: 1.2cm, bottom: 1cm))
#show underline: set underline(offset: 2pt)
#show link: set text(fill: blue)
#show link: underline
#show figure.caption: it => it.body
#set par(first-line-indent: 1.4em)
#align(center)[
  #show heading.where(level: 1): set text(size: 20pt)
  #show heading.where(level: 2): set text(size: 18pt)
  #set text(weight: "semibold", size: 14pt)
  = Waveguides and Antennas [EE3140]
  == Programming Assignment - 1
  #v(0.4em)
  #grid(
    column-gutter: 2em,
    columns: 2,
    [*EE24BTECH11201*], [*Shiva Varun Adicherla*],
  )
  #v(-0.4em)
  #line(length: 100%)
]

#set heading(numbering: "A.1")
#let pagebridge = align(right)[
  #v(1.2em)
  _(Continued in the next page)_
  #pagebreak()
]
#let to-sci(val, digits: 2) = {
  if val == 0 { return $0$ }
  let exponent = calc.floor(calc.log(calc.abs(val)))
  let mantissa = calc.round(val / calc.pow(10, exponent), digits: digits)
  $#mantissa times 10^(#exponent)$
}
#v(-1.1em)
= Part
==
#let dt = 5E-13
#let T = 15E-9
The spatial division used in the simulation is given by:
$ Delta z = "Domain Size"/(N-1) => (6 m)/2000 approx 0.3 "cms" $
This is used because between N points we'll have "N-1" $Delta z$, which we are assuming
will span the whole domain. The *Odd* number of points were chosen because it allows for a symmetric
division of the domain into 2 equal halves.

The simulation time step I've used is *$d t = #to-sci(dt) s$*, along with simulating for *$15$ns*, because,
if we look at the source, it's basically a guassian so we start with the exponential at $-(5/2)^2$
and we'll land back to the same at *$t=10"ns"$*, so dividing the time I want to simulate with the time step yields:
$ "Number of Steps" = T/(d t) => #to-sci(T / dt) $

Now for the spatial resolution, which I've chosen to interpret in the following way:
$ "Spatial Resolution" = lambda/(Delta z) $
$lambda$ isn't hard to calculate and is given by $lambda = c_0 times f^(-1) => 30$ cms, and out spatial resolution will be given by:
$ "Spatial Resolution" = 30/0.3 approx 100 $

The courant condition from my understanding is essentially a part of establishing *causality* in the simualtion,
and we achieve this by making sure that the information flow only happens to the adjacent cells/points.
$ max(V_"in the simulation") times d t <= Delta z ==> d t <= (Delta z)/max(V_"in the sim") $
Since this is free space, we can approximate the speed of light $c_0 approx #to-sci(3E8)$ m/s, so our $d t$'s limit will be:
$ d t <= (3 times 10^(-3))/(3 times 10^8) => d t <= 10"ps" $
Which I've considered and satisfied, by taking $d t = #to-sci(5E-13)$ .

==
I've animated the evolution of electric and magnetic fields in gnuplot and the gifs and the source code
are in my #link("https://github.com/ShivaVarunAdicherla/EE3140-PA", "Github repo"), but some of the images are attached below.
#grid(
  columns: 2,
  column-gutter: 0.3em,
  image("./PartA/plot1.pdf"), image("./PartA/plot2.pdf"),
)

==
To estimate the propagation velocity, I've essentially fixed a level whose position I'm going to check in the next instant.
Illustrated by the following figure.
#figure(
  image("Velocity Illustration.png", width: 50%),
  caption: [The red line represents the level I fixed, the solid line the current E, dashed the future E],
)

A simple warning though, this approach requires us to either find the exact values of the intersections which
we can guarentee by *Linearly Interpolating* between the points in which we expect the intersection to be. And
this too only is *accurate at the zero crossing points* because the slope is maximum there. Having a low slope
at the point of interest will cause issues with numerical precision because the target value is already close.

The velocity then is given by the following:
$ V = (delta z)/(d t) $
#align(center, rect(
  ```c
    double target = g->Egrid[i];
                  ... // This is not the complete code, just an outline
      double e1 = g->Egrid_prev[j];
      double e2 = g->Egrid_prev[j + 1];

      // Find the region in which the target lives
      if ((e1 <= target && target <= e2) || (e2 <= target && target <= e1)) {
        double alpha = 0.0;
        double diff = e2 - e1;
        if (fabs(diff) > 1e-13) { // Failsafe if Slope is really small
          alpha = (target - e1) / diff;
        }
                  ...
        }
      }
    }
                  ...
    double dz = ((double)i - alpha) * g->dx;
    return dz / dt;
  }
  ```,
))
After choosing a proper point, and using this method, I've calculated the velocity to be:
$ V => c_0 = 299757617.346131 "m/s" (delta = 0.0116 %) $

Determining the Intrinsic Impedance is straight forward and can be evaluated anywhere except for the
points at which H is very low, and lead to numerical errors, it is given by:
$ eta = E/H $
#align(center, rect(outset: 0.2em)[
  ```c
  printf("Calculated Impedance is :%lf Ohm\n",
          fabs(2 * g.Egrid[1620] / (g.Hgrid[1620] + g.Hgrid[1619])));

  ```])

Since we have the #underline[H and E field's grid staggered according to the Yee's Method, we'll have to account for this
  by avergaing out the values of H adjacent to the E grid point], Which I'm going in the code above.

Here, I've chosen a point where the Electric and Magnetic fields were at the maximum, to avoid the said numerical
error.

$ eta = E/H => 376.7605" "Omega " "(delta = 0.605%) $

==
The electric field undergoes throught the polarity inversion because of the boundary condtions imposed by
the maxwell's equations, which essentially state that the #underline[tangential electric field at the surface
  of PEC boundary must be 0, so the reflected Electric field is forced to have a polarity opposite to that of the incident Electric Field.]
As a result of this, we'll also have some surface current on the PEC boundary which also forces the Magnetic field inside the PEC
to be 0 too.

#text("Program output", size: 13pt, weight: "bold")
```sh
❯ ls
- makefile  parta  parta.c  plot1.pdf  plot2.pdf  plot.gif
❯ ./parta
- calculated velocity is :299757617.346131 m/s
- calculated impedance is :376.760500 ohm
```

= Part
#v(0.4em)
For this part of the question, to avoid additional reflections from the source (The source I've used in the
previous part is a hard source), I've changed the source type here to a *soft source*, since there isn't a constraint mentioned in the
question sheet.
==

Here the identification of Incident and Transimitted waves is straight forward, but the identification
of reflected wave is not as straight forward. As given in the note of the question, we'll have the reference
incident wave from the simulation we have had in Part A or the wave that was generated in the other direction as the reference
wave in the free space. #underline(offset: 2pt)[Subtracting the complete wave we have in this
  simulation with the reference wave we had, will essentially yeild *us the pure reflected wave for a brief Interval.*]

The pictures are attached below:
#align(center + horizon, image("./PartB/SeperateWave-NW.pdf", height: 27%))
#align(center, image("./PartB/SeperateWave-TI.pdf", height: 27%))
#align(center + horizon, image("./PartB/SeperateWave-TR.pdf", height: 27%))

The figure labelled net wave is the RAW simulation data, and other $2$ figures contain the #underline[transmitted wave in the
  RED coloured region (Basically marks the Dielectric with *$epsilon_r = 3.2$*) and the Reflected/Transmitted part in the White/Nocolor Region according
  to the figure label.]

#align(center, block(inset: 0.4em, stroke: 1pt)[
  ```c
     if (i == 23000) { // Timestep at which I'm capturing the snapshot
      double tempE[g.N];
      double tempH[g.N];
      // Reflected Wave
      for (int i = 1; i < 2001; i++) {
        // E field is flipped since I'm using the wave that was
        // propagating in the other direction.
        tempE[3000 + i] = g.Egrid[3000 + i] - g.Egrid[3000 - i];
        tempH[3000 + i] = g.Hgrid[3000 + i] + g.Hgrid[3000 - i];
      }
      for (int i = 5001; i < 6001; i++) {
        tempE[i] = g.Egrid[i];
        tempH[i] = g.Hgrid[i];
      }
      gplotSetOption(&g1, "set output \"SeperateWave-TR.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Transmitted and Reflected Wave\"");
      gplotPBDSEM(&g1, tempE, tempH, "E", "H", g.N);
      // Incident wave
      for (int i = 3001; i < 5000; i++) {
        tempE[i] = g.Egrid[6000 - i];
        tempH[i] = -g.Hgrid[6000 - i];
      }
      gplotSetOption(&g1, "set output \"SeperateWave-TI.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Transmitted and Incident Wave\"");
      gplotPBDSEM(&g1, tempE, tempH, "E", "H", g.N);
      // Complete wave
      gplotSetOption(&g1, "set output \"SeperateWave-NW.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Net Wave\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "E", "H", g.N);
    }
  ```
])

==
Similar to the question in Part A, we can extract the value of the velocity using the same method, which yeilds us:
$ V = 167516463.117863 "m/s" (delta = 0.04%) $
$ "We can verify this with" V = 1/(sqrt(epsilon mu)) => (c_0)/(sqrt(3.2)) = 167589078.80743 "m/s" $

==

For calculating the Transmission and Reflection co-efficients, I've basically taken time instants at which
the target part *of the wave is isolated*,#underline[ then identified the peak/Extrema value of Electric field over the target range
  which I've found out by inspection.] I've chosen the time instants depicted in the following graphs/figures.

#grid(
  columns: 2,
  figure(image("./PartB/plot1B.pdf")), figure(image("./PartB/plot2B.pdf")),
)

The figure on the left is/was used to extract the peak value of Electric Field in Incident wave and Reflected Wave,
where as the figure in the right was used to extract the peak Value of Electric Field in the dielectric.
#underline[The main idea is that the *peak stays at the same spot with respect to the whole wave*.] Then dividing the
Electric Field Amplitudes will yeild us the required Co-efficients as they are defined as such.
$
  Gamma = E_"Reflected"/E_"Incident" => (eta_2 -eta_1)/(eta_2 + eta_1) "     " Tau = E_"Transmitted"/E_"Incident" => (2eta_2)/(eta_2+eta_1)
$
Calculating the values of $Gamma$ and $Tau$ from the given information, we get:
$
  Gamma = ((1/sqrt(3.2) - 1)/(1/sqrt(3.2)+1)) => #(calc.round((1 / calc.sqrt(3.2) - 1) / (1 / calc.sqrt(3.2) + 1), digits: 6))
  "     " Tau = ((2/sqrt(3.2))/(1/sqrt(3.2)+1)) => #(calc.round((2 / calc.sqrt(3.2)) / (1 / calc.sqrt(3.2) + 1), digits: 6))
$

Extracted values of $Gamma$ and $T$ from field amplitude samples are:
$ Gamma = -0.282924" "(delta = 0.022%) "    " Tau = 0.716769" "(delta = 0.051%) $

The sign extraction of the co-efficients was possible because we had asymetric waveform with regards to the
maximum/minimum values. I've exploited this fact to extract the sign, and the code responsible and it's explaination is given in the next page.
#pagebridge
#align(center, block(stroke: 1pt, inset: 0.4em)[
  ```c
  if (i == 16000) // TimeStep at which I'm simpling the Peaks
      incidentMax = (findMax(g.Egrid, g.Egrid + g.N - 1) >
                     fabs(findMin(g.Egrid, g.Egrid + g.N - 1))) // If the +ve peak is greater
                        ? findMax(g.Egrid, g.Egrid + g.N - 1) // Take the amplitude as +ve Max
                        : findMin(g.Egrid, g.Egrid + g.N - 1); // Else we take the -ve Peak
  if (i == 33333) {
    printf("The velocity of the EM Wave in the Dielectric is: %lf\n",
           returnVelocity(&g, 5020, dt, 2));
    reflectedMax = (findMax(g.Egrid + 3000, g.Egrid + 4000) >
                    fabs(findMin(g.Egrid + 3000, g.Egrid + 4000)))
                       ? findMax(g.Egrid + 3000, g.Egrid + 4000)
                       : findMin(g.Egrid + 3000, g.Egrid + 4000);
  }
  if (i == 39999) {
    transmittedMax = (findMax(g.Egrid + 5001, g.Egrid + 6000) >
                      fabs(findMin(g.Egrid + 5001, g.Egrid + 6000)))
                         ? findMax(g.Egrid + 5001, g.Egrid + 6000)
                         : findMin(g.Egrid + 5001, g.Egrid + 6000);
  }
  ```
])
==
#text("Assumption:", size: 13pt, weight: "semibold") Since the question doesn't mention if the poynting vector plot should be with respect to the time
or with respect to the $z$, I'm going to consider the plot to be along time.

The calculation of poynting vector is straight forward, with us *just taking $E times H$*
and since $E$ is along $arrow(j)$ and $H$ along $arrow(x)$, we'll have to compensate with a
*-ve sign.* (Since $arrow(j) times arrow(i) = -arrow(k)$)

#figure(
  image("./PartB/PoyntingVec.pdf", width: 75%),
)

Identification of power flow inside the simulation is easy - #underline[we can just look at the sign of poynting vector,
  If the sign is +ve, the power flow is towards $+arrow(z)$ and when it's -ve it's towards $-arrow(z)$].
\
We can prove/find out the average power flow through integrating the poynting flux at a point,
which we can give by just doing the following:
$ sum_tau E_k [tau] = sum_tau E_(l)[tau]"     " "Where" k != l $


#align(center, block(stroke: 1pt, inset: 0.4em)[
  ```c
  // Book-keeping for Part 4
  poyntingVec1[i] = -(g.Egrid[4996]) * (g.Hgrid[4996] + g.Hgrid[4997]) / 2;
  // printf("%lf\n", poyntingVec1[i]);
  poyntingInt1 -= (g.Egrid[4996]) * (g.Hgrid[4996] + g.Hgrid[4997]) / 2;
  poyntingInt2 -= (g.Egrid[5004]) * (g.Hgrid[5004] + g.Hgrid[5005]) / 2;
  ```
])

So essentially, #underline[if the values of Time integral of poynting vector at 2 (different) points is the
  same], we'll have proved that Energy conservation indeed holds true at the interface between free-space and
dielectric.

#text("Program output", size: 13pt, weight: "bold")
```sh
❯ ls
- Makefile  partB  partB.c  plot1B.pdf  plot2B.pdf  plot.gif  PoyntingVec.pdf  SeperateWave-NW.pdf  SeperateWave.pdf  SeperateWave-TI.pdf  SeperateWave-TR.pdf
❯ ./partB
- The velocity of the EM Wave in the Dielectric is: 167516463.117863
- Reflection Co-efficient is: -0.282924
- Transmission Co-efficient is: 0.716769
- Time Integral of Poynting Vector before and after the dielectric boundary: 68.069187 68.044646
```

And we can indeed see that the *time average of poynting vector basically yeilded the approximately the same values,
hence confirming that the interface is lossless.*

= Part
#v(0.4em)
I really don't need to extend the simulation framework, since we have projected the wave into orthogonal
basis, so I can instead run 2 simulations independently (Maxwell's equations are essentially linear here)
and combine the results.

==

I've chosen an observation point around $z = 0.6m$, and used $2$ 1D simulations, and combined the results to get the
resulting polarisation trajectory. Since we can't have anything that's interactive in this pdf/report, I'm omitting
the time evolution of the polarisation graph. It is however available on my git repository, that's given at the end of this
report.

==

#grid(
  columns: 4,
  column-gutter: 0.4em,
  row-gutter: 1em,
  figure(image("./PartC/PlotPhi0.000.pdf"), caption: [$phi = 0" rad"$]),
  figure(image("./PartC/PlotPhi0.785.pdf"), caption: [$phi = pi/4" rad"$]),
  figure(image("./PartC/PlotPhi1.571.pdf"), caption: [$phi = pi/2" rad"$]),
  figure(image("./PartC/PlotPhi2.356.pdf"), caption: [$phi = (3pi)/4" rad"$]),

  figure(image("./PartC/PlotPhi3.142.pdf"), caption: [$phi = pi" rad"$]),
  figure(image("./PartC/PlotPhi3.927.pdf"), caption: [$phi = (5pi)/4" rad"$]),
  figure(image("./PartC/PlotPhi4.712.pdf"), caption: [$phi = (3pi)/2" rad"$]),
  figure(image("./PartC/PlotPhi5.498.pdf"), caption: [$phi = (7pi)/4" rad"$]),
)

I obtained the above figures after sweeping the $phi$ over $[0,2pi)$ in steps of
$pi/4$. The perspective of the above plots is of the view when the wave is moving
out of the paper/screen.

==

As discussed in the class, we can describe this with the following equation.
$ (x/(E_x))^2 + (y/E_y)^2 - (2 x y)/(E_x E_y) = sin^2(phi) $
Below are the figures which vary both $phi$ and Amplitudes of the waves.
#grid(
  columns: 2,
  column-gutter: 1em,
  figure(image("./PartC/PlotPhi0.785-A2.pdf"), caption: [$phi = pi/4" rad"$]),
  figure(image("./PartC/PlotPhi1.571-A2.pdf"), caption: [$phi = pi/2" rad"$]),
)
The convention used to describe the direction of polarisation is also discussed in the class and is given #underline(offset: 2pt)[by the hand
  which is able to twist in the direction of rotation polarisation while the thumb of the said hand is in the
  direction of propagation of the EM wave].

#v(1em)
#par(first-line-indent: 0em)[
  #text("GitHub Repository", weight: "bold", size: 16pt)
]
#v(0.4em)
The github repository can be accessed #link("https://github.com/ShivaVarunAdicherla/EE3140-PA", "here").
It will include a readme for this Programming assignment holding all the animated Gifs, and
the requested source code.

#v(1em)
#line(length: 100%)
#align(center)[
  #text("END", weight: "bold", size: 14pt)
]
