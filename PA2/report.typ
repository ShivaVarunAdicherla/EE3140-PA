#set page(paper: "a4", margin: (
  top: 1.5cm,
  left: 1.2cm,
  right: 1.2cm,
  bottom: 1cm,
))
#show underline: set underline(offset: 2pt)
#show link: set text(fill: blue)
#show link: underline
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

#set heading(numbering: "1.A")
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
// The 1st Question
= Part

#pagebreak()
// The 2nd Question
= Part

#pagebreak()
// The 3rd Question
= Part
==
As discussed in the class, if we assume that $omega$ is high enough to satisfy the conditions: $R << omega L$ and $G << omega C$, then we can consider that transmission line, a *_low loss_* transmission
line case, and the following quantities can be calculated by:

$
  Z_0 = sqrt((R+j omega L)/(G + j omega C)) approx sqrt(L/C) & ==> Z_0 = sqrt((250 times 10^(-9))/(100 times 10^(-12))) = 50 Omega \
  v_p = 1/sqrt(L C) & ==> v_p = 1/sqrt(25 times 10^(-18)) = 2 times 10^8 m s^(-1)
$

The above phase/propagation velocity is physical and valid, since it's lower than the intrinsic velocity of EM waves in Free Space / Vacuum ($c_0 approx 3 times 10^8 m s^(-1)$).

==
To determine the values of Reflection co-efficient $Gamma$ at the given 2 points, we:
$
              Gamma & = (Z_L - Z_0)/(Z_L + Z_0) \
    Gamma_"At Load" & = (100 - 50)/(100+50) ==> 1/3 \
    "In the Switch" & "is closed condition" \
  Gamma_"At Source" & = (25 - 50)/(25 + 50) ==> -1/3 \
    "In the Switch" & "is open condition" \
  Gamma_"At Source" & = (infinity - 50)/(infinity + 50) ==> 1 \
$

Some other important quantities that we are going to need are:
$
                   T_("Single Way") = d/v_p & => 20/(2 times 10^8) = 100 n s \
  V_("Initial") = 1 times Z_0 / (Z_0 + R_s) & => 50 / 75 = 2/3 V \
$

Also according to the question, we have the switch closed at $t=0s$ (Inferred it's was open before this), and then opens again at $t=250n s$.

In the case where we were calculating the $Gamma$ at the load and source, we are essentially seeing an impedance dis-continuity, and had to satisfy the general voltage and current equations
$V(z) = V_0^+ exp(-j beta z) + V_0^- exp(j beta z)$ and $I(z) = I_0^+ exp(-j beta z) - I_0^- exp(j beta z)$ continuity at the impedance mis-matched boundaries/z locations, along with satisfying $V(z) =
Z(z) I(z)$. When the switch opens, we'll have the current at the location of switch to be zero. To calculate the extra voltage pulse we get when opening the switch,
we solve the boundary condition in the following way.

$
  I(0, 250 n s) = 2/(3 Z_0) -2/(9 Z_0) -2/(27 Z_0) + V_2^+/Z_0 = 0 #h(3em) ("The 2/9 Voltage pulse is moving in opposite direction")
$

Then the Value of $V_2^+ = -10/27$

== Bounce Diagrams and Waveforms
#figure(
  image(
    "BoundsDiagram.png",
    width: 100%,
  ),
  caption: [The bounds diagram of the given scenario],
)

#import "@preview/simple-plot:1.1.0": line-plot, plot

#figure(
  plot(
    xmin: 0,
    xmax: 750,
    ymin: -0.2,
    ymax: 1.2,
    width: 17cm,
    height: 3.5cm,
    show-grid: true,
    xtick-step: 100,
    xlabel: [Time (ns)],
    ylabel: [Voltage (V)],

    line-plot(
      (
        (0, 0),
        (100, 0),
        (100, 0.899),
        (300, 0.899),
        (300, 0.790),
        (350, 0.790),
        (350, 0.296),
        (500, 0.296),
        (500, 0.263),
        (550, 0.263),
        (550, 0.098),
        (700, 0.098),
        (700, 0.087),
        (750, 0.087),
      ),
      stroke: blue + 1.2pt,
    ),
  ),
  caption: [The $V(t)$ observed at the Load],
)

#figure(
  plot(
    xmin: 0,
    xmax: 750,
    ymin: -0.2,
    ymax: 1.2,
    width: 17cm,
    height: 3.5cm,
    show-grid: true,
    xtick-step: 100,
    xlabel: [Time (ns)],
    ylabel: [Voltage (V)],

    line-plot(
      (
        (0, 0),
        (0, 0.6666),
        (200, 0.6666),
        (200, 0.814),
        (250, 0.814),
        (250, 0.444),
        (400, 0.444),
        (400, 0.395),
        (450, 0.395),
        (450, 0.148),
        (600, 0.148),
        (600, 0.131),
        (650, 0.131),
        (650, 0.049),
        (750, 0.049),
      ),
      stroke: blue + 1.2pt,
    ),
  ),
  caption: [The $V(t)$ observed at the Source],
)

The desmos versions of these graphs are available #link("https://www.desmos.com/calculator/nqwhf1hsvi", "at"), and we can also interactively view the voltages. However since we'll have
infinite reflections in this case, they've been truncated to about 5 in the calculations.

== FTDT Simulation
For this part, I've hacked/modified my previous FTDT simulator used in the previous assignment. All the changes have been documented in the source code as comments.
=== Lossless condition
The gist of this FTDT simulation is:
$
  (partial V)/(partial z) & = -L (partial I)/(partial t) - R I #h(3em) ("R is kept 0 in lossless condition") \
  (partial I)/(partial z) & = -C (partial V)/(partial t) - G V #h(3em) ("G is kept 0 in lossless condition") \
$
Which is solved by:
#pagebridge
#align(center)[
  #rect(
    ```c
    cCon = (dt / (grid->dx * Inductance));
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] -= cCon * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    cCon = (dt / (grid->dx * Capacitance));
    for (int j = 0; j < grid->N; j++) {
      grid->Egrid[j] -= cCon * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]

#figure(
  image(
    "PartC/Lossless_Voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

#let rows = csv("PartC/lossless.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.47,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}


=== Lossy condition
#v(0.5em)
#align(center)[
  #rect(
    ```c
    cCon = (dt / (2.0 * Inductance));
    double cCon1 = (1.0 - Resistance * cCon) / (1.0 + Resistance * cCon);
    double cCon2 = (dt / (Inductance * grid->dx)) / (1.0 + Resistance * cCon);
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] =
          cCon1 * grid->Hgrid[j] - cCon2 * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    cCon = (dt / (2.0 * Capacitance));
    cCon1 = (1.0 - Conductance * cCon) / (1.0 + Conductance * cCon);
    cCon2 = (dt / (Capacitance * grid->dx)) / (1.0 + Conductance * cCon);
    for (int j = 1; j < grid->N - 1; j++) {
      grid->Egrid[j] =
          cCon1 * grid->Egrid[j] - cCon2 * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]

#figure(
  image(
    "PartC/Lossy_voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

The Theoretical data in the following plots is still the lossless case, and is just there for Comparison. And as expected the peak/max voltages that we see in the simulation are lower than expected,
while keeping the exact same shape. This essentially means the transmission line is still dispersion-less.

#let rows = csv("PartC/lossy.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.47,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}


=== Dispersive - Dielectric Model
$
  epsilon_r (omega) = epsilon_oo + frac(epsilon_s - epsilon_oo, 1 + j omega tau) \
  C'(omega) = C_oo' + frac(Delta C', 1 + j omega tau), quad Delta C' = C_oo' (frac(epsilon_s - epsilon_oo, epsilon_oo)) \
  tau frac(partial P', partial t) + P'(t) = Delta C' V(t) \
  frac(partial V, partial z) = - L' frac(partial I, partial t) - R' I \
  frac(partial I, partial z) = - C_oo' frac(partial V, partial t) - frac(partial P', partial t) - G' V \
  v_(p, oo) = frac(1, sqrt(L' C_oo')), quad v_(p, 0) = frac(1, sqrt(L' (C_oo' + Delta C')))
$
Which is solved by:
#align(center)[
  #rect(
    ```c
    const double dC_p = Capacitance * (eps_s - eps_inf) / eps_inf;
    const double alpha = (2.0 * tau - dt) / (2.0 * tau + dt);
    const double beta_p = (2.0 * dt * dC_p) / (2.0 * tau + dt);
    const double cB = beta_p / (2.0 * dt);
    const double cP = (1.0 - alpha) / dt;

    cCon = (dt / (2.0 * Inductance));
    double cCon1 = (1.0 - Resistance * cCon) / (1.0 + Resistance * cCon);
    double cCon2 = (dt / (Inductance * grid->dx)) / (1.0 + Resistance * cCon);
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] =
          cCon1 * grid->Hgrid[j] - cCon2 * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    memcpy(grid->Hgrid_prev, grid->Egrid_prev, sizeof(double) * grid->N);

    for (int j = 1; j < grid->N - 1; j++) {
      double dI = (grid->Hgrid[j] - grid->Hgrid[j - 1]) / grid->dx;
      double V_old = grid->Egrid[j];
      double P_old = grid->Hgrid_prev[j];

      double V_new = (V_old * (Capacitance / dt - cB - Conductance / 2.0) - dI +
                      cP * P_old) /
                      (Capacitance / dt + cB + Conductance / 2.0);

      grid->Egrid_prev[j] = alpha * P_old + beta_p * (V_new + V_old) / 2.0;
      grid->Egrid[j] = V_new;
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]

#figure(
  image(
    "PartC/Dielectric_Voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

#let rows = csv("PartC/dielectric.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.27,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}

And as we can see the output from the dielectric dispersion model is much more smeared, since the phase velocities for different frequencies aren't the same anymore, like in the case of Lossless or
Low-Loss case.

=== Comparison between the Theoretical/Analytical and Simulated results

As we can see from the above graphs, there's *severe ringing visible in the output Voltage and Current* which is the result of us *not having a smooth excitation source*, which will not contain the
infinite series of harmonics that can't exactly be represented in the finite discrete grid we have, and hence we can't expect the sharp rise and fall times in the voltage. #underline[A way to reduce
  the effect of this would be to use more finer meshing] but that is prohibitively expensive as we'd need much more compute time, and reality also won't let us have the ideal theoretical response we are
hoping for.

On the plus side, we see that there's a *very nice match between the expected amplitudes of Lossless Simulation and Analytical data*.

#v(1em)
#par(first-line-indent: 0em)[
  #text("GitHub Repository", weight: "bold", size: 16pt)
]
#v(0.2em)
The github repository can be accessed #link("https://github.com/ShivaVarunAdicherla/EE3140-PA", "here").
It will include a readme for this Programming assignment holding all the animated Gifs, and
the requested source code.

#par(first-line-indent: 0em)[
  #text("AI Disclosure", weight: "bold", size: 16pt)
]
#v(0.2em)
There was use of Generative AI (ChatGPT) to solve the Questions - 1 (Understanding the FEM and writing code for it, in python), and Question - 2 (For the smith Chart diagrams and the plots) and
Question - 3 (Idea to solve the Dielectric Model part), Code completely written by hand for this one.

#v(1em)
#line(length: 100%)
#align(center)[
  #text("END", weight: "bold", size: 14pt)
]

