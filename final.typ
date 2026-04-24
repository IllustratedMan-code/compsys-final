#import "@local/assignments:1.0.0": conf
#import "@preview/subpar:0.2.2"
#show: doc => conf(title: "Computational Systems Biology Final Assignment", doc)



@f1 (the replication of Figure 2 in Dovzhenok et al.) was the easiest to create, notably, as long as the initial conditions were above a certain threshold, they eventually settled into this occilatory behavior. I used the `DifferentialEquations` package in Julia for @f1 and it was very straightforward. Plotting is done using the `Plots` package.

#subpar.grid(
    columns: 2,
    figure(
        image("Figure-2A.svg"),
        caption: []
    ), <f1a>,
    figure(
        image("Figure-2B.svg"),
        caption: []
    ), <f1b>,
    figure(
        image("Figure-2C.svg"),
        caption: []
    ), <f1c>,
    label: <f1>,
)

#pagebreak()

@f2 was not so easy to create. I used the `BifurcationKit` package for this and followed the general procedure for these plots. First, I found the Hopf points, then used `BifurcationKit` to simulate the model using parameter values along a given range. Since there is no Hopf point for @f2c, I instead sampled a simulation (as in @f1) as the "initial guess" for `BifurcationKit`.  Finding the Hopf points was not very intuitive as whatever algorithm BifurcationKit uses seems to be very sensitive to initial conditions. So if $k_03$ for example was not in a particular range, then the Hopf point would not be found, even though it exists.  I found out that if I boosted the initial conditions, by setting the initial concentrations to large values, the algorithm became less sensitive to these parameter changes, something which was necessary to reproduce @f3.

#subpar.grid(
    columns: 2,
    figure(
        image("Figure-3A.svg"),
        caption: []
    ), <f2a>,
    [#figure(
        image("Figure-3B.svg"),
        caption: []
    ) <f2b>],
    [#figure(
        image("Figure-3C.svg"),
        caption: []
        ) <f2c>],
    figure(
        image("Figure-3D.svg"),
        caption: []
    ),
    label: <f2>,
    caption: [Reproduction of Figure 3]
)
#grid( columns: 2,
    [
        @f3 uses essentially the same process to reproduce the figure, except the black lines, as in @f3c, needed to be extracted from the function that found the Hopf points. Luckily, these values were fairly easy to find in the generated object.

        To note, I needed to set all the initial concentrations to 10 to get the full range of the black lines in @f3b and @f3c for some reason, otherwise the Hopf points would fail to be calculated. 
    ],

[#subpar.grid(
    columns: 1,
    figure(
        image("Figure-4A.svg", width:100%),
        caption: []
    ), <f3a>,
    [#figure(
        image("Figure-4B.svg", width:100%),
        caption: []
    ) <f3b>],
    [#figure(
        image("Figure-4C.svg", width:100%),
        caption: []
    ) <f3c>],
    label: <f3>,
    caption: [Reproduction of Figure 4]
)
]
)

= Code

Nicer code link: https://illustratedman-code.github.io/compsys-final/

#let code = read("notebook.jl")
#let lines = code.split("\n").slice(0, 457)
#let code = lines.join("\n")

#raw(lang: "julia", block: true, code)

