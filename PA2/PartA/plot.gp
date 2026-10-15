# Some of the code here is almost a straight copy of examples in gnuplotting.org
if (!exists("OUTDIR")) OUTDIR = "results"
LONG  = OUTDIR."/longitudinal.dat"
TRANS = OUTDIR."/transverse.dat"
WALL  = OUTDIR."/boundary.dat"

stats WALL using 1 name "LX" nooutput
stats WALL using 2 name "LY" nooutput

# Find extrema of the wall geometry
stats WALL using 3 name "WX" nooutput
stats WALL using 4 name "WY" nooutput

# Combine the extrema from both datasets
XMIN = (LX_min < WX_min ? LX_min : WX_min)
XMAX = (LX_max > WX_max ? LX_max : WX_max)
YMIN = (LY_min < WY_min ? LY_min : WY_min)
YMAX = (LY_max > WY_max ? LY_max : WY_max)

# Set the ranges automatically
set xrange [1.5*XMIN:XMAX*1.5]
set yrange [1.5*YMIN:YMAX*1.5]

FDUMMY = 0
LONG_EXT = sprintf( \
    "< (cat '%s'; printf '%%.16g %%.16g %%.16g\\n' \
    %g %g %g %g %g %g %g %g %g %g %g %g)", \
    LONG, \
    XMIN*5, YMIN*5, FDUMMY, \
    XMIN*5, YMAX*5, FDUMMY, \
    XMAX*5, YMIN*5, FDUMMY, \
    XMAX*5, YMAX*5, FDUMMY)

set palette defined (0 "#00008f", 0.25 "#0055ff", 0.5 "#00ffff", 0.75 "#ffff00", 1 "#ff0000")

set terminal pdfcairo size 6in,3in enhanced font "Sans,11"
set size ratio -1
set output OUTDIR."/heatmap.pdf"
set title "Re longitudinal field at z = 0"
set xlabel "x (m)"; set ylabel "y (m)"
set cblabel "Re{psi}"
set pm3d interpolate 0,0
set dgrid3d 200,200 qnorm 2
set view map
unset key
splot LONG_EXT using 1:2:3 with pm3d, \
      WALL using 1:2:($3-$1):($4-$2) with vectors nohead lw 2 lc "black"
unset dgrid3d
unset pm3d

set xrange [1.2*XMIN:XMAX*1.2]
set yrange [1.2*YMIN:YMAX*1.2]

set key inside right top
set offset graph 0.01, graph 0.01, graph 0.01, graph 0.01
set output OUTDIR."/transverse_E.pdf"
set title "Transverse E_t (colour = |E_t|)"
set xlabel "x (m)"; set ylabel "y (m)"
set cblabel "|E_t| (arrow units)"
plot TRANS using ($1-$3/2):($2-$4/2):($3):($4):($7) \
       with vectors head size screen 0.012,7,60 filled lc palette title "E_t"

set output OUTDIR."/transverse_H.pdf"
set title "Transverse H_t (colour = |H_t|)"
set cblabel "|H_t| (arrow units)"
plot TRANS using ($1-$5/2):($2-$6/2):($5):($6):($8) \
       with vectors head size screen 0.012,7,60 filled lc palette title "H_t", \
     WALL using 1:2:($3-$1):($4-$2) with vectors nohead lw 1.5 lc "black" title "PEC"

set output OUTDIR."/transverse.pdf"
set title "Transverse E_t (red) and H_t (blue)"
unset cblabel
plot TRANS using ($1-$3/2):($2-$4/2):($3):($4) \
       with vectors head size screen 0.012,7,60 filled lc "red" title "E_t", \
     TRANS using ($1-$5/2):($2-$6/2):($5):($6) \
       with vectors head size screen 0.012,7,60 filled lc "blue" title "H_t", \
     WALL using 1:2:($3-$1):($4-$2) with vectors nohead lw 1.5 lc "black" title "PEC"
