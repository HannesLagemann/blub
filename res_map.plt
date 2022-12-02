load './style.plt'
cd working_directory = system("dirname ".ARG0)."/"
lw_size=0
l_type=1
pt_size=1.5
TicsFontSize="Times-Roman,18pt"
LabelFontSize="Times-Roman,20pt"
CornerLabelFontSize="Times-Roman,50pt"
set key font "".TicsFontSize

set style line 1 lt -1 lw 1.5 lc rgbcolor "dark-khaki" pt 1 ps 1
set style line 2 lt -1 lw 1.5 lc rgbcolor "dark-green" pt 2 ps 1
set style line 3 lt -1 lw 1.5 lc rgbcolor "blue" pt 3 ps 1
set style line 4 lt -1 lw 1.5 lc rgbcolor "coral" pt 4 ps 1
set style line 5 lt -1 lw 1.5 lc rgbcolor "magenta" pt 5 ps 1
set style line 6 lt -1 lw 1.5 lc rgbcolor "cyan" pt 6 ps 1
set style line 7 lt -1 lw 1.5 lc rgbcolor "purple" pt 7 ps 1
set style line 8 lt -1 lw 1.5 lc rgbcolor "orange" pt 8 ps 1
set style line 9 lt -1 lw 1.5 lc rgbcolor "light-blue" pt 9 ps 1
set style line 10 lt -1 lw 1.5 lc rgbcolor "red" pt 10 ps 1
set style line 11 lt -1 lw 1.5 lc rgbcolor "black" pt 11 ps 1

#set contour base
set view map; set size square
set sample 11; set isosamples 11
set pm3d map at bstbst
set palette
set colorbox size 0.01, 0.01
set lmargin 0
set pm3d flush begin
set view map
set size ratio 1.0
set bmargin 2.0
set lmargin 2.0
set rmargin 2.0
set tmargin 2.0

n=1
#set label "(c): {/Symbol b}=0" font "".CornerLabelFontSize at screen 0.0, 0.95
set xlabel "({/Symbol w}^{(D)}-{/Symbol w}^{(0)})/2 {/Symbol p} [MHz]"  font "".LabelFontSize
set ylabel "t [ns]" font "".LabelFontSize
set cblabel "empty"
set cbrange[0:1]
set grid xtics lc "white"
set grid ytics lc "white"

name=system("ls res_map.dat")
name=name[1:strlen(name)-4]
file=name.".dat"
fac=1000
offset=45.0

set title sprintf("{/Symbol w^{(0)}}/2 {/Symbol p}= %.3F [GHz]",offset)
set cblabel "p^{(0,0,0)}(t,{/Symbol w}^{(D)})" font "".LabelFontSize
suf=name.".eps"
set output suf
splot file every n u ($1-offset)*fac:($2):($3) title ""
unset label
q
