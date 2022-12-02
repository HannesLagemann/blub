load './style.plt'
cd working_directory = system("dirname ".ARG0)."/"

MarkerFontSize="Times-Roman,20pt"
TicsFontSize="Times-Roman,30pt"
LabelFontSize="Times-Roman,30pt"
CornerLabelFontSize="Times-Roman,50pt"

set xlabel "{/Symbol j}/2 {/Symbol p}"
set ylabel "(E_{~z{.4-}}({/Symbol j})-E_{0}({/Symbol j}))/2 {/Symbol p} [GHz]"
lw_size=1.5
pt_size=1.2
set key samplen 0.01
set key font ",11.5"
set key above
set key vertical maxrows 2
set key horizontal maxcols 7
set xtics

titles = "'~z{.4-} = (0,0,0)'      '~z{.4-} = (0,0,1)'      '~z{.4-} = (0,1,0)'      '~z{.4-} = (0,0,2)'      '~z{.4-} = (0,1,1)'      '~z{.4-} = (0,2,0)'      '~z{.4-} = (0,0,3)'      '~z{.4-} = (0,1,2)'      '~z{.4-} = (0,2,1)'      '~z{.4-} = (0,3,0)' "
set tmargin 3.0
set key at screen 0.54, 0.95 center vertical height 1.0 box lw 1.75 maxrows 2 width 1.0 #font "".TicsFontSize
name="spec"
file=name.".dat"
suf=name.".eps"
set output suf

n=1
y_hight=0.885
n=50

unset grid
plot for [i=3:16] file every 1 using ($1):(column(i)) title "" axes x1y1 with l lt 1 lw lw_size lc rgb Colors[i-2],\
     for [i=3:16] file every n using ($1):(column(i)) title word(titles, i) axes x1y1 with p pt LptType[(i-2)%5+1] pointsize pt_size lc rgb Colors[i-2]
set output
q
