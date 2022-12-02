#set format x "10^{%L}"
#set format y "10^{%L}"
set term postscript eps enhanced dashed color 24 linewidth 3  font "Times-Roman,24pt"
set grid
set key font "Times-Roman,27pt"
set termoption dashed

lw_size=0
l_type=1
pt_size=2.5

TicsFontSize="Times-Roman,30pt"
LabelFontSize="Times-Roman,35pt"

array Colors[20]
Colors[1]="blue"
Colors[2]="forest-green"
Colors[3]="dark-red"
Colors[4]="dark-violet"
Colors[5]="gray30"
Colors[6]="dark-red"
Colors[7]="blue"
Colors[8]="forest-green"
Colors[9]="brown"
Colors[10]="dark-violet"
Colors[11]="dark-pink"
Colors[12]="black"
Colors[13]="dark-cyan"
Colors[14]="dark-orange"
Colors[15]="gray30"
Colors[16]="brown"
Colors[17]="dark-violet"
Colors[18]="dark-pink"
Colors[19]="blue"
Colors[20]="dark-cyan"






array LptType[5]
LptType[1]=7
LptType[2]=5
LptType[3]=11
LptType[4]=13
LptType[5]=15



array LtType[5]
LtType[1]=1
LtType[2]=4
LtType[3]=2
LtType[4]=1
LtType[5]=6


##{/Symbol d}
