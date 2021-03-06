#https://rpubs.com/tjmahr/sem_diagrammer  <-- cool example of how to extract a diagram directly from a Lavaan object.


library(DiagrammeR)


grViz("
digraph {

graph [layout = neato,
       overlap = true,
       outputorder = edgesfirst]

node [shape = circle ]

crtFac     [pos = '0,3!', label = ' CRT ', shape = circle, fontsize=14 width=1.,fixedsize=true]
allPolFac  [pos = '0,0!', label = ' All \n politics ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwnatFac  [pos = '-0.5,-2.5!', label = ' M&W \n NatDiff ', shape = circle, fontsize=14 width=1.,fixedsize=true]
fmFac  [pos = '-4,1.2!', label = ' Free \n market ', shape = circle, fontsize=14 width=1.,fixedsize=true]
evoFac  [pos = '4,-1.2!', label = ' Evolution ', shape = circle, fontsize=14 width=1.,fixedsize=true]
camFac  [pos = '4,0!', label = ' Reject \n CAM ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwevoFac  [pos = '1,-3!', label = ' M&W \n Evo ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwequFac  [pos = '-2,-2!', label = ' M&W \n Equal ', shape = circle, fontsize=14 width=1.,fixedsize=true]
religFac  [pos = '-4,-1.2!', label = ' Religiosity ', shape = circle, fontsize=14 width=1.,fixedsize=true]
vaxFac  [pos = '4,1.2!', label = ' Vax ', shape = circle, fontsize=14 width=1.,fixedsize=true]
consFac  [pos = '-4,0!', label = ' Conser- \n vatism ', shape = circle, fontsize=14 width=1.,fixedsize=true]

allPolFac  ->  fmFac  [label = ' 0.66 ' fontsize=18]
allPolFac  ->  consFac  [label = ' 0.70 ' fontsize=18]
allPolFac  ->  religFac  [label = ' 0.59 ' fontsize=18]
fmFac  ->  vaxFac  [label = ' -0.08 ' fontsize=18]
mwequFac  ->  vaxFac  [label = '  0.15 ' fontsize=18]
crtFac  ->  vaxFac  [label = '  0.21 '  fontsize=18]
religFac  ->  evoFac  [label = ' -0.31 ' fontsize=18]
mwevoFac  ->  evoFac  [label = '  0.86 ' fontsize=18]
mwnatFac  ->  evoFac  [label = ' -0.78 ' fontsize=18]
crtFac  ->  evoFac  [label = '  0.25 ' fontsize=18]
allPolFac  ->  camFac  [label = ' -0.20 ' fontsize=18]
crtFac  ->  camFac  [label = '  0.17 ' fontsize=18]
}
")
