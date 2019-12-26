#https://rpubs.com/tjmahr/sem_diagrammer  <-- cool example of how to extract a diagram directly from a Lavaan object.


library(DiagrammeR)
grViz("
digraph SEM {

graph [layout = neato,
       overlap = true,
       outputorder = edgesfirst]

node [shape = rectangle]

a [pos = '-4,1!', label = 'e1', shape = circle]
b [pos = '-3,1!', label = 'ind_1']
c [pos = '-3,0!', label = 'ind_2']
d [pos = '-3,-1!', label = 'ind_3']
e [pos = '-1,0!', label = 'latent a', shape = ellipse]
f [pos = '1,0!', label = 'latent b', shape = ellipse]
g [pos = '1,1!', label = 'e6', shape = circle]
h [pos = '3,1!', label = 'ind_4']
i [pos = '3,-1!', label = 'ind_5']
j [pos = '4,1!', label = 'e4', shape = circle]
k [pos = '4,-1!', label = 'e5', shape = circle]

a->b
e->b [label = '0.6']
e->c [label = '0.6']
e->d [label = '0.6']

e->f [label = '0.321', headport = 'w']
g->f [tailport = 's', headport = 'n']

d->c [dir = both]

f->h [label = '0.6', tailport = 'ne', headport = 'w']
f->i [label = '0.6']

j->h
k->i

}
")

