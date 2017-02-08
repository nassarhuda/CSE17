import PyPlot
using Compose
using Colors


A0 =[
#    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0#1
     0 0 1 0 1 0 0 0 0 0 0 1 0 0 0#2
     1 1 0 1 0 0 0 0 0 0 0 0 0 0 0#3
     0 0 1 0 1 1 0 0 0 0 0 0 0 0 0#4
     0 1 0 1 0 1 0 0 0 0 0 0 0 0 0#5
     0 0 0 1 1 0 1 0 1 0 1 0 0 0 0#6
     0 0 0 0 0 1 0 1 0 0 0 0 0 0 0#7
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0#8
     0 0 0 0 0 1 0 0 0 1 0 1 0 0 0#9
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0#10
     0 0 0 0 0 1 0 0 0 0 0 0 1 1 1#11
     0 1 0 0 0 0 0 0 1 0 0 0 1 0 0#12
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0#13
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 1#14
     0 0 0 0 0 0 0 0 0 0 1 0 0 1 0#15
     ]
 
xyA0 =[
    1.4000    2.5000#1
    1.3000    2.1000#2
    1.4000    2.3000#3
    1.5000    2.1000#4
    1.4000    1.9000#5
    1.5000    1.7000#6
    1.3000    1.6000#7
    1.2000    1.2000#8
    1.6000    1.9000#9
    1.6500    2.1000#10
    1.6000    1.5000#11
    1.7000    1.7500#12
    1.6600    1.6000#13
    1.5000    1.3000#14
    1.7000    1.0000#15
    ]

n = size(A0,1)
assert(isequal(A0,A0'))
B0 =[
#    1 2 3 4 5 6 7 9 0 1 2 3 4 5
     0 0 1 0 0 0 0 0 0 0 0 0 0 0#1
     0 0 1 0 1 0 0 0 0 0 0 0 0 0#2
     1 1 0 1 0 0 0 0 0 0 0 0 0 0#3
     0 0 1 0 1 1 0 0 1 0 0 0 0 0#4
     0 1 0 1 0 1 1 0 0 0 0 0 0 0#5
     0 0 0 1 1 0 0 1 0 1 0 0 0 0#6
     0 0 0 0 1 0 0 0 0 0 0 0 0 0#7
     0 0 0 0 0 1 0 0 1 0 1 0 0 0#9
     0 0 0 1 0 0 0 1 0 0 0 0 1 0#10
     0 0 0 0 0 1 0 0 0 0 0 1 1 1#11
     0 0 0 0 0 0 0 1 0 0 0 1 0 0#12
     0 0 0 0 0 0 0 0 0 1 1 0 0 0#13
     0 0 0 0 0 0 0 0 1 1 0 0 0 1#14
     0 0 0 0 0 0 0 0 0 1 0 0 1 0#15
     ];
 
xyB0 =[
    1.4000    2.5000#1
    1.3000    2.1000#2
    1.4000    2.3000#3
    1.5000    2.1000#4
    1.4000    1.9000#5
    1.5000    1.7000#6
    1.3000    1.6000#7
    1.5700    1.9000#9
    1.6500    2.1000#10
    1.6000    1.5000#11
    1.7000    1.7500#12
    1.6600    1.6000#13
    1.5000    1.3000#14
    1.7000    1.0000#15
    ];

n = size(B0,1)
assert(isequal(B0,B0'))

A = A0
B = B0
xyA = xyA0
xyB = xyB0

xyA[:,1] = xyA0[:,1] - 1
xyA[:,2] = xyA0[:,2] - .5

xmin,xmax = extrema([xyA[:,1];xyB[:,1]])
ymin,ymax = extrema([xyA[:,2];xyB[:,2]])


keepA = [1,2,3,4,5,6,9,11,12,14,15]
keepB = [1,2,3,4,5,6,9,10,11,13,14]

extraA = [9,12,14]
extraB = [5,12,14]
X = zeros(Float64,length(keepA)+length(extraA),2)
Y = zeros(Float64,length(keepA)+length(extraA),2)

X[:,1] = [xyA[keepA,1];xyA[extraA,1]]
X[:,2] = [xyB[keepB,1];xyB[extraB,1]]

Y[:,1] = [xyA[keepA,2];xyA[extraA,2]]
Y[:,2] = [xyB[keepB,2];xyB[extraB,2]]

points = [[(X[i,1], Y[i,1]), (X[i,2], Y[i,2])] for i in 1:size(X,1)]
all_edges = line(points)
lw = 0.5
opacity = 0.5
graph = compose(context(),all_edges, linewidth(lw), fill(nothing), stroke(colorant"black"), strokeopacity(opacity),
     strokedash([1mm, 1mm]),)

area = UnitBox(xmin,ymax,xmax-xmin,ymin-ymax,leftpad=2mm,rightpad=2mm,toppad=2mm,bottompad=2mm)

opacity = 0.99

n = size(A0,1)
ei,ej = findnz(A)
lx = [xyA[ei,1]';xyA[ej,1]']
ly = [xyA[ei,2]';xyA[ej,2]']
len = size(lx,2)
points = [[(lx[1,i], ly[1,i]), (lx[2,i], ly[2,i])] for i in 1:len];
all_edges = line(points)
radius = 4.5pt
graph2 = compose(context(),all_edges, linewidth(lw), fill(nothing), stroke(colorant"green"),
    strokeopacity(opacity))
all_nodes = circle(xyA[:,1],xyA[:,2],[radius])
plot2 = compose(context(), all_nodes, fill(colorant"green"), stroke(nothing))

n = size(B0,1)
ei,ej = findnz(B)
lx = [xyB[ei,1]';xyB[ej,1]']
ly = [xyB[ei,2]';xyB[ej,2]']
len = size(lx,2)
points = [[(lx[1,i], ly[1,i]), (lx[2,i], ly[2,i])] for i in 1:len];
all_edges = line(points)
radius = 4.5pt
graph3 = compose(context(),all_edges, linewidth(lw), fill(nothing), stroke(colorant"red"),
    strokeopacity(opacity))
all_nodes = circle(xyB[:,1],xyB[:,2],[radius])
plot3 = compose(context(), all_nodes, fill(colorant"red"), stroke(nothing))

fig = compose(context(units=area), graph, graph2, plot2, graph3, plot3)

img = PDF("figures/netalign.pdf",150mm,150mm)
draw(img,fig)