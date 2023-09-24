
ball=sphere(pos=vector(0,1,0), radius=0.50, color=color.white, make_trail=True)
ground= box(pos=vector(0,0,0), size=vector(200,0.1,0.80), color=color.green)
g1=graph(xtitle="t [s]", ytitle="y [m]")
f1 = gcurve(color=color.red)
g2=graph(xtitle="t [s]", ytitle="y [m]")
f2 = gcurve(color=color.blue)
g=vector(0,-9.8,0)
ball.m=60
v0=274
tetha=pi/2
ball.p=ball.m*v0*vector(cos(tetha),sin(tetha),0)
vscale=.001
varrow=arrow(pos=ball.pos,axis=vscale*ball.p, color=color.blue)

ball2=sphere(pos=vector(50,1,0), radius=0.50, color=color.blue, make_trail=True)
ball2.m=60
ball2.p=ball2.m*v0*vector(cos(tetha),sin(tetha),0)
A=pi*ball.radius**2
pho=.8
c=0.30

t=0
dt=0.01

while ball.pos.y>0.9:
  rate(100)
  F=ball.m*g
  F2=ball2.m*g -.5*c*pho*A*mag(ball2.p)**2*norm(ball2.p)/ball2.m**2
  ball.p= ball.p + F*dt
  ball2.p=ball2.p + F2*dt
  ball.pos= ball.pos + ball.p*dt/ball.m
  ball2.pos= ball2.pos + ball2.p*dt/ball2.m
  #varrow.pos=ball.pos
  #varrow.axis=vscale*ball.v
  t= t +dt
  
  f1.plot(t,ball.pos.y)
  f2.plot(t,ball2.pos.y)