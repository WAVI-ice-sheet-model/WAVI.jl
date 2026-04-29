
xx= simulation.model.grid.xxh/1000
yy = simulation.model.grid.yyh/1000
s = simulation.model.fields.gh.s
u = simulation.model.fields.gh.u
v = simulation.model.fields.gh.v
Φ = simulation.model.fields.g3d.Φ
arrow_length = 100/1000
n_arrows = 50
ss = randperm(length(xx))[1:n_arrows]
clim = (0,1.0)
cmap = colormap("Blues", logscale = true)
p1=heatmap(xx[:,1],yy[1,:],Φ[:,:,1]',clim = clim,xlabel="x (km)",ylabel="y (km)", colormap = cmap, framestyle = :box)
p2=heatmap(xx[:,1],yy[1,:],Φ[:,:,2]',clim = clim,xlabel="x (km)",ylabel="y (km)", colormap = cmap, framestyle = :box)
p3=heatmap(xx[:,1],yy[1,:],Φ[:,:,3]',clim = clim,xlabel="x (km)",ylabel="y (km)", colormap = cmap, framestyle = :box)
p4=heatmap(xx[:,1],yy[1,:],Φ[:,:,4]',clim = clim,xlabel="x (km)",ylabel="y (km)", colormap = cmap, framestyle = :box)
p5=contour(xx[:,1],yy[1,:],s', levels = 0:100:2000,xlabel="x (km)",ylabel="y (km)", framestyle = :box)
p5 = quiver!(p5,xx[ss],yy[ss],quiver=(arrow_length*u[ss], arrow_length*v[ss]), color = :red, arrow = Plots.arrow(:closed, :head, 0.01, 0.01))
plot(p5,p4,p3,p2,p1;layout=(5,1),size=(700,1000),left_margin=10Plots.mm,right_margin=10Plots.mm)
