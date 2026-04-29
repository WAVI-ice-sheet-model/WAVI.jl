"""
advect3D!(field::AbstractArray{T,3},model::AbstractModel{T,N}) where {T,N}

Advect a 3D field using semi-Lagrangian advection and plug-flow assumption
"""
function advect3D!(field::AbstractArray{T,3},model::AbstractModel{T,N}) where {T,N}
@unpack gh = model.fields
@unpack nx,ny,x0,y0,dx,dy,σ,nσ = model.grid
@unpack dt = model.params

x = range(start=(x0+0.5dx),step=dx,length=nx)
y = range(start=(y0+0.5dy),step=dy,length=ny)
grid_coords = (x,y,σ)
itp = interpolate(grid_coords,field,Gridded(Linear()))

for k=1:nσ
    for j = 1:ny
        for i = 1:nx
            if gh.mask[i,j]
               σ_dot_material = -(σ[k]*gh.accumulation[i,j]+(one(σ[k])-σ[k])*gh.basal_melt[i,j])./gh.h[i,j]
               x_source = x[i] - dt*gh.u[i,j]
               y_source = y[j] - dt*gh.v[i,j]
               σ_source = σ[k] - dt* σ_dot_material
               x_source = clamp(x_source,minimum(x),maximum(x))
               y_source = clamp(y_source,minimum(y),maximum(y))
               σ_source = clamp(σ_source,minimum(σ),maximum(σ))
               field[i,j,k]=itp(x_source,y_source,σ_source)
            end
        end
    end
end

return model

end