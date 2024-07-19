using DynamicalSystems
using GLMakie

function three_body!(du, u, p, t)

	# position, velocity
	x1, y1, x2, y2, x3, y3,
	vx1, vy1, vx2, vy2, vx3, vy3 = u

	# mass
	m1, m2, m3 = p

	# velocity
	du[1] = vx1
	du[2] = vy1

	du[3] = vx2
	du[4] = vy2

	du[5] = vx3
	du[6] = vy3

	# distance from body to body
	r12 = sqrt((x1 - x2)^2 + (y1 - y2)^2)
	r13 = sqrt((x1 - x3)^2 + (y1 - y3)^2)
	r23 = sqrt((x2 - x3)^2 + (y2 - y3)^2)

	# acceleration (set G = 1)
	du[7]  = -m2 * ((x1 - x2) / (r12^3)) - m3 * ((x1 - x3) / (r13^3))
	du[8]  = -m2 * ((y1 - y2) / (r12^3)) - m3 * ((y1 - y3) / (r13^3))

	du[9]  = -m3 * ((x2 - x3) / (r23^3)) - m1 * ((x2 - x1) / (r12^3))
	du[10] = -m3 * ((y2 - y3) / (r23^3)) - m1 * ((y2 - y1) / (r12^3))

	du[11] = -m1 * ((x3 - x1) / (r13^3)) - m2 * ((x3 - x2) / (r23^3))
	du[12] = -m1 * ((y3 - y1) / (r13^3)) - m2 * ((y3 - y2) / (r23^3))
end

function dotsize(m)
	return maximum((10 * (log10(m) + 1), 10.0))
end


no_animar = false
masses = [100.0, 0.1, 0.1]

v1 = 11.0
v2 = 9.0
initial_velocities = [0.0, 0.0,
                      0.0, v1,
                      0.0, -v2]

initial_positions = [0.05, 0.0, masses[1]/(v1^2), 0.0, -masses[1]/(v2^2), 0.0]

# initial_velocities = randn(6)/10
# initial_positions = [-1, 0, 1, 0, 0, sqrt(5)]

initial_momentum = [sum(masses .* initial_velocities[1:2:end]),
                    sum(masses .* initial_velocities[2:2:end])]

initial_velocities[1:2:end] = initial_velocities[1:2:end] .- initial_momentum[1] ./ sum(masses)
initial_velocities[2:2:end] = initial_velocities[2:2:end] .- initial_momentum[2] ./ sum(masses)

initial_conditions = [initial_positions; initial_velocities]
ds = ContinuousDynamicalSystem(three_body!, initial_conditions, masses)

tspan = (0.0, 1000.0)

tr, t = trajectory(ds, 10.0;Δt = 0.001)

if no_animar
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, tr[:, 1], tr[:, 2], label = "Body 1")
    lines!(ax, tr[:, 3], tr[:, 4], label = "Body 2")
    lines!(ax, tr[:, 5], tr[:, 6], label = "Body 3")
    fig
else

	point1 = Observable(Point2f[])
	point2 = Observable(Point2f[])
	point3 = Observable(Point2f[])
	planets = Observable([Point2f(initial_positions[1], initial_positions[2]),
						  Point2f(initial_positions[3], initial_positions[4]),
						  Point2f(initial_positions[5], initial_positions[6])])

	colors = Observable(Int[])

	set_theme!(theme_black())

	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1)
	scatter!(planets, markersize = dotsize.(masses), color = [:yellow, :brown, :brown])
	l = lines!(ax, point1, color = colors, colormap = :inferno, transparency = true)
	lines!(ax, point2, color = colors, colormap = :inferno, transparency = true)
	lines!(ax, point3, color = colors, colormap = :inferno, transparency = true)

	ax.xlabel = L"$x$"
	ax.ylabel = L"$y$"
				

	x1, y1 = tr[:, 1], tr[:, 2]
	x2, y2 = tr[:, 3], tr[:, 4]
	x3, y3 = tr[:, 5], tr[:, 6]

	final_frame = 10000
	trail = 250
	skip_every = 4
	record(fig, "3-body.mp4", 1:skip_every:final_frame; framerate=60) do frame
		planets[] = [Point2f(x1[frame], y1[frame]),
		Point2f(x2[frame], y2[frame]),
		Point2f(x3[frame], y3[frame])]
		for i in 1:skip_every
			push!(point1[], Point2f(x1[frame+i], y1[frame+i]))
			push!(point2[], Point2f(x2[frame+i], y2[frame+i]))
			push!(point3[], Point2f(x3[frame+i], y3[frame+i]))
			push!(colors[], frame)
		end
		if frame > trail*skip_every
			for i in 1:skip_every
				popfirst!(point1[])
				popfirst!(point2[])
				popfirst!(point3[])
				popfirst!(colors[])
			end
		end
		notify(point1); notify(point2); notify(point3); notify(colors) # tell points and colors that their value has been updated
		#Rescale axis:
		trail_index = maximum((frame-trail*skip_every, 1))
		l.colorrange = (trail_index, frame) # update plot attribute directly

		autolimits!(ax)

	end
end