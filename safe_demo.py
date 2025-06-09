import projections_safe as ps

A, D = 0.0, 20.0
lon, lat = [0, 180], [-70, -20]

# 包装后会返回 (x, y, status)
x, y, st = ps.gnomonic_projection(A, D, lon, lat)
print("status =", st, "x,y =", x, y)

# 调用逆向
lon2, lat2, st2 = ps.inverse_gnomonic_projection(A, D, x, y)

