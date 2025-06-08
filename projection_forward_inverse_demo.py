from projections import (
    gnomonic_projection,
    inverse_gnomonic_projection,
    azimuthal_equidistant_projection,
    inverse_azimuthal_equidistant_projection,
    orthographic_projection,
    inverse_orthographic_projection,
)

def print_result(name, lon, lat, x, y, lon2, lat2):
    print(f"==== {name} ====")
    print("原始经纬度:      ", lon, lat)
    print("投影后平面坐标:  ", x, y)
    print("逆投影回经纬度:  ", lon2, lat2)
    print()

def main():
    center_lon, center_lat = 0.0, 20.0   # 投影中心 (A, D)

    # 方式一：输入绝对经纬度
    lon_abs = [0.0, 180.0]
    lat_abs = [10.0, -20.0]

    # 方式二：输入相对坐标（即与(A, D)点的平面坐标，例如用gnomonic_projection等算出来的ξ、η等）
    # 你可以先用一组绝对经纬度反投影到平面坐标（模拟“测量”），再用这些平面坐标测试逆变换
    rel_x = [0.5, -0.5]
    rel_y = [0.1, -0.1]

    print("【方式一：绝对经纬度 → 投影 → 逆投影】")
    for name, proj, invproj in [
        ("Gnomonic", gnomonic_projection, inverse_gnomonic_projection),
        ("Azimuthal Equidistant", azimuthal_equidistant_projection, inverse_azimuthal_equidistant_projection),
        ("Orthographic", orthographic_projection, inverse_orthographic_projection),
    ]:
        x, y = proj(center_lon, center_lat, lon_abs, lat_abs)
        lon2, lat2 = invproj(center_lon, center_lat, x, y)
        print_result(name, lon_abs, lat_abs, x, y, lon2, lat2)

    print("【方式二：相对坐标(平面) → 逆投影 → 投影】")
    for name, proj, invproj in [
        ("Gnomonic", gnomonic_projection, inverse_gnomonic_projection),
        ("Azimuthal Equidistant", azimuthal_equidistant_projection, inverse_azimuthal_equidistant_projection),
        ("Orthographic", orthographic_projection, inverse_orthographic_projection),
    ]:
        # 先用逆变换（平面→球面），再用正变换（球面→平面），应该能回到原始平面坐标
        lon2, lat2 = invproj(center_lon, center_lat, rel_x, rel_y)
        x2, y2 = proj(center_lon, center_lat, lon2, lat2)
        print_result(name, lon2, lat2, x2, y2, lon2, lat2)
        print(f"输入的相对坐标:   {rel_x} {rel_y}")
        print(f"正逆投影后坐标:   {x2} {y2}\n")

if __name__ == "__main__":
    main()

