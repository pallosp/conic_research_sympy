def ConicCenter(conic):
    x, y, z = conic.row(0).cross(conic.row(1))
    return (x / z, y / z)
