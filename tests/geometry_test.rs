use yfoil::geometry::read_geometry_from_file;
use yfoil::geometry::InvalidGeometryError;
mod _utils;

#[test]
fn panic_on_erroneous_point() {
    let result_or_error =
        read_geometry_from_file("tests/fixtures/aerofoil_with_invalid_last_point.json");
    assert_error!(
        result_or_error,
        InvalidGeometryError::BeginAndEndAtTrailingEdgeError()
    );
}
