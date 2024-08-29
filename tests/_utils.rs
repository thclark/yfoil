
pub mod utils {
    #[macro_export]
    macro_rules! assert_error {
        ($result_or_error:expr, $expected_error:expr) => {
            // Check it's actually an error first
            if !$result_or_error.is_err() {
                panic!("assert_err failed: expression is not an error")
            }

            // Having confirmed it's an error, check it produces the required string
            let actual = $result_or_error.unwrap_err().to_string();
            let expected = $expected_error.to_string();
            if expected != actual {
                panic!("assert_err failed: error message does not match teh expected error\nexpected:   `{}`\nactual:     `{}`", expected, actual)
            }
        }
    }
}