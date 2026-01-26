"""=================================================================================================


Michael Gribskov     11 June 2025
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
test.describe("Example Test Cases")
emulator = Emulator(1, 10, 5)
test.assert_equals(solve(emulator), 5)
emulator = Emulator(1, 10, 8)
test.assert_equals(solve(emulator), 8)
emulator = Emulator(2, 10, 39)
test.assert_equals(solve(emulator), 39)
emulator = Emulator(2, 16, 123)
test.assert_equals(solve(emulator), 123)

