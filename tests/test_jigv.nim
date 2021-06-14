import unittest
import pedfile
import ../jigv


suite "jigv suite":
  test "that first affected works":

    var samples = parse_ped("tests/a.ped")
    check samples.first_affected_or_zero == 2

  test "that with no affected sample index is zero":
    var samples = parse_ped("tests/a.ped")
    samples[2].affected = false
    check samples.first_affected_or_zero == 0
