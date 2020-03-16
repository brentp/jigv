version       = "0.0.1"
author        = "Brent Pedersen"
description   = "igv.js server"
license       = "MIT"


# Dependencies

requires "hts >= 0.3.4", "jester"
requires "argparse >= 0.7.0"
installExt = @["nim"]

skipDirs = @["tests"]

