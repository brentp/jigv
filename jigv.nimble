version       = "0.1.2"
author        = "Brent Pedersen"
description   = "igv.js static files"
license       = "MIT"


# Dependencies

requires "hts >= 0.3.16", "https://github.com/brentp/pedfile >= 0.0.3"
requires "argparse == 0.10.1", "zippy"
installExt = @["nim"]

skipDirs = @["tests"]
