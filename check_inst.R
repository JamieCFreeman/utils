#
# 2020-4-10 JCF
#

#######################
# Helper functions
#######################

#
# Check package install, return TRUE if exists, FALSE if does not
#

check_inst <- function(pkg) {
 
   nzchar( system.file(package = pkg) )

  }


# Source:
# https://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages
