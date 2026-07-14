# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

# Files (mostly Makefiles) to be copied to the build directories
# NB: make.depend files and Makefiles including them are copied in makedeps.sh, not here

AC_DEFUN([X_AC_QE_MAKE_TREE], [

AC_CONFIG_FILES([install/extlibs_makefile:$srcdir/install/extlibs_makefile])
AC_CONFIG_FILES([install/install_utils:$srcdir/install/install_utils])
AC_CONFIG_FILES([install/plugins_makefile:$srcdir/install/plugins_makefile])
AC_CONFIG_FILES([install/plugins_list:$srcdir/install/plugins_list])
AC_CONFIG_FILES([install/tldeps:$srcdir/install/tldeps])

AC_CONFIG_FILES([Makefile:$srcdir/Makefile])

AC_CONFIG_FILES([COUPLE/Makefile:$srcdir/COUPLE/Makefile])
AC_CONFIG_FILES([COUPLE/src/Makefile:$srcdir/COUPLE/src/Makefile])

AC_CONFIG_FILES([FFTXlib/Makefile:$srcdir/FFTXlib/Makefile])
AC_CONFIG_FILES([FFTXlib/tests/Makefile:$srcdir/FFTXlib/tests/Makefile])

AC_CONFIG_FILES([PIOUD/Makefile:$srcdir/PIOUD/Makefile])
AC_CONFIG_FILES([PIOUD/src/Makefile:$srcdir/PIOUD/src/Makefile])
AC_CONFIG_FILES([PIOUD/src/make.depend:$srcdir/PIOUD/src/make.depend])

AC_CONFIG_FILES([NEB/Makefile:$srcdir/NEB/Makefile])

AC_CONFIG_FILES([XSpectra/Makefile:$srcdir/XSpectra/Makefile])

AC_CONFIG_FILES([XClib/Makefile:$srcdir/XClib/Makefile])

AC_CONFIG_FILES([Modules/Makefile:$srcdir/Modules/Makefile])

AC_CONFIG_FILES([atomic/Makefile:$srcdir/atomic/Makefile])

AC_CONFIG_FILES([PP/Makefile:$srcdir/PP/Makefile])
AC_CONFIG_FILES([PP/simple_transport/src/Makefile:$srcdir/PP/simple_transport/src/Makefile])

AC_CONFIG_FILES([QEHeat/Makefile:$srcdir/QEHeat/Makefile])

AC_CONFIG_FILES([GWW/Makefile:$srcdir/GWW/Makefile])
AC_CONFIG_FILES([GWW/util/Makefile:$srcdir/GWW/util/Makefile])
AC_CONFIG_FILES([GWW/minpack/Makefile:$srcdir/GWW/minpack/Makefile])

AC_CONFIG_FILES([HP/Makefile:$srcdir/HP/Makefile])

AC_CONFIG_FILES([PW/Makefile:$srcdir/PW/Makefile])

AC_CONFIG_FILES([PWCOND/Makefile:$srcdir/PWCOND/Makefile])

AC_CONFIG_FILES([KS_Solvers/Makefile:$srcdir/KS_Solvers/Makefile])

AC_CONFIG_FILES([CPV/Makefile:$srcdir/CPV/Makefile])

# TG: TDDFPT/ColorCalculator/Install?
AC_CONFIG_FILES([TDDFPT/Makefile:$srcdir/TDDFPT/Makefile])

# TG: autotest.inc?
AC_CONFIG_FILES([UtilXlib/tests/Makefile:$srcdir/UtilXlib/tests/Makefile])

AC_CONFIG_FILES([LAXlib/tests/Makefile:$srcdir/LAXlib/tests/Makefile])

AC_CONFIG_FILES([PHonon/Makefile:$srcdir/PHonon/Makefile])

AC_CONFIG_FILES([KCW/Makefile:$srcdir/KCW/Makefile])

AC_CONFIG_FILES([EPW/Makefile:$srcdir/EPW/Makefile])
AC_CONFIG_FILES([EPW/irobjs/Makefile:$srcdir/EPW/irobjs/Makefile])

AC_CONFIG_FILES([archive/README.md:$srcdir/archive/README.md])

]
)
