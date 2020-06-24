AC_DEFUN([PUFF_FIND_NETCDF],
[# make precious variables for ./configure --help
AC_ARG_VAR(NETCDF_LIB,location of netCDF library)
AC_ARG_VAR(NETCDF_INC,location of netCDF header files)
# add user-specified netCDF locations if specified
# otherwise add /usr/local/ locations if they exist
if test "$NETCDF_LIB"; then
  if test -d "$NETCDF_LIB"; then
    LDFLAGS="$LDFLAGS -L$NETCDF_LIB"
  else
    echo "*** NETCDF_LIB location \"$NETCDF_LIB\" does not exist. ***"
  fi
else
  if test -d "/usr/local/lib"; then
    LDFLAGS="$LDFLAGS -L/usr/local/lib "
  fi
fi
if test "$NETCDF_INC"; then
  if test -d "$NETCDF_INC"; then
    CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"
  else
    echo "*** NETCDF_INC location \"$NETCDF_INC\" does not exist. ***"
  fi
dnl else
dnl   if test -d "/usr/local/include"; then
 dnl   CPPFLAGS="$CPPFLAGS -I/usr/local/include "
dnl  fi
fi
])

