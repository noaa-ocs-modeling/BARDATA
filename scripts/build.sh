#!/bin/bash

###########################################################################
### Author:  Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
###
### Version - 1.0 Mon Nov 16 2020
###########################################################################


###====================
# Make sure that the current working directory is in the PATH
[[ ! :$PATH: == *:".":* ]] && export PATH="${PATH}:."


# Get the directory where the script is located
if [[ $(uname -s) == Darwin ]]; then
#  readonly scrDIR="$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)"
  readonly scrNAME="$( grealpath -s "${BASH_SOURCE[0]}" )"
  readonly scrDIR="$(cd "$(dirname "${scrNAME}" )" && pwd -P)"
else
#  readonly scrDIR="$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)"
  readonly scrNAME="$( realpath -s "${BASH_SOURCE[0]}" )"
  readonly scrDIR="$(cd "$(dirname "$(realpath -s "${BASH_SOURCE[0]}")" )" && pwd -P)"
fi

lst="${scrDIR}/functions_build ${scrDIR}/scripts/functions_build functions_build "
funcs=
for ilst in ${lst}
do
  if [ -f "${ilst:-}" ]; then
    funcs="${ilst}"
    break
  fi
done

if [ -n "${funcs:+1}" ]; then
  source "${funcs}"
else
  echo " ### ERROR :: in ${scrNAME}"
  echo "     Cannot load the required file: functions_build"
  echo "     Exiting now ..."
  echo
  exit 1
fi

unset ilst funcs
###====================


############################################################
### BEG:: SYSTEM CONFIGURATION
############################################################

# Call ParseArgs to get the user input.
ParseArgs "${@}"

# Set the variables for this script
getEnvVars ${scrDIR}

# Check if the user supplied valid components
checkComponents

# Get the compilers to use for this project compilation
getCompilerNames "${COMPILER}"

############################################################
### END:: SYSTEM CONFIGURATION
############################################################


############################################################
### START THE CALCULATIONS
############################################################

if [ ${CLEAN:-0} -ne 0 ]; then
  if [ ${CLEAN:-0} -gt 0 ]; then
    echo "User requested to only clean the project files and then exit."
  else
    echo "User requested to clean the project files."
  fi

  for icomp in ${COMPONENT}
  do
   [ "${icomp}" == "OGCMDL" ] &&  compileOGCMDL distclean
   [ "${icomp}" == "NUOPC"  ] &&  compileNUOPC distclean
  done

  if [ ${CLEAN:-0} -gt 0 ]; then
    exit 0
  fi
fi


##########
# Generate some flags to pass to CMake
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${OPT_TYPE}"
[ -n "${VERBOSE:+1}" ]         && CMAKE_FLAGS+=" -DCMAKE_VERBOSE_MAKEFILE=TRUE"
[ -n "${ADD_CMAKE_FLAGS:+1}" ] && CMAKE_FLAGS+=" ${ADD_CMAKE_FLAGS}"
##########


for icomp in ${COMPONENT}
do
 [ "${icomp}" == "OGCMDL" ] &&  compileOGCMDL compile
 [ "${icomp}" == "NUOPC"  ] &&  compileNUOPC compile
done


exit 0
