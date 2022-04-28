#!/usr/bin/env bash

if [ -e /.dockerenv ]; then
  set -eo pipefail
  set  -o errtrace

  trap "ERROR: There was an error, details to follow" ERR
fi

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
#  readonly scrDIR="$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}" )" )" && pwd -P)"
  readonly scrNAME="$( grealpath -s "${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}" )"
  readonly scrDIR="$(cd "$(dirname "${scrNAME}" )" && pwd -P)"
else
#  readonly scrDIR="$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}" )" )" && pwd -P)"
  readonly scrNAME="$( realpath -s "${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}" )"
  readonly scrDIR="$(cd "$(dirname "${scrNAME}" )" && pwd -P)"
fi

lst="${scrDIR:+${scrDIR}/}functions_build ${scrDIR:+${scrDIR}/}scripts/functions_build functions_build"
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
  [ $? -ne 0 ] && exit 1
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

for icomp in ${COMPONENT}
do
  [ "${icomp}" == "OGCMDL" ] && compileOGCMDL compile
  [ "${icomp}" == "NUOPC"  ] &&  compileNUOPC compile
done

exit 0
