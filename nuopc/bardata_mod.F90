!----------------------------------------------------------------
!               M O D U L E   B A R D A T _ M O D
!----------------------------------------------------------------
!> @file bardat_mod.F90
!>
!>
!> @brief
!>
!>
!> @details
!>
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

MODULE Bardat_Mod

  USE MPI
  USE ESMF
  USE NUOPC
  USE NetCDF

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN = 16), PARAMETER, PUBLIC :: MODEL_NAME = 'BARDATA'
  CHARACTER(LEN =  3), PARAMETER, PUBLIC :: MODEL_PRFX = 'bar'

  INTEGER, PARAMETER :: IMISSV = -999999

  PUBLIC :: ReadConfig
  PUBLIC :: InitBarData, ReadBarData
  PUBLIC :: ToLowerCase, ToUpperCase

  PUBLIC :: create_parallel_esmf_grid_from_griddata
  PUBLIC :: construct_griddata_from_netcdf


  ! For reading the NetCDF file
  INTEGER,                          PUBLIC :: TSTRLEN, NT, NX, NY, NYY, NYYY
  CHARACTER(LEN = 16), ALLOCATABLE, PUBLIC :: TIMES(:)
  REAL(ESMF_KIND_R8),  ALLOCATABLE, PUBLIC :: lon(:), lat(:), lonC(:), latC(:), latS(:)
  REAL(ESMF_KIND_R8),  ALLOCATABLE, PUBLIC :: BPGX(:, :, :), BPGY(:, :, :), SigTS(:, :, :), &
                                              MLD(:, :, :), NB(:, :, :), NM(:, :, :)
  REAL(ESMF_KIND_R8),  ALLOCATABLE, PUBLIC :: KE(:, :, :), CDisp(:, :, :), DispX(:, :, :), DispY(:, :, :)


  INTEGER :: ncID
  INTEGER :: timeDimID, strDimID, nxDimID, nyDimID, nyyDimID, nyyyDimID
  INTEGER :: timeDimVAL, strDimVAL, nxDimVAL, nyDimVAL, nyyDimVAL, nyyyDimVAL

  INTEGER :: timeVarID, lonVarID, latVarID, loncVarID, latcVarID, latsVarID
  INTEGER :: bpgxVarID, bpgyVarID, sigtsVarID, mldVarID, nbVarID, nmVarID
  INTEGER :: keVarID, cdispVarID, dispxVarID, dispyVarID


  CHARACTER(LEN = 256)  :: modopt_dir, modopt_nam, modopt_grd
  CHARACTER(LEN = 256)  :: FILE_NAME
  CHARACTER(LEN = 2048) :: info

  ! Type GridData_T for a structured grid
  TYPE, PUBLIC :: GridData_T
    TYPE(ESMF_VM)            :: vm
    INTEGER                  :: maxIndex(2)
    TYPE(ESMF_GridConn_Flag) :: connflagDim1(2)
    TYPE(ESMF_GridConn_Flag) :: connflagDim2(2)
    TYPE(ESMF_CoordSys_Flag) :: coordSys 
  END TYPE GridData_T


  CONTAINS


  !----------------------------------------------------------------
  !  S U B R O U T I N E   I N I T B A R D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Extracts all the static data from the "baroclinic" data file.
  !>
  !> @details
  !>   Examines the input file for valid and/or missing data dimensions
  !>   and data variables. Upon success, extracts all the static data from
  !>   the "baroclinic" data file.
  !>
  !> @param
  !>
  !----------------------------------------------------------------
  SUBROUTINE InitBarData()

    USE NetCDF, ONLY: nf90_open, nf90_close, NF90_NOWRITE, &
                      nf90_inq_dimid, nf90_inquire_dimension, &
                      nf90_inq_varid, nf90_get_var

    IMPLICIT NONE

    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // '_mod:InitBarData)'

    LOGICAL :: THERE
    INTEGER :: rc

    FILE_NAME = TRIM(ADJUSTL(modopt_dir)) // '/' // TRIM(ADJUSTL(modopt_nam))
    PRINT *, ' FILE_NAME  > ', FILE_NAME

    INQUIRE (FILE = FILE_NAME, EXIST = THERE)
    IF (.NOT. THERE) STOP TRIM(ADJUSTL(MODEL_NAME)) // ' NetCDF file does not exist!'

    ncID = 0
    ! Open the file.
    CALL CheckErr(nf90_open(TRIM(FILE_NAME), NF90_NOWRITE, ncID))

    !==============================
    !===== Get all dimension IDs and their values
    CALL CheckErr(nf90_inq_dimid(ncID, 'time',   timeDimID))
    CALL CheckErr(nf90_inq_dimid(ncID, 'strlen', strDimID))
    CALL CheckErr(nf90_inq_dimid(ncID, 'NX',     nxDimID))
    CALL CheckErr(nf90_inq_dimid(ncID, 'NY',     nyDimID))
    CALL CheckErr(nf90_inq_dimid(ncID, 'NYY',    nyyDimID))
    CALL CheckErr(nf90_inq_dimid(ncID, 'NYYY',   nyyyDimID))

    CALL CheckErr(nf90_inquire_dimension(ncID, timeDimID, len = timeDimVAL))
    CALL CheckErr(nf90_inquire_dimension(ncID, strDimID,  len = strDimVAL))
    CALL CheckErr(nf90_inquire_dimension(ncID, nxDimID,   len = nxDimVAL))
    CALL CheckErr(nf90_inquire_dimension(ncID, nyDimID,   len = nyDimVAL))
    CALL CheckErr(nf90_inquire_dimension(ncID, nyyDimID,  len = nyyDimVAL))
    CALL CheckErr(nf90_inquire_dimension(ncID, nyyyDimID, len = nyyyDimVAL))

    !==============================
    !===== Get all variable IDs
    CALL CheckErr(nf90_inq_varid(ncID, 'time', timeVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'lon',  lonVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'lat',  latVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'lonc', loncVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'latc', latcVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'lats', latsVarID))

    CALL CheckErr(nf90_inq_varid(ncID, 'BPGX',  bpgxVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'BPGY',  bpgyVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'SigTS', sigtsVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'MLD',   mldVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'NB',    nbVarID))
    CALL CheckErr(nf90_inq_varid(ncID, 'NM',    nmVarID))

    !CALL CheckErr(nf90_inq_varid(ncID, 'KE',    keVarID))
    !CALL CheckErr(nf90_inq_varid(ncID, 'CDisp', cdispVarID))
    !CALL CheckErr(nf90_inq_varid(ncID, 'DispX', dispxVarID))
    !CALL CheckErr(nf90_inq_varid(ncID, 'DispY', dispyVarID))

    !==============================
    !===== Allocate variable arrays
    IF (.NOT. ALLOCATED(TIMES)) ALLOCATE(TIMES(1:timeDimVAL))
    IF (.NOT. ALLOCATED(lon))   ALLOCATE(lon(1:nxDimVAL))
    IF (.NOT. ALLOCATED(lat))   ALLOCATE(lat(1:nyDimVAL))
    IF (.NOT. ALLOCATED(lonC))  ALLOCATE(lonC(1:nxDimVAL))
    IF (.NOT. ALLOCATED(latC))  ALLOCATE(latC(1:nyyDimVAL))
    IF (.NOT. ALLOCATED(latS))  ALLOCATE(latS(1:nyyyDimVAL))

    !==============================
    !===== Read some variables
    CALL CheckErr(nf90_get_var(ncID, timeVarID, TIMES))
    CALL CheckErr(nf90_get_var(ncID, lonVarID,  lon))
    CALL CheckErr(nf90_get_var(ncID, latVarID,  lat))
    CALL CheckErr(nf90_get_var(ncID, loncVarID, lonC))
    CALL CheckErr(nf90_get_var(ncID, latcVarID, latC))
    CALL CheckErr(nf90_get_var(ncID, latsVarID, latS))

    TSTRLEN = strDimVAL
    NT      = strDimVAL
    NX      = nxDimVAL
    NY      = nyDimVAL
    NYY     = nyyDimVAL
    NYYY    = nyyyDimVAL

    WRITE (info, *) subname, ' --- init ' // TRIM(ADJUSTL(MODEL_NAME)) // ' NetCDF file  --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE InitBarData

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   R E A D B A R D A T A
  !----------------------------------------------------------------
  !> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
  !>
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadBarData(currTime)

    IMPLICIT NONE

    TYPE(ESMF_Time), INTENT(IN)   :: currTime

    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // '_mod:ReadBarData)'

    TYPE(ESMF_Time)         :: dataTime
    TYPE(ESMF_TimeInterval) :: dTime

    INTEGER                 :: start(3), count(3)
    REAL(ESMF_KIND_R8)      :: deltaTime(NT)

    INTEGER                 :: YY, MM, DD, HH, MN, SS
    INTEGER                 :: dtYY, dtMM, dtDD, dtHH, dtMN, dtSS
    INTEGER                 :: iCnt, it, rc

    rc = ESMF_SUCCESS


    DO iCnt = 1, NT
      CALL SplitDateTimeString(PreProcessDateTimeString(TIMES(iCnt)), YY, MM, DD, HH, MN, SS)

      CALL ESMF_TimeSet(dataTime, yy = YY, mm = MM, dd = DD, h = HH, m = MN, s = SS, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__,  &
                             file = __FILE__)) &
       RETURN  ! bail out

       dTime = currTime - dataTime
       CALL ESMF_TimeIntervalGet(dTime, d = dtDD, h = dtHH, m = dtMN, s = dtSS, rc = rc)
       IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                              line = __LINE__,  &
                              file = __FILE__)) &
         RETURN  ! bail out

      deltaTime(iCnt) = dtDD + dtHH / 24.0 + dtMN / 1440.0  + dtSS / 86400.0
    END DO

    it = MINLOC(ABS(deltaTime), 1)

    WRITE (info, *) subname, TRIM(ADJUSTL(MODEL_NAME)) // ' file time > ', TRIM(TIMES(it))
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

    !==============================
    !===== Allocate variable arrays
    IF (.NOT. ALLOCATED(BPGX))  ALLOCATE(BPGX(1:NX, 1:NY, 1))
    IF (.NOT. ALLOCATED(BPGY))  ALLOCATE(BPGY(1:NX, 1:NYY, 1))
    IF (.NOT. ALLOCATED(SigTS)) ALLOCATE(SigTS(1:NX, 1:NY, 1))
    IF (.NOT. ALLOCATED(MLD))   ALLOCATE(MLD(1:NX, 1:NY, 1))
    IF (.NOT. ALLOCATED(NB))    ALLOCATE(NB(1:NX, 1:NY, 1))
    IF (.NOT. ALLOCATED(NM))    ALLOCATE(NM(1:NX, 1:NY, 1))

    !IF (.NOT. ALLOCATED(KE))    ALLOCATE(KE(1:NX, 1:NY, 1))
    !IF (.NOT. ALLOCATED(CDisp)) ALLOCATE(CDisp(1:NX, 1:NYYY, 1))
    !IF (.NOT. ALLOCATED(DispX)) ALLOCATE(DispX(1:NX, 1:NYYY, 1))
    !IF (.NOT. ALLOCATED(DispY)) ALLOCATE(DispY(1:NX, 1:NYYY, 1))

    !==============================
    !===== Read the variables
    start = (/1, 1, it/)
    count = (/NX, NY, 1/)

    BPGX  = 0.0
    BPGY  = 0.0
    SigTS = 0.0
    MLD   = 0.0
    NB    = 0.0
    NM    = 0.0
!    KE    = 0.0
!    CDisp = 0.0
!    DispX = 0.0
!    DispY = 0.0

    CALL CheckErr(nf90_get_var(ncID, bpgxVarID,  BPGX,  start, count))
    CALL CheckErr(nf90_get_var(ncID, bpgyVarID,  BPGY,  start, (/NX, NYY, 1/)))
    CALL CheckErr(nf90_get_var(ncID, sigtsVarID, SigTS, start, count))
    CALL CheckErr(nf90_get_var(ncID, mldVarID,   MLD,   start, count))
    CALL CheckErr(nf90_get_var(ncID, nbVarID,    NB,    start, count))
    CALL CheckErr(nf90_get_var(ncID, nmVarID,    NM,    start, count))

    !CALL CheckErr(nf90_get_var(ncID, keVarID,    KE,    start, count))
    !CALL CheckErr(nf90_get_var(ncID, cdispVarID, CDisp, start, (/NX, NYYY, 1/)))
    !CALL CheckErr(nf90_get_var(ncID, dispxVarID, DispX, start, (/NX, NYYY, 1/)))
    !CALL CheckErr(nf90_get_var(ncID, dispyVarID, DispY, start, (/NX, NYYY, 1/)))

    WRITE (info, *) subname, ' --- read ' // TRIM(ADJUSTL(MODEL_NAME)) // ' NetCDF file  --- '
    !PRINT *, info
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE ReadBarData

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   C H E C K E R R
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Checks for error in a NetCDF procedure call.
  !>
  !> @details
  !>
  !>
  !> @param[in]
  !>   status    The status 
  !> @param[in]
  !>   message   A message string to be included in the output message
  !>             Usually the name of the calling procedure (optional)
  !> @param[in]
  !>   ncID      The ID of the NetCDF file. If present CheckErr will try
  !>             to close the file (optional)
  !> @param[in]
  !>   varID     The ID of the variable in the NetCDF file. If present CheckErr
  !>             will output the ID of the variable (optional)
  !>
  !----------------------------------------------------------------
  SUBROUTINE CheckErr(status, message, ncID, varID)

    USE NetCDF, ONLY: nf90_strerror, nf90_noerr, nf90_close

    IMPLICIT NONE

    ! (should include name of calling procedure)
    INTEGER, INTENT(IN)                      :: status
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: message

    ! (Provide this argument if you want "CheckErr" to try to close
    ! the file.)
    INTEGER, INTENT(IN), OPTIONAL            :: ncID
    INTEGER, INTENT(IN), OPTIONAL            :: varID

    ! Variable local to the procedure:
    INTEGER status_close

    IF (status /= nf90_noerr) THEN
      IF (PRESENT(message)) PRINT *, message, ":"
      IF (PRESENT(varID)) PRINT *, "varid = ", varID
      PRINT *, TRIM(nf90_strerror(status))
      IF (PRESENT(ncID)) THEN
        ! Try to close, to leave the file in a consistent state:
        status_close = nf90_close(ncID)
        IF (status_close /= nf90_noerr) THEN
          PRINT *, "nf90_close:"
          PRINT *, TRIM(nf90_strerror(status_close))
        END IF
      END IF
      STOP
    END IF

  END SUBROUTINE CheckErr

!================================================================================

  !-----------------------------------------------------------------------
  !- Sub !!!????
  !-----------------------------------------------------------------------
  !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
  !> this function extracts the scalars and arrays required for construction of a
  !> MeshData_T object.
  !> After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
  !> or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
  !> @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
  !> and \c peCount of the \c MPI_Communicator.
  !> @param global_fort14_dir This is the directory path (relative to the executable
  !> or an absolute path) which contains the global \c fort.14 file (not the fort.14
  !> after decomposition).
  !> @param the_data This is the output MeshData_T object.
  !>

  !> \details As the name of this function suggests, this funciton creates a parallel
  !> ESMF_Mesh from MeshData_T object. This function should be called collectively by
  !> all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
  !> should be called prior to calling this function.
  !> \param theData This the input MeshData_T object.
  !> \param out_esmf_grid This is the ouput ESMF_Mesh object.
  SUBROUTINE create_parallel_esmf_grid_from_griddata(varX, varY, theData, out_esmf_grid)

    IMPLICIT NONE

    TYPE(ESMF_Grid),  INTENT(OUT)                 :: out_esmf_grid
    REAL(ESMF_KIND_R8), INTENT(IN)                :: varX(:), varY(:)
    TYPE(GridData_T), INTENT(IN)                  :: theData
    INTEGER                                       :: rc

    ! THis fuction is 31.6.9. Create a Grid with user set edge connections and
    ! a regular distribution
    ! maxIndex specifies the dimension of grid: the upper extent of the grid
    ! array; 
    ! connflagDim1/2 = ESMF_GRIDCONN_NONE (default) without presence.
    ! coordSys = ESMF_COORDSYS_SPH_DEG (default) if not specified.
    ! But where did the coordinate information comes in??
    ! =============== another method to create grid code below: =================!
    ! out_esmf_grid=ESMF_GridCreate(maxIndex=theData%maxIndex, connflagDim1=theData%connflagDim1, &
    !           connflagDim2=theData%connflagDim2, coordSys=theData%coordSys,  rc=rc)

    ! call a function to put in the coordinate for this grid:
    out_esmf_grid = CreateGrid_ModelGrid(theData%maxIndex(1), theData%maxIndex(2), varX, varY, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

  END SUBROUTINE create_parallel_esmf_grid_from_griddata

!================================================================================

  SUBROUTINE construct_griddata_from_netcdf(nDimX, nDimY, theData)

    IMPLICIT NONE

    INTEGER, INTENT(IN)             :: nDimX, nDimY
    TYPE(GridData_T), INTENT(INOUT) :: theData

    theData%maxIndex        = (/nDimX, nDimY/)
    theData%connflagDim1(1) = ESMF_GRIDCONN_NONE
    theData%connflagDim1(2) = ESMF_GRIDCONN_NONE
    theData%connflagDim2(1) = ESMF_GRIDCONN_NONE
    theData%connflagDim2(2) = ESMF_GRIDCONN_NONE
    theData%coordSys        = ESMF_COORDSYS_SPH_DEG ! spherical grid in degrees

  END SUBROUTINE

!================================================================================

  FUNCTION CreateGrid_ModelGrid(nVarX, nVarY, varX, varY, rc)

    IMPLICIT NONE

    TYPE(ESMF_Grid)                :: CreateGrid_ModelGrid

    INTEGER, INTENT(IN)            :: nVarX, nVarY
    REAL(ESMF_KIND_R8), INTENT(IN) :: varX(:), varY(:)
    INTEGER, INTENT(OUT), OPTIONAL :: rc

    REAL(ESMF_KIND_R8), POINTER    :: coordX(:, :), coordY(:, :)
    INTEGER                        :: i, j

    CHARACTER(LEN = *), PARAMETER  :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // &
                                                '_mod: CreateGrid_ModelGrid)'

    rc = ESMF_SUCCESS

    CreateGrid_ModelGrid = ESMF_GridCreateNoPeriDim(name = "ModelGrid", &
                           minIndex  = (/1, 1/), &
                           maxIndex  = (/nVarX, nVarY/), &
                           indexflag = ESMF_INDEX_GLOBAL, &
                           rc=rc)

    ! Add coordinates
    call ESMF_GridAddCoord(CreateGrid_ModelGrid, &
                           staggerloc = ESMF_STAGGERLOC_CENTER, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        LINE = __LINE__, &
        FILE = __FILE__)) &
        RETURN  ! bail out

    CALL ESMF_GridGetCoord(CreateGrid_ModelGrid, coordDim = 1,  &
                           staggerloc = ESMF_STAGGERLOC_CENTER, &
                           farrayPtr = coordX, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        LINE = __LINE__, &
        FILE = __FILE__)) &
        RETURN  ! bail out

    CALL ESMF_GridGetCoord(CreateGrid_ModelGrid, coordDim = 2, &
                           staggerloc = ESMF_STAGGERLOC_CENTER, &
                           farrayPtr = coordY, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        LINE = __LINE__, &
        FILE = __FILE__)) &
        RETURN  ! bail out

    ! Set coordinates
    DO i = 1, nVarX
      DO j = 1, nVarY
        coordX(i, j) = varX(i)
        coordY(i, j) = varY(j)
      END DO
    END DO

  WRITE(info,*) subname,' --- ' // TRIM(ADJUSTL(MODEL_NAME)) // ' CreateGrid_ModelGrid completed --- '
  CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)
     
  END FUNCTION CreateGrid_ModelGrid

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   R E A D C O N F I G
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   This subroutine reads the program's main control file.
  !>
  !> @details
  !>   Reads the control file of the program and it is repeatedly calling GetLineRecord
  !>   to process each line in the file. Upon successful processing of the line,
  !>   it sets the relevant program parameters and variables. This subroutine is called
  !>   first as it is required by the subsequent model run. \n
  !>   The control file (default filename pahm_control.in) contains all required
  !>   settings (user configured) required to run the program. Most of the settings
  !>   have default values, in case the user hasn't supplied a value.
  !>
  !> @param[in]
  !>   inpFile   The full pathname of the input file
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadConfig(OptPrefix, ConfFile)

    IMPLICIT NONE

    ! prefix of the field name in the config file
    CHARACTER(LEN = *), INTENT(IN)                   :: OptPrefix
    ! config file name, default config.rc
    CHARACTER(ESMF_MAXPATHLEN), INTENT(IN), OPTIONAL :: ConfFile

    ! Local variables
    CHARACTER(ESMF_MAXPATHLEN) :: fname ! config file name
    TYPE(ESMF_Config)          :: cf    ! the Config itself
    INTEGER                    :: rc

    rc = ESMF_SUCCESS

    !Initiate reading resource file
    cf = ESMF_ConfigCreate(rc = rc)  ! Create the empty Config
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! Name the Resource File
    IF (PRESENT(ConfFile)) THEN
      fname = TRIM(ADJUSTL(ConfFile))
    ELSE
      fname = "config.rc"
    END IF

    CALL ESMF_ConfigLoadFile(cf, fname, rc = rc) ! Load the Resource File
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_ConfigGetAttribute(cf, modopt_dir,                              &
                                 label = TRIM(ADJUSTL(OptPrefix)) // '_dir:', &
                                 default = './',                              &
                                 rc = rc)

    CALL ESMF_ConfigGetAttribute(cf, modopt_nam,                              &
                                 label = TRIM(ADJUSTL(OptPrefix)) // '_nam:', &
                                 default = 'moddata.nc',                      &
                                 rc = rc)

    CALL ESMF_ConfigGetAttribute(cf, modopt_grd,                              &
                                 label = TRIM(ADJUSTL(OptPrefix)) // '_grd:', &
                                 default = 'mod_grd.nc',                      &
                                 rc = rc)

    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! Destroy the Config
    CALL ESMF_ConfigDestroy(cf, rc = rc)

  END SUBROUTINE ReadConfig

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   T O  L O W E R  C A S E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert a string to lower-case.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpString   The input string
  !>
  !> @return
  !>   outString: The ouput string in lower case
  !>
  !----------------------------------------------------------------
  PURE FUNCTION ToLowerCase(inpString) RESULT(outString)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: inpString

    INTEGER, PARAMETER        :: DUC = ICHAR('A') - ICHAR('a')
    CHARACTER(LEN(inpString)) :: outString
    CHARACTER                 :: ch
    INTEGER                   :: i

    DO i = 1, LEN(inpString)
      ch = inpString(i:i)
      IF ((ch >= 'A') .AND. (ch <= 'Z')) ch = CHAR(ICHAR(ch) - DUC)
      outString(i:i) = ch
    END DO

    RETURN

  END FUNCTION ToLowerCase

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   T O  U P P E R  C A S E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert a string to upper-case.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpString   The input string
  !>
  !> @return
  !>   outString: The ouput string in upper case
  !>
  !----------------------------------------------------------------
  PURE FUNCTION ToUpperCase(inpString) RESULT(outString)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: inpString

    INTEGER, PARAMETER        :: DUC = ICHAR('A') - ICHAR('a')
    CHARACTER(LEN(inpString)) :: outString
    CHARACTER                 :: ch
    INTEGER                   :: i

    DO i = 1, LEN(inpString)
      ch = inpString(i:i)
      IF ((ch >= 'a') .AND. (ch <= 'z')) ch = CHAR(ICHAR(ch) + DUC)
      outString(i:i) = ch
    END DO

    RETURN

  END FUNCTION ToUpperCase

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S P L I T  D A T E  T I M E  S T R I N G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Splits a date string into components.
  !>
  !> @details
  !>   This subroutine splits the string inDate (YYYYMMDDhhmmss) in six integers that is,
  !>   "iYear (YYYY)", "iMonth (MM)", "iDay (DD)", "iHour (hh)", "iMin (mm)" and "iSec (ss)".
  !>
  !> @param[in]
  !>   inDateTime  The input date string: YYYYMMDDhhmmss
  !> @param[out]
  !>   iYear       The year (YYYY, integer, 1582 <= YYYY, output)
  !> @param[out]
  !>   iMonth      The month of the year (MM, integer, 1 <= MM <=12, output)
  !> @param[out]
  !>   iDay        The day of the month (DD, integer, 1 <= DD <=31, output)
  !> @param[out]
  !>   iHour       The hour of the day (hh, integer, 0 <= hh <= 23, output)
  !> @param[out]
  !>   iMin        The minute of the hour (mm, integer, 0 <= mm <= 59, output)
  !> @param[out]
  !>   iSec        The second of the minute (ss, integer, 0 <= ss <= 59, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE SplitDateTimeString(inDateTime, iYear, iMonth, iDay, iHour, iMin, iSec)

    IMPLICIT NONE

    ! Global Variables
    CHARACTER(LEN=*), INTENT(IN)   :: inDateTime
    INTEGER, INTENT(OUT)           :: iYear, iMonth, iDay, iHour, iMin, iSec

    ! Local Variables
    CHARACTER(LEN=LEN(inDateTime)) :: tmpDateStr
    INTEGER                        :: errIO

    !----- START CALCULATIONS -----

    tmpDateStr = PreProcessDateTimeString(inDateTime)

    IF (TRIM(tmpDateStr) == '') THEN
      iYear  = IMISSV
      iMonth = 0
      iDay   = 0
      iHour  = 0
      iMin   = 0
      iSec   = 0

      RETURN
    END IF

    READ(tmpDateStr(1:4), '(I4.4)', IOSTAT=errIO) iYear
      IF ((errIO /= 0) .OR. (iYear < 1582)) iYear = IMISSV

    READ(tmpDateStr(5:6), '(I2.2)', IOSTAT=errIO) iMonth
      IF ((errIO /= 0) .OR. (iMonth < 1) .OR. (iMonth > 12)) iMonth = 0

    READ(tmpDateStr(7:8), '(I2.2)', IOSTAT=errIO) iDay
      IF ((errIO /= 0) .OR. (iDay < 0) .OR. (iDay > MonthDays(iYear, iMonth))) iDay = 0

    READ(tmpDateStr(9:10), '(I2.2)', IOSTAT=errIO) iHour
      IF ((errIO /= 0) .OR. (iHour < 0) .OR. (iHour >= 23)) iHour = 0

    READ(tmpDateStr(11:12), '(I2.2)', IOSTAT=errIO) iMin
      IF ((errIO /= 0) .OR. (iMin < 0) .OR. (iMin >= 60)) iMin = 0

    READ(tmpDateStr(13:14), '(I2.2)', IOSTAT=errIO) iSec
      IF ((errIO /= 0) .OR. (iSec < 0) .OR. (iSec >= 60)) iSec = 0

  END SUBROUTINE SplitDateTimeString

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   P R E  P R O C E S S  D A T E  T I M E  S T R I N G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Pre-processes an arbitrary date string.
  !>
  !> @details
  !>   This function returns a date/time string in the format YYYYMMDDhhmmss by
  !>   removing all non-numeric characters from the string.
  !>
  !> @param[in]
  !>   inDateTime  The input date string
  !>
  !> @return
  !>   myValOut    The string datetime as an integer in the form: YYYYMMDDhhmmss
  !>
  !----------------------------------------------------------------
  FUNCTION PreProcessDateTimeString(inDateTime) Result(myValOut)

    IMPLICIT NONE

    ! Global Variables
    CHARACTER(LEN=*), INTENT(IN)   :: inDateTime
    CHARACTER(LEN=LEN(inDateTime)) :: myValOut

    ! Local Variables
    CHARACTER(LEN=1)               :: c
    INTEGER                        :: i, iPos

    !----- START CALCULATIONS -----

    myValOut = ' '
    iPos = 1

    DO i = 1, LEN(inDateTime)
      c = inDateTime(i:i)
      IF ((48 <= ichar(c)) .AND. (ichar(c) <= 57)) THEN
        myValOut(iPos:iPos) = c
        iPos = iPos + 1
      ENDIF
    END DO

    RETURN

  END FUNCTION PreProcessDateTimeString

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   L E A P  Y E A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Checks for a leap year.
  !>
  !> @details
  !>   This function tries to determine if a Gregorian year (>= 1582) 
  !>   is a leap year or not.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !>
  !> @return
  !>   myVal .TRUE. if it is a leap year or .FALSE. otherwise
  !>
  !----------------------------------------------------------------
  LOGICAL FUNCTION LeapYear(iYear) RESULT(myVal)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iYear

    !----- START CALCULATIONS -----

    IF (iYear < 1582) THEN
      myVal = .FALSE.

      RETURN
    END IF

    ! ADCIRC uses the construct leap = (iYear / 4) * 4 == iYear
    ! to determine if a year is a leap year. This produces wrong
    ! results, example while 1700, 1900, 2100 are not leap years,
    ! the above construct determines that these years are leap years.
    ! Needs to be fixed.

    IF ((MOD(iYear, 100) /= 0) .AND. (MOD(iYear, 4) == 0)) THEN
      myVal = .TRUE.
    ELSE IF (MOD(iYear, 400) == 0) THEN
      myVal = .TRUE.
    ELSE
      myVal = .FALSE.
    END IF

    RETURN
  END FUNCTION LeapYear

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   Y E A R  D A Y S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the days of the year.
  !>
  !> @details
  !>   This function calculates the number of calendar days of a 
  !>   Gregorian year (>= 1582).
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !>
  !> @return
  !>   myVal     The days of the year (365 or 366)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION YearDays(iYear) RESULT(myVal)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iYear

    !----- START CALCULATIONS -----

    myVal = 365
    IF (LeapYear(iYear)) myVal = 366

    RETURN
  END FUNCTION YearDays

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   M O N T H  D A Y S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the days in the month of the year.
  !>
  !> @details
  !>   This function calculates the number of calendar days in a month
  !>   of a Gregorian year (>= 1582). In case of an error, the value
  !>   IMISSV (-999999) is returned.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <= 12)
  !>
  !> @return
  !>   myVal     The days of the month
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION MonthDays(iYear, iMonth) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN) :: iYear, iMonth

    ! Local variables
    INTEGER :: leap, monLen(12, 2)

    !----- START CALCULATIONS -----

    IF ((iYear < 1582) .OR. (iMonth < 1) .OR. (iMonth > 12)) THEN
      myVal = IMISSV

      RETURN
    END IF

    ! Initialize lenghts of months:
    monLen = RESHAPE((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,     &
                        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),  &
                     (/ 12, 2 /))

    leap = 1
    IF (LeapYear(iYear)) leap = 2

    myVal = monLen(iMonth, leap)

    RETURN
  END FUNCTION MonthDays

!================================================================================

END MODULE Bardat_Mod
