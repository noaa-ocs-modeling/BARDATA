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

  !-----------------------------------------------------------------------------
  ! ADCIRC mesh utility
  !-----------------------------------------------------------------------------
  USE MPI
  USE ESMF
  USE NUOPC
  USE NetCDF

  !USE MESH   , only: np,ne,nm,slam,sfea
  !USE GLOBAL , only: IMAP_EL_LG,NODES_LG
  !USE GLOBAL , only: ETA2, UU2, VV2  ! Export water level and velocity fileds to bardata model
  !USE GLOBAL,  ONLY: RSNX2, RSNY2    ! Import bardata 2D forces from bardata model
  !USE SIZES  , only: ROOTDIR

  IMPLICIT NONE

  CHARACTER(LEN = 280)  :: modopt_dir, modopt_nam, modopt_grd
  CHARACTER(LEN = 280)  :: FILE_NAME
  CHARACTER(LEN = 2048) :: info

  ! info for reading bardata netcdf file
  INTEGER               :: nnode, nelem, ntime, noel
  REAL(ESMF_KIND_R8), ALLOCATABLE     :: LONS(:), LATS(:), TIMES(:)
  INTEGER, ALLOCATABLE     :: TRI(:, :)
  INTEGER, ALLOCATABLE     :: TRID(:, :)
  REAL(ESMF_KIND_R8), ALLOCATABLE     :: UWND(:, :), VWND(:, :), PRES(:, :)
  !netcdf vars
  INTEGER :: ncid, NOD_dimid, rec_dimid, ELM_dimid, NOE_dimid
  INTEGER :: LON_varid, LAT_varid, rec_varid, tri_varid
  INTEGER :: UWND_varid, VWND_varid, PRES_varid

  !> \author Ali Samii - 2016
    !! See: https://github.com/samiiali
    !! \brief This object stores the data required for construction of a parallel or serial
    !! ESMF_Mesh from <tt>fort.14, fort.18, partmesh.txt</tt> files.
    !!
  TYPE meshdata
    !> \details vm is an ESMF_VM object.  ESMF_VM is just an ESMF virtual machine class,
      !! which we will use to get the data about the local PE and PE count.
    TYPE(ESMF_VM)                      :: vm
    !> \details This array contains the node coordinates of the mesh. For
      !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
      !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
    REAL(ESMF_KIND_R8), ALLOCATABLE    :: NdCoords(:)
    !> \details This array contains the elevation of different nodes of the mesh
    REAL(ESMF_KIND_R8), ALLOCATABLE    :: bathymetry(:)
    !> \details Number of nodes present in the current PE. This is different from the
      !! number of nodes owned by this PE (cf. NumOwnedNd)
    INTEGER(ESMF_KIND_I4)              :: NumNd
    !> \details Number of nodes owned by this PE. This is different from the number of
      !! nodes present in the current PE (cf. NumNd)
    INTEGER(ESMF_KIND_I4)              :: NumOwnedNd
    !> \details Number of elements in the current PE. This includes ghost elements and
      !! owned elements. However, we do not bother to distinguish between owned
      !! element and present element (as we did for the nodes).
    INTEGER(ESMF_KIND_I4)              :: NumEl
    !> \details Number of nodes of each element, which is simply three in 2D ADCIRC.
    INTEGER(ESMF_KIND_I4)              :: NumND_per_El
    !> \details Global node numbers of the nodes which are present in the current PE.
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: NdIDs(:)
    !> \details Global element numbers which are present in the current PE.
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: ElIDs(:)
    !> \details The element connectivity array, for the present elements in the current PE.
      !! The node numbers are the local numbers of the present nodes. All the element
      !! connectivities are arranged in this one-dimensional array.
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: ElConnect(:)
    !> \details The number of the PE's which own each of the nodes present this PE.
      !! This number is zero-based.
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: NdOwners(:)
    !> \details An array containing the element types, which are all triangles in our
      !! application.
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: ElTypes(:)
    !> \details This is an array, which maps the indices of the owned nodes to the indices of the present
      !! nodes. For example, assume we are on <tt>PE = 1</tt>, and we have four nodes present, and the
      !! first and third nodes belong to <tt>PE = 0</tt>. So we have:
      !! \code
      !! NumNd = 4
      !! NumOwnedNd = 2
      !! NdOwners = (/0, 1, 0, 1/)
      !! NdIDs = (/2, 3, 5, 6/)
      !! owned_to_present = (/2, 4/)    <-- Because the first node owned by this PE is actually
      !!                                    the second node present on this PE, and so on.
      !! \endcode
    INTEGER(ESMF_KIND_I4), ALLOCATABLE :: owned_to_present_nodes(:)
  END TYPE meshdata


  CONTAINS


!-----------------------------------------------------------------------
!- Sub !!!????
!-----------------------------------------------------------------------
  SUBROUTINE NCDF_InitBarData()

    IMPLICIT NONE

    CHARACTER(LEN = *), PARAMETER :: NOD_NAME = "node"
    CHARACTER(LEN = *), PARAMETER :: NOE_NAME = "noel"
    CHARACTER(LEN = *), PARAMETER :: ELM_NAME = "element"
    CHARACTER(LEN = *), PARAMETER :: LAT_NAME = "latitude"
    CHARACTER(LEN = *), PARAMETER :: LON_NAME = "longitude"
    CHARACTER(LEN = *), PARAMETER :: REC_NAME = "time"
    CHARACTER(LEN = *), PARAMETER :: UWND_NAME = "uwnd"
    CHARACTER(LEN = *), PARAMETER :: VWND_NAME = "vwnd"
    CHARACTER(LEN = *), PARAMETER :: PRES_NAME = "P"
    CHARACTER(LEN = *), PARAMETER :: TRI_NAME = "tri"

    CHARACTER(LEN = 140)          :: units
    CHARACTER(LEN = *), PARAMETER :: subname = '(bardata_mod:NCDF_InitBarData)'

    logical :: THERE
    INTEGER :: lat, lon, i, iret, rc, num

    FILE_NAME = TRIM(modopt_dir)//'/'//TRIM(modopt_nam)
    PRINT *, ' FILE_NAME  > ', FILE_NAME
    INQUIRE (FILE = FILE_NAME, EXIST = THERE)
    IF (.NOT. THERE) stop 'BARDATA netcdf grdfile does not exist!'

    ncid = 0
    ! Open the file.
    CALL CheckErr(nf90_open(TRIM(FILE_NAME), NF90_NOWRITE, ncid))

    ! Get ID of unlimited dimension
    !CALL CheckErr( nf90_inquire(ncid, unlimitedDimId = rec_dimid) )

    ! Get ID of limited dimension
    CALL CheckErr(nf90_inq_dimid(ncid, REC_NAME, rec_dimid))
    CALL CheckErr(nf90_inq_dimid(ncid, NOD_NAME, NOD_dimid))
    CALL CheckErr(nf90_inq_dimid(ncid, ELM_NAME, ELM_dimid))
    CALL CheckErr(nf90_inq_dimid(ncid, NOE_NAME, NOE_dimid))

    ! How many values of "nodes" are there?
    CALL CheckErr(nf90_inquire_dimension(ncid, NOD_dimid, len = nnode))
    CALL CheckErr(nf90_inquire_dimension(ncid, ELM_dimid, len = nelem))
    CALL CheckErr(nf90_inquire_dimension(ncid, NOE_dimid, len = noel))
    ! What is the name of the unlimited dimension, how many records are there?
    CALL CheckErr(nf90_inquire_dimension(ncid, rec_dimid, len = ntime))

    !PRINT *,  ' nelem  > ',nelem , ' noel  > ' ,noel,  ' ntime > ',ntime

    ! Get the varids of the pressure and temperature netCDF variables.
    CALL CheckErr(nf90_inq_varid(ncid, LAT_NAME, LAT_varid))
    CALL CheckErr(nf90_inq_varid(ncid, LON_NAME, LON_varid))
    CALL CheckErr(nf90_inq_varid(ncid, REC_NAME, rec_varid))
    CALL CheckErr(nf90_inq_varid(ncid, UWND_NAME, UWND_varid))
    CALL CheckErr(nf90_inq_varid(ncid, VWND_NAME, VWND_varid))
    CALL CheckErr(nf90_inq_varid(ncid, PRES_NAME, PRES_varid))
    CALL CheckErr(nf90_inq_varid(ncid, TRI_NAME, TRI_varid))

    !allocate vars
    IF (.NOT. allocated(LATS)) allocate (LATS(1:nnode))
    IF (.NOT. allocated(LONS)) allocate (LONS(1:nnode))
    IF (.NOT. allocated(TIMES)) allocate (TIMES(1:ntime))
    IF (.NOT. allocated(TRI)) allocate (TRI(1:noel, 1:nelem))
    IF (.NOT. allocated(TRID)) allocate (TRID(1:noel, 1:nelem))
    ! read vars
    CALL CheckErr(nf90_get_var(ncid, LAT_varid, LATS))
    CALL CheckErr(nf90_get_var(ncid, LON_varid, LONS))
    CALL CheckErr(nf90_get_var(ncid, rec_varid, TIMES))
    !CALL CheckErr(nf90_get_var(ncid, UWND_varid, UWND  ))
    !TODO: Why the order is other way???? Might change the whole forcing fields!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT <<<<<
    ! plot input and out put to be sure we are not scrambling the data. the same for HWRF netcdf file
    CALL CheckErr(nf90_get_var(ncid, TRI_varid, TRI, start = (/1, 1/), count = (/noel, nelem/)))
    !TRI = int( TRID )

    !DO num = 1,10
    !    PRINT *,  "TRI", TRI(1,num), TRI(2,num), TRI(3,num)
    !END DO

    WRITE (info, *) subname, ' --- init bardata netcdf file  --- '
    !PRINT *, info
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   N C D F _ R E A D B A R D A T A
  !----------------------------------------------------------------
  !> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
  !>
  !> Extracts the mesh data from the mesh  file using PaHM's ReadMesh()
  !! subroutine and populates the MeshData object "theData" accordingly.
  !!
  !----------------------------------------------------------------
  SUBROUTINE NCDF_ReadBarData(currTime)

    IMPLICIT NONE

    TYPE(ESMF_Time), intent(in)     :: currTime
    TYPE(ESMF_Time)                :: refTime
    TYPE(ESMF_TimeInterval)        :: dTime

    CHARACTER(LEN = 140)          :: units
    CHARACTER(LEN = *), PARAMETER :: subname = '(bardata_mod:NCDF_ReadBarData)'
    INTEGER, PARAMETER :: NDIMS = 2
    INTEGER    :: start(NDIMS), count(NDIMS)
    logical    :: THERE
    REAL       :: delta_d_all(ntime), delta_d_ref
    !INTEGER   :: dimids(NDIMS)

    character  :: c1, c2, c3, c4, c5, c6, c7
    INTEGER    :: yy, mm, dd, hh, min, ss
    INTEGER    :: d_d, d_h, d_m, d_s
    INTEGER    :: lat, lon, it, iret, rc

    rc = ESMF_SUCCESS

    !units = "days since 1990-01-01 00:00:00"
    CALL CheckErr(nf90_get_att(ncid, rec_varid, 'units', units))
    READ (units, '(a4,a7,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)', iostat = iret) &
      c1, c2, yy, c3, mm, c4, dd, c5, hh, c6, min, c7, ss

    IF (iret .ne. 0) THEN
      PRINT *, 'Fatal error: A non valid time units string was provided'
      STOP 'bardata_mod: NCDF_ReadBarData'
    END IF

    CALL ESMF_TimeSet(refTime, yy = yy, mm = mm, dd = dd, h = hh, m = min, s = ss, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dTime = currTime - refTime
    CALL ESMF_TimeIntervalGet(dTime, d = d_d, h = d_h, m = d_m, s = d_s, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    delta_d_ref = d_d + d_h / 24.0 + d_m / (24.0 * 60.0) + d_s / (24.0 * 3600.0)

    DO it = 1, ntime
      delta_d_all(it) = delta_d_ref - TIMES(it)
    END DO

    it = minloc(abs(delta_d_all), dim = 1)

    IF (abs(delta_d_all(it)) > 7200.) THEN
      WRITE (info, *) subname, ' --- STOP BARData: Time dif is GT 2 hours ---  '
      CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)
      STOP ' --- STOP BARData: Time dif is GT 2 hours ---  '
    END IF

    !it = 1
    !PRINT *, 'bardata file it index > ',it

    WRITE (info, *) subname, 'bardata file it index > ', it
    !PRINT *, info
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

    CALL ESMF_TimePrint(refTime, preString = "BARDATA refTime =  ", rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !WRITE(info,*) ' Read BARData netcdf:',it,
    !CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

    !alocate vars
    IF (.NOT. allocated(UWND)) allocate (UWND(1:nnode, 1))
    IF (.NOT. allocated(VWND)) allocate (VWND(1:nnode, 1))
    IF (.NOT. allocated(PRES)) allocate (PRES(1:nnode, 1))

    start = (/1, it/)
    count = (/nnode, 1/)  !for some reason the order here is otherway around?!

    !PRINT *, start+count
    !PRINT *,size(UWND(ntime,:))
    CALL CheckErr(nf90_get_var(ncid, UWND_varid, UWND, start, count))
    CALL CheckErr(nf90_get_var(ncid, VWND_varid, VWND, start, count))
    CALL CheckErr(nf90_get_var(ncid, PRES_varid, PRES, start, count))

    !PRINT *,FILE_NAME , '   HARD CODED for NOWWWW>>>>>     Time index from bardata file is > ', it, UWND(1:10,1)
    WRITE (info, *) subname, ' --- read BARData netcdf file  --- '
    !PRINT *, info
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE
  !-----------------------------------------------------------------------
  !- Sub !!!????
  !-----------------------------------------------------------------------

  SUBROUTINE construct_meshdata_from_netcdf(the_data)

    IMPLICIT NONE

    INTEGER                               :: i1
    INTEGER, PARAMETER                    :: dim1 = 2, spacedim = 2, NumND_per_El = 3
    TYPE(meshdata), intent(inout)         :: the_data
    the_data % NumEl = nelem
    the_data % NumNd = nnode
    allocate (the_data % NdIDs(the_data % NumNd))
    allocate (the_data % ElIDs(the_data % NumEl))
    allocate (the_data % NdCoords(dim1 * the_data % NumNd))
    allocate (the_data % bathymetry(the_data % NumNd))
    allocate (the_data % ElConnect(NumND_per_El * the_data % NumEl))
    allocate (the_data % NdOwners(the_data % NumNd))
    allocate (the_data % ElTypes(the_data % NumEl))
    allocate (the_data % owned_to_present_nodes(the_data % NumNd))

    DO i1 = 1, the_data % NumNd, 1
      the_data % NdIDs(i1) = i1
      the_data % NdCoords((i1 - 1) * dim1 + 1) = LONS(i1)
      the_data % NdCoords((i1 - 1) * dim1 + 2) = LATS(i1)
    END DO
    DO i1 = 1, the_data % NumEl, 1
      the_data % ElIDs(i1) = i1
      the_data % ElConnect((i1 - 1) * NumND_per_El + 1) = TRI(1, i1)
      the_data % ElConnect((i1 - 1) * NumND_per_El + 2) = TRI(2, i1)
      the_data % ElConnect((i1 - 1) * NumND_per_El + 3) = TRI(3, i1)
    END DO
    !We have only one node therefore:
    the_data % NdOwners = 0                  !process 0 owns all the nodes
    the_data % NumOwnedND = the_data % NumNd   !number of nodes = number of owned nodes
    the_data % owned_to_present_nodes = the_data % NdIDs

    the_data % ElTypes = ESMF_MESHELEMTYPE_TRI

    close (14)
  END SUBROUTINE

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
  !>   ncid      The ID of the NetCDF file. If present CheckErr will try
  !>             to close the file (optional)
  !> @param[in]
  !>   varid     The ID of the variable in the NetCDF file. If present CheckErr
  !>             will output the ID of the variable (optional)
  !>
  !----------------------------------------------------------------
  SUBROUTINE CheckErr(status, message, ncid, varid)

    USE NetCDF, ONLY: nf90_strerror, nf90_noerr, nf90_close

    IMPLICIT NONE

    ! (should include name of calling procedure)
    INTEGER, INTENT(IN)                      :: status
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: message

    ! (Provide this argument if you want "CheckErr" to try to close
    ! the file.)
    INTEGER, INTENT(IN), OPTIONAL            :: ncid
    INTEGER, INTENT(IN), OPTIONAL            :: varid

    ! Variable local to the procedure:
    INTEGER status_close

    IF (status /= nf90_noerr) THEN
      IF (PRESENT(message)) PRINT *, message, ":"
      IF (PRESENT(varid)) PRINT *, "varid = ", varid
      PRINT *, TRIM(nf90_strerror(status))
      IF (PRESENT(ncid)) THEN
        ! Try to close, to leave the file in a consistent state:
        status_close = nf90_close(ncid)
        ! (do not call "nf95_close", we do not want to recurse)
        IF (status_close /= nf90_noerr) THEN
          PRINT *, "nf90_close:"
          PRINT *, TRIM(nf90_strerror(status_close))
        END IF
      END IF
      STOP 1
    END IF

  END SUBROUTINE CheckErr

!================================================================================

  !-----------------------------------------------------------------------
  !- Sub !!!????
  !-----------------------------------------------------------------------
  !> \author Ali Samii - 2016
    !! See: https://github.com/samiiali
  !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
    !! this function extracts the scalars and arrays required for construction of a
    !! meshdata object.
    !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
    !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
    !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
    !! and \c peCount of the \c MPI_Communicator.
    !! @param global_fort14_dir This is the directory path (relative to the executable
    !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
    !! after decomposition).
    !! @param the_data This is the output meshdata object.
    !!

  !> \details As the name of this function suggests, this funciton creates a parallel
    !! ESMF_Mesh from meshdata object. This function should be called collectively by
    !! all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
    !! should be called prior to calling this function.
    !! \param the_data This the input meshdata object.
    !! \param out_esmf_mesh This is the ouput ESMF_Mesh object.
  SUBROUTINE create_parallel_esmf_mesh_from_meshdata(the_data, out_esmf_mesh)

    IMPLICIT NONE

    TYPE(ESMF_Mesh), intent(out)                  :: out_esmf_mesh
    TYPE(meshdata), intent(in)                    :: the_data
    INTEGER, PARAMETER                            :: dim1 = 2, spacedim = 2, NumND_per_El = 3
    INTEGER                                       :: rc
    out_esmf_mesh = ESMF_MeshCreate(parametricDim = dim1, spatialDim = spacedim, &
                                    nodeIDs = the_data % NdIDs, nodeCoords = the_data % NdCoords, &
                                    nodeOwners = the_data % NdOwners, elementIDs = the_data % ElIDs, &
                                    elementTypes = the_data % ElTypes, elementConn = the_data % ElConnect, &
                                    rc = rc)

    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

  END SUBROUTINE
  !-----------------------------------------------------------------------
  !- Sub !!!????
  !-----------------------------------------------------------------------

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

    CHARACTER(ESMF_MAXPATHLEN), INTENT(IN), OPTIONAL :: ConfFile  ! config file name
    CHARACTER(LEN = *), INTENT(IN)                  :: OptPrefix ! config file name

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
                                 default = 'moddata.nc',                      &
                                 rc = rc)

    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! Destroy the Config
    CALL ESMF_ConfigDestroy(cf, rc = rc)

  END SUBROUTINE ReadConfig

!================================================================================

END MODULE Bardat_Mod
