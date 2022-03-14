!----------------------------------------------------------------
!               M O D U L E   B A R D A T _ C A P
!----------------------------------------------------------------
!> @file bardat_cap.F90
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

MODULE BarData

  USE MPI
  USE ESMF
  USE NUOPC
  USE NUOPC_Model, &
    model_routine_SS => SetServices, &
    model_label_SetClock => label_SetClock, &
    model_label_Advance => label_Advance, &
    model_label_CheckImport => label_CheckImport, &
    model_label_Finalize => label_Finalize

  USE Bardat_Mod, only: meshdata
  USE Bardat_Mod, only: create_parallel_esmf_mesh_from_meshdata
  !USE Bardat_Mod, only: atm_int,atm_num,atm_den
  USE Bardat_Mod, only: UWND, VWND, PRES
  USE Bardat_Mod, only: ReadConfig

  !read from netcdf file
  USE Bardat_Mod, only: NCDF_InitBarData, NCDF_ReadBarData
  USE Bardat_Mod, only: construct_meshdata_from_netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC SetServices

  TYPE fld_list_type
    CHARACTER(LEN = 64) :: stdname
    CHARACTER(LEN = 64) :: shortname
    CHARACTER(LEN = 64) :: unit
    LOGICAL             :: assoc    ! is the farrayPtr associated with internal data
    REAL(ESMF_KIND_R8), DIMENSION(:), POINTER :: farrayPtr
  END TYPE fld_list_type

  INTEGER, PARAMETER    :: fldsMax = 100
  INTEGER               :: fldsToWav_num = 0
  TYPE(fld_list_type)   :: fldsToWav(fldsMax)
  INTEGER               :: fldsFrATM_num = 0
  TYPE(fld_list_type)   :: fldsFrATM(fldsMax)

  TYPE(meshdata), save  :: mdataInw, mdataOutw
  INTEGER, SAVE         :: iwind_test = 0
  CHARACTER(LEN = 2048) :: info
  INTEGER               :: dbrc     ! temporary debug rc value


  CONTAINS


  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   S E T S E R V I C E S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   NUOPC SetService method is the only public entry point.
  !>
  !> @details
  !>   SetServices registers all of the user-provided subroutines
  !>   in the module with the NUOPC layer.
  !>
  !> @param[in]
  !>   model   An ESMF_GridComp object
  !> @param[out]
  !>     rc  Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE SetServices(model, rc)

    IMPLICIT NONE

    TYPE(ESMF_GridComp)  :: model
    INTEGER, INTENT(OUT) :: rc

    ! Local Variables
    TYPE(ESMF_VM)                 :: vm
    CHARACTER(LEN = *), PARAMETER :: subname = '(BARDATA:SetServices)'

    rc = ESMF_SUCCESS

    ! readd config file
    CALL ReadConfig('bar')

    ! the NUOPC model component will register the generic methods
    CALL NUOPC_CompDerive(model, model_routine_SS, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! set entry point for methods that require specific implementation
    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                                 phaseLabelList = (/"IPDv00p1"/), userRoutine = InitializeP1, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                                 phaseLabelList = (/"IPDv00p2"/), userRoutine = InitializeP2, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! attach specializing method(s)
    !CALL NUOPC_CompSpecialize(model, specLabel = model_label_SetClock, &
    !  specRoutine = SetClock, rc = rc)
    !IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
    !  line = __LINE__, &
    !  file = __FILE__)) &
    !  RETURN  ! bail out
    CALL NUOPC_CompSpecialize(model, specLabel = model_label_Advance, &
                              specRoutine = ModelAdvance, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !CALL NUOPC_CompSpecialize(model, specLabel = model_label_CheckImport, &
    !  specPhaseLabel = "RunPhase1", specRoutine = CheckImport, rc = rc)
    !IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
    !  line = __LINE__, &
    !  file = __FILE__)) &
    !  RETURN  ! bail out

    CALL NUOPC_CompSpecialize(model, specLabel = model_label_Finalize, &
                              specRoutine = BARDATA_model_finalize, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL NCDF_InitBarData()

    WRITE (info, *) subname, ' --- Read BARDATA info from file --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

    WRITE (info, *) subname, ' --- adc SetServices completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE SetServices

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   I N I T I A L I Z E P 1
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Initialization subroutine called by NUOPC.
  !>
  !> @details
  !>   The purpose s to set which version of the
  !>   Initialize Phase Definition (IPD) to use.
  !>
  !> @param[in]
  !>   model         An ESMF_GridComp object
  !> @param[in]
  !>   importState   An ESMF_State object for import fields
  !> @param[in]
  !>   exportState   An ESMF_State object for export fields
  !> @param[in]
  !>   clock         An ESMF_Clock object
  !> @param[out]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE InitializeP1(model, importState, exportState, clock, rc)

    IMPLICIT NONE

    TYPE(ESMF_GridComp)  :: model
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: clock
    INTEGER, INTENT(OUT) :: rc

    ! Local Variables
    INTEGER              :: num, i
    CHARACTER(LEN = *), PARAMETER  :: subname = '(BARDATA:AdvertiseFields)'

    rc = ESMF_SUCCESS

    CALL BARDATA_FieldsSetup()
!

    DO num = 1, fldsToWav_num
      !PRINT *,  "fldsToWav_num  ", fldsToWav_num
      !PRINT *,  "fldsToWav(num)%shortname  ", fldsToWav(num)%shortname
      !PRINT *,  "fldsToWav(num)%stdname  ", fldsToWav(num)%stdname

      WRITE (info, *) subname, "fldsToWav(num)%shortname  ", fldsToWav(num) % shortname
      CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)
    END DO

    CALL BARDATA_AdvertiseFields(importState, fldsToWav_num, fldsToWav, rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

!----------------------------------------------------------------
    DO num = 1, fldsFrATM_num
      !PRINT *,  "fldsFrATM_num  ", fldsFrATM_num
      !PRINT *,  "fldsFrATM(num)%shortname  ", fldsFrATM(num)%shortname
      !PRINT *,  "fldsFrATM(num)%stdname  ", fldsFrATM(num)%stdname
      WRITE (info, *) subname, "fldsFrATM(num)%stdname  ", fldsFrATM(num) % stdname
      CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)

    END DO
!
    CALL BARDATA_AdvertiseFields(exportState, fldsFrATM_num, fldsFrATM, rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    WRITE (info, *) subname, ' --- initialization phase 1 completed --- '
    !PRINT *,      subname,' --- initialization phase 1 completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)
!
  END SUBROUTINE InitializeP1

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   I N I T I A L I Z E P 2
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Initialization subroutine called by NUOPC to realize import and export fields.
  !>
  !> @details
  !> The fields to import and export are stored in the fldsTo* and fldsFr*
  !> arrays, respectively.  Each field entry includes the standard name,
  !> information about whether the field's grid will be provided by the cap,
  !> and optionally a pointer to the field's data array. Currently, all fields
  !> are defined on the same mesh defined by the cap.
  !>
  !> @param[in]
  !>   model         An ESMF_GridComp object
  !> @param[in]
  !>   importState   An ESMF_State object for import fields
  !> @param[in]
  !>   exportState   An ESMF_State object for export fields
  !> @param[in]
  !>   clock         An ESMF_Clock object
  !> @param[out]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE InitializeP2(model, importState, exportState, clock, rc)

    IMPLICIT NONE

    TYPE(ESMF_GridComp)  :: model
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: clock, driverClock
    INTEGER, INTENT(OUT) :: rc

    ! local variables
    TYPE(ESMF_TimeInterval) :: BARDATATimeStep
    TYPE(ESMF_Field)        :: field
    !Saeed added
    TYPE(meshdata)               :: mdataw
    TYPE(ESMF_Mesh)              :: ModelMesh, meshIn, meshOut
    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_Time)              :: startTime
    INTEGER                      :: localPet, petCount
    CHARACTER(LEN = *), PARAMETER   :: subname = '(BARDATA:RealizeFieldsProvidingGrid)'

    rc = ESMF_SUCCESS

    !PRINT *,"BARDATA ..1.............................................. >> "
    !> \details Get current ESMF VM.
    CALL ESMF_VMGetCurrent(vm, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !PRINT *,"BARDATA ..2.............................................. >> "
    ! Get query local pet information for handeling global node information
    CALL ESMF_VMGet(vm, localPet = localPet, petCount = petCount, rc = rc)
    ! CALL ESMF_VMPrint(vm, rc = rc)

    !PRINT *,localPet,"< LOCAL pet, BARDATA ..3.............................................. >> "
    !! Assign VM to mesh data type.
    mdataw % vm = vm

    CALL construct_meshdata_from_netcdf(mdataw)

    CALL create_parallel_esmf_mesh_from_meshdata(mdataw, ModelMesh)
    !

    CALL ESMF_MeshWRITE(ModelMesh, filename = "bardata_mesh.nc", rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    meshIn = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataInw = mdataw
    mdataOutw = mdataw

    !PRINT *,"..................................................... >> "
    !PRINT *,"NumNd", mdataw%NumNd
    !PRINT *,"NumOwnedNd", mdataw%NumOwnedNd
    !PRINT *,"NumEl", mdataw%NumEl
    !PRINT *,"NumND_per_El", mdataw%NumND_per_El
    !PRINT *,"NumOwnedNd mdataOutw", mdataOutw%NumOwnedNd

    CALL BARDATA_RealizeFields(importState, meshIn, mdataw, fldsToWav_num, fldsToWav, "BARDATA import", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
!
    CALL BARDATA_RealizeFields(exportState, meshOut, mdataw, fldsFrATM_num, fldsFrATM, "BARDATA export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !Init BARDATA
!    ! query Component for the driverClock
!    CALL NUOPC_ModelGet(model, driverClock = driverClock, rc = rc)
!    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!      line = __LINE__, &
!      file = __FILE__)) &
!      RETURN  ! bail out

    ! get the start time and current time out of the clock
!    CALL ESMF_ClockGet(driverClock, startTime = startTime, rc = rc)
!    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!      line = __LINE__, &
!      file = __FILE__)) &
!      RETURN  ! bail out

!    CALL NCDF_ReadBarData(startTime)

    WRITE (info, *) subname, ' --- initialization phase 2 completed --- '
    !PRINT *,      subname,' --- initialization phase 2 completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = dbrc)

  END SUBROUTINE InitializeP2

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   B A R D A T A _ A D V E R T I S E F I E L D S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Advertises a set of fields in an ESMF_State object.
  !>
  !> @details
  !>   Advertises a set of fields in an ESMF_State object by calling
  !>   NUOPC_Advertise in a loop.
  !>
  !> @param[inout]
  !>   state         The ESMF_State object in which to advertise the fields
  !> @param[in]
  !>   nfields       The number of fields
  !> @param[inout]
  !>   field_defs    An array of fld_list_type listing the fields to advertise
  !> @param[inout]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE BARDATA_AdvertiseFields(state, nfields, field_defs, rc)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(INOUT)             :: state
    INTEGER, INTENT(IN)                          :: nfields
    TYPE(fld_list_type), INTENT(INOUT)          :: field_defs(:)
    INTEGER, INTENT(INOUT)                      :: rc

    INTEGER                                     :: i
    CHARACTER(LEN = *), PARAMETER  :: subname = '(BARDATA:BARDATA_AdvertiseFields)'

    rc = ESMF_SUCCESS

    DO i = 1, nfields
      !PRINT *, 'Advertise: '//TRIM(field_defs(i)%stdname)//'---'//TRIM(field_defs(i)%shortname)
      CALL ESMF_LogWrite('Advertise: '//TRIM(field_defs(i) % stdname), ESMF_LOGMSG_INFO, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out

      CALL NUOPC_Advertise(state, &
                           standardName = field_defs(i) % stdname, &
                           name = field_defs(i) % shortname, &
                           rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END DO
    !PRINT *,      subname,' --- IN   --- '

  END SUBROUTINE BARDATA_AdvertiseFields

!================================================================================


  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   B A R D A T A _ F I E L D S S E T U P
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>
  !>
  !> @details
  !>
  !>
  !> @param
  !>
  !----------------------------------------------------------------
  SUBROUTINE BARDATA_FieldsSetup

    IMPLICIT NONE

    INTEGER                        :: rc
    CHARACTER(LEN = *), PARAMETER  :: subname = '(BARDATA:BARDATA_FieldsSetup)'

    !--------- import fields to BARDATA  -------------

    !--------- export fields from BARDATA -------------
    CALL fld_list_add(num = fldsFrATM_num, fldlist = fldsFrATM, stdname = "air_pressure_at_sea_level", shortname = "pmsl")
    CALL fld_list_add(num = fldsFrATM_num, fldlist = fldsFrATM, stdname = "inst_zonal_wind_height10m", shortname = "izwh10m")
    CALL fld_list_add(num = fldsFrATM_num, fldlist = fldsFrATM, stdname = "inst_merid_wind_height10m", shortname = "imwh10m")

    WRITE (info, *) subname, ' --- Passed--- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE BARDATA_FieldsSetup

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   F L D _ L I S T _ A D D
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>
  !>
  !> @details
  !>
  !>
  !> @param
  !>
  !----------------------------------------------------------------
  SUBROUTINE fld_list_add(num, fldlist, stdname, data, shortname, unit)

    IMPLICIT NONE

    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    INTEGER, INTENT(INOUT)  :: num
    TYPE(fld_list_type), INTENT(INOUT)  :: fldlist(:)
    CHARACTER(LEN = *), INTENT(IN)     :: stdname
    REAL(ESMF_KIND_R8), DIMENSION(:), OPTIONAL, target :: data
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: shortname
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: unit

    ! local variables
    INTEGER :: rc
    CHARACTER(LEN = *), PARAMETER :: subname = '(BARDATA:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    IF (num > fldsMax) THEN
      CALL ESMF_LogWrite(TRIM(subname)//": ERROR num gt fldsMax "//TRIM(stdname), &
                         ESMF_LOGMSG_ERROR, line = __LINE__, file = __FILE__, rc = rc)
      RETURN
    END IF

    fldlist(num) % stdname = TRIM(stdname)
    IF (PRESENT(shortname)) THEN
      fldlist(num) % shortname = TRIM(shortname)
    ELSE
      fldlist(num) % shortname = TRIM(stdname)
    END IF

    IF (PRESENT(data)) THEN
      fldlist(num) % assoc = .true.
      fldlist(num) % farrayPtr => data
    ELSE
      fldlist(num) % assoc = .FALSE.
    END IF

    IF (PRESENT(unit)) THEN
      fldlist(num) % unit = unit
    END IF

    WRITE (info, *) subname, ' --- Passed--- '

  END SUBROUTINE fld_list_add

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   B A R D A T A _ R E A L I Z E F I E L D S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Adds a set of fields to an ESMF_State object.
  !>
  !> @details
  !> Each field is wrapped in an ESMF_Field object.  Memory is either allocated
  !> by ESMF or an existing BARDATA pointer is referenced.
  !>
  !> @param[inout]
  !>   state         The ESMF_State object to add fields to
  !> @param[in]
  !>   mesh/grid     The ESMF_Grid object on which to define the fields
  !> @param[in]
  !>   nfields       The number of fields
  !> @param[inout]
  !>   field_defs    Array of fld_list_type indicating the fields to add
  !> @param[in]
  !>   tag           Used to output to the log
  !> @param[out]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE BARDATA_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(INOUT)             :: state
    TYPE(ESMF_Mesh), INTENT(IN)                 :: mesh
    TYPE(meshdata)                              :: mdata
    INTEGER, INTENT(IN)                         :: nfields
    TYPE(fld_list_type), INTENT(INOUT)          :: field_defs(:)
    CHARACTER(LEN = *), INTENT(IN)              :: tag
    INTEGER, INTENT(INOUT)                      :: rc

    TYPE(ESMF_Field)                            :: field
    INTEGER                                     :: i
    CHARACTER(LEN = *), PARAMETER  :: subname = '(BARDATA:BARDATA_RealizeFields)'

    rc = ESMF_SUCCESS

    DO i = 1, nfields
      field = ESMF_FieldCreate(name = field_defs(i) % shortname, mesh = mesh, &
                               typekind = ESMF_TYPEKIND_R8, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out

      IF (NUOPC_IsConnected(state, fieldName = field_defs(i) % shortname)) THEN

        CALL NUOPC_Realize(state, field = field, rc = rc)
        IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                               line = __LINE__, &
                               file = __FILE__)) &
          RETURN  ! bail out

        CALL ESMF_LogWrite(subname//tag//" Field "//field_defs(i) % stdname//" is connected.", &
                           ESMF_LOGMSG_INFO, &
                           line = __LINE__, &
                           file = __FILE__, &
                           rc = dbrc)

        !PRINT *,      subname,' --- Connected --- '

      ELSE
        CALL ESMF_LogWrite(subname//tag//" Field "//field_defs(i) % stdname//" is not connected.", &
                           ESMF_LOGMSG_INFO, &
                           line = __LINE__, &
                           file = __FILE__, &
                           rc = dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !IF (associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        CALL ESMF_StateRemove(state, (/field_defs(i) % shortname/), rc = rc)
        IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                               line = __LINE__, &
                               file = __FILE__)) &
          RETURN  ! bail out
        !PRINT *,      subname,' --- Not-Connected --- '
        !PRINT *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '
      END IF
    END DO

    WRITE (info, *) subname, ' --- OUT--- '
    !PRINT *,      subname,' --- OUT --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)

  END SUBROUTINE BARDATA_RealizeFields

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   M O D E L A D V A N C E
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Called by NUOPC to advance the model by a single timestep.
  !    TODO: check! this is not what we want.
  !>
  !> @details
  !>   This subroutine copies field data out of the cap import state and into the
  !>   model internal arrays. Then it calls *_Run to make a NN timesteps.
  !>   Finally, it copies the updated arrays into the cap export state.
  !>
  !> @param[in]
  !>   model   An ESMF_GridComp object
  !> @param[out]
  !>     rc  Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE ModelAdvance(model, rc)

    IMPLICIT NONE

    TYPE(ESMF_GridComp)  :: model
    INTEGER, INTENT(OUT) :: rc

    ! local variables
    TYPE(ESMF_Clock)              :: clock
    TYPE(ESMF_State)              :: importState, exportState
    TYPE(ESMF_Time)               :: currTime
    TYPE(ESMF_TimeInterval)       :: timeStep
    CHARACTER(LEN = *), PARAMETER    :: subname = '(BARDATA:ModelAdvance)'
    !tmp vector
    REAL(ESMF_KIND_R8), POINTER   :: tmp(:)

    !imports

    ! exports
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_uwnd(:)
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_vwnd(:)
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_pres(:)

    TYPE(ESMF_StateItem_Flag)     :: itemType
    TYPE(ESMF_Mesh)               :: mesh
    TYPE(ESMF_Field)              :: lfield
    CHARACTER(LEN = 128)            :: fldname, timeStr
    INTEGER                       :: i1
    ! local variables for Get methods
    INTEGER :: YY, MM, DD, H, M, S

    rc = ESMF_SUCCESS
    ! query the Component for its clock, importState and exportState
    CALL NUOPC_ModelGet(model, modelClock = clock, importState = importState, &
                        exportState = exportState, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

    CALL ESMF_ClockPrint(clock, options = "currTime", &
                         preString = "------>Advancing BARDATA from: ", rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_ClockGet(clock, currTime = currTime, timeStep = timeStep, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimePrint(currTime + timeStep, &
                        preString = "------------------BARDATA-------------> to: ", rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy = YY, mm = MM, dd = DD, h = H, m = M, s = S, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    PRINT *, "BARDATA currTime = ", YY, "/", MM, "/", DD, " ", H, ":", M, ":", S
    WRITE (info, *) "BARDATA currTime = ", YY, "/", MM, "/", DD, " ", H, ":", M, ":", S
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)

    CALL ESMF_TimeGet(currTime, timeStringISOFrac = timeStr, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !-----------------------------------------
    !   IMPORT
    !-----------------------------------------

    !-----------------------------------------
    !   EXPORT
    !-----------------------------------------
    !update uwnd, vwnd, pres from nearset time in bardata netcdf file
    !TODO: update file name!!!!
    CALL NCDF_ReadBarData(currTime)

    !pack and send exported fields
    ALLOCATE (tmp(mdataOutw % NumOwnedNd))

    ! >>>>> PACK and send UWND
    CALL State_GetFldPtr(ST = exportState, fldname = 'izwh10m', fldptr = dataPtr_uwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    iwind_test = iwind_test + 1
    !fill only owned nodes for tmp vector
    DO i1 = 1, mdataOutw % NumOwnedNd, 1
      tmp(i1) = UWND(mdataOutw % owned_to_present_nodes(i1), 1)
      !tmp(i1) = iwind_test  * i1 / 100000.0
      !tmp(i1) = -3.0
    END DO
    !assign to field
    dataPtr_uwnd = tmp
    !----------------------------------------
    ! >>>>> PACK and send VWND
    CALL State_GetFldPtr(ST=exportState, fldname = 'imwh10m', fldptr = dataPtr_vwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !fill only owned nodes for tmp vector
    DO i1 = 1, mdataOutw % NumOwnedNd, 1
      tmp(i1) = VWND(mdataOutw % owned_to_present_nodes(i1), 1)
      !tmp(i1) = 15.0
    END DO
    !assign to field
    dataPtr_vwnd = tmp
    !----------------------------------------
    ! >>>>> PACK and send PRES
    CALL State_GetFldPtr(ST = exportState, fldname = 'pmsl', fldptr = dataPtr_pres, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !fill only owned nodes for tmp vector
    DO i1 = 1, mdataOutw % NumOwnedNd, 1
      tmp(i1) = PRES(mdataOutw % owned_to_present_nodes(i1), 1)

      IF (abs(tmp(i1)) > 1e11) THEN
        STOP '  dataPtr_pmsl > mask1 > in BARDATA ! '
      END IF
      !tmp(i1) = 1e4
    END DO
    !assign to field
    dataPtr_pres = tmp
    !----------------------------------------

    !! TODO:  not a right thing to do. we need to fix the grid object mask <<<<<<
    !where(dataPtr_uwnd > 3e4) dataPtr_uwnd = 0.0
    !where(dataPtr_vwnd > 3e4) dataPtr_vwnd = 0.0

  END SUBROUTINE ModelAdvance

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   S T A T E _ G E T F L D P T R
  !----------------------------------------------------------------
  !  author 
  !>
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !>
  !> @param[in]
  !>   ST         The ESMF_State object
  !> @param[in]
  !>   fldName    The name of the fields
  !> @param[in]
  !>   fldPtr     A pointer to 1D array
  !> @param[out]
  !>   rc         Return code
  !----------------------------------------------------------------
  SUBROUTINE State_GetFldPtr(ST, fldname, fldptr, rc, dump, timeStr)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(IN) :: ST
    CHARACTER(LEN = *), INTENT(IN) :: fldname
    REAL(ESMF_KIND_R8), POINTER, INTENT(IN) :: fldptr(:)
    INTEGER, INTENT(OUT), OPTIONAL :: rc
    LOGICAL :: dump
    CHARACTER(LEN = *), INTENT(INOUT), OPTIONAL :: timeStr

    ! local variables
    TYPE(ESMF_Field) :: lfield
    INTEGER :: lrc
    CHARACTER(LEN = *), PARAMETER :: subname = '(BARDATA:State_GetFldPtr)'

    CALL ESMF_StateGet(ST, itemName = TRIM(fldname), field = lfield, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_FieldGet(lfield, farrayPtr = fldptr, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    IF (PRESENT(rc)) rc = lrc

    IF (dump) THEN
      IF (.NOT. PRESENT(timeStr)) timeStr = "_"
      CALL ESMF_FieldWRITE(lfield, &
                           fileName='moddata_field_'//TRIM(fldname)//TRIM(timeStr)//'.nc', &
                           rc = rc, overwrite = .TRUE.)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END IF

  END SUBROUTINE State_GetFldPtr

!================================================================================

  SUBROUTINE CheckImport(model, rc)

    IMPLICIT NONE

    TYPE(ESMF_GridComp)   :: model
    INTEGER, INTENT(OUT)  :: rc

    ! This is the routine that enforces correct time stamps on import Fields

    ! local variables
    TYPE(ESMF_Clock)        :: driverClock
    TYPE(ESMF_Time)         :: startTime, currTime
    TYPE(ESMF_State)        :: importState
    TYPE(ESMF_Field)        :: field
    LOGICAL                 :: atCorrectTime

    rc = ESMF_SUCCESS
    RETURN

!    ! query Component for the driverClock
!    CALL NUOPC_ModelGet(model, driverClock = driverClock, rc = rc)
!    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!      line = __LINE__, &
!      file = __FILE__)) &
!      RETURN  ! bail out
!
!    ! get the start time and current time out of the clock
!    CALL ESMF_ClockGet(driverClock, startTime = startTime, &
!      currTime = currTime, rc = rc)
!    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!      line = __LINE__, &
!      file = __FILE__)) &
!      RETURN  ! bail out

  END SUBROUTINE CheckImport

!================================================================================

  !-----------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling BARDATA_Final.
  !!
  !! @param gcomp the ESMF_GridComp object
  !! @param rc return code
  SUBROUTINE BARDATA_model_finalize(gcomp, rc)

    IMPLICIT NONE

    ! input arguments
    TYPE(ESMF_GridComp)  :: gcomp
    INTEGER, INTENT(OUT) :: rc

    ! local variables
    TYPE(ESMF_Clock)     :: clock
    TYPE(ESMF_Time)                        :: currTime
    CHARACTER(LEN = *), PARAMETER  :: subname = '(BARDATA:bardata_model_finalize)'

    rc = ESMF_SUCCESS

    WRITE (info, *) subname, ' --- finalize called --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)

    CALL NUOPC_ModelGet(gcomp, modelClock = clock, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_ClockGet(clock, currTime = currTime, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    WRITE (info, *) subname, ' --- finalize completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)

  END SUBROUTINE BARDATA_model_finalize

!================================================================================

END MODULE BarData
