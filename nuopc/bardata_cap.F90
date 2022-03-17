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

  USE Bardat_Mod

  IMPLICIT NONE

  PRIVATE


  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE State_GetFldPtr
    MODULE PROCEDURE State_GetFldPtr
    MODULE PROCEDURE State_GetFldPtr2D
  END INTERFACE State_GetFldPtr

  PUBLIC SetServices

  TYPE FldList_T
    CHARACTER(LEN = 64) :: stdName
    CHARACTER(LEN = 64) :: shortName
    CHARACTER(LEN = 64) :: unit
    LOGICAL             :: assoc    ! is the farrayPtr associated with internal data
    LOGICAL             :: connected !++ GML 2011-11-11 for judge if the pressure is connected or not
    REAL(ESMF_KIND_R8), DIMENSION(:), POINTER :: farrayPtr
  END TYPE FldList_T

  INTEGER, PARAMETER :: maxFlds = 100
  INTEGER            :: fldsToWav_num = 0
  TYPE(FldList_T)    :: fldsToWav(maxFlds)
  INTEGER            :: numFldsFromBAR = 0
  TYPE(FldList_T)    :: fldsFromBAR(maxFlds)

  TYPE(GridData_T), SAVE :: mDataInw, mDataOutw
  INTEGER, SAVE          :: iwind_test = 0
  CHARACTER(LEN = 2048)  :: info
  INTEGER                :: dbrc     ! temporary debug rc value


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
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':SetServices)'

    rc = ESMF_SUCCESS

    ! Read config file
    CALL ReadConfig(MODEL_PRFX)

    ! the NUOPC model component will register the generic methods
    CALL NUOPC_CompDerive(model, model_routine_SS, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! set entry point for methods that require specific implementation
    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,   &
                                 phaseLabelList = (/"IPDv00p1"/), &
                                 userRoutine = InitializeP1, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__,  &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,   &
                                 phaseLabelList = (/"IPDv00p2"/), &
                                 userRoutine = InitializeP2, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__,  &
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
                              specRoutine = ModelFinalize, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL InitBarData()

    WRITE (info, *) subname, ' --- ' // TRIM(ADJUSTL(MODEL_NAME)) // ' SetServices completed --- '
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
    CHARACTER(LEN = *), PARAMETER  :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':InitializeP1)'

    rc = ESMF_SUCCESS

    CALL FieldsSetup()

!   DO num = 1, numFldsFromBAR
!     !PRINT *,  "numFldsFromBAR  ", numFldsFromBAR
!     !PRINT *,  "fldsFromBAR(num)%shortName  ", fldsFromBAR(num)%shortName
!     !PRINT *,  "fldsFromBAR(num)%stdName  ", fldsFromBAR(num)%stdName
!     WRITE (info, *) subname, "fldsFromBAR(num)%stdName  ", fldsFromBAR(num) % stdName
!     CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)

!   END DO

    CALL AdvertiseFields(exportState, numFldsFromBAR, fldsFromBAR, rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    WRITE (info, *) subname, ' --- initialization phase 1 completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = dbrc)

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
  !> and optionally a pointer to the field's data array.
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
    TYPE(ESMF_Field)              :: field
    !Saeed added
    TYPE(GridData_T)              :: mData, mDataY, mDataYY, mDataYYY
    TYPE(ESMF_Grid)               :: modelGrid, modelGridY, modelGridYY, modelGridYYY
    TYPE(ESMF_VM)                 :: vm
    TYPE(ESMF_Time)               :: startTime
    INTEGER                       :: localPet, petCount
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':InitializeP2)'

    rc = ESMF_SUCCESS

    ! Get current ESMF VM.
    CALL ESMF_VMGetCurrent(vm, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! Get query local pet information for handeling global node information
    CALL ESMF_VMGet(vm, localPet = localPet, petCount = petCount, rc = rc)
    ! CALL ESMF_VMPrint(vm, rc = rc)

    ! Assign VM to grid data type.
    mData % vm = vm
    mDataY % vm = vm
    mDataYY % vm = vm
    mDataYYY % vm = vm

    ! lon: lon(NX), lat: lat(NY), lonC: lonC(NX), latC: latC(NYY), lats: latS(NYYY)
    CALL construct_griddata_from_netcdf(NX, NY, mData)
    CALL construct_griddata_from_netcdf(NX, NY, mDataY)
    CALL construct_griddata_from_netcdf(NX, NYY, mDataYY)
    CALL construct_griddata_from_netcdf(NX, NYYY, mDataYYY)

    CALL create_parallel_esmf_grid_from_griddata(lon, lat, mData, modelGrid)
    CALL create_parallel_esmf_grid_from_griddata(lonC, lat, mDataY, modelGridY)
    CALL create_parallel_esmf_grid_from_griddata(lon, latC, mDataYY, modelGridYY)
    CALL create_parallel_esmf_grid_from_griddata(lon, latS, mDataYYY, modelGridYYY)
                                      
    !==== BPGX
    CALL RealizeAField(exportState, modelGridY, numFldsFromBAR, fldsFromBAR, 'bpgx', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    !==== BPGY
    CALL RealizeAField(exportState, modelGridYY, numFldsFromBAR, fldsFromBAR, 'bpgy', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    !==== SigTS
    CALL RealizeAField(exportState, modelGrid, numFldsFromBAR, fldsFromBAR, 'sigts', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    !==== MLD
    CALL RealizeAField(exportState, modelGrid, numFldsFromBAR, fldsFromBAR, 'mld', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    !==== NB
    CALL RealizeAField(exportState, modelGrid, numFldsFromBAR, fldsFromBAR, 'nb', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out
    !==== NM
    CALL RealizeAField(exportState, modelGrid, numFldsFromBAR, fldsFromBAR, 'nm', &
                       TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

!   !==== KE
!   CALL RealizeAField(exportState, modelGrid, numFldsFromBAR, fldsFromBAR, 'ke', &
!                      TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out
!   !==== CDisp
!   CALL RealizeAField(exportState, modelGridYYY, numFldsFromBAR, fldsFromBAR, 'cdisp', &
!                      TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out
!   !==== DispX
!   CALL RealizeAField(exportState, modelGridYYY, numFldsFromBAR, fldsFromBAR, 'dispx', &
!                      TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out
!   !==== DispY
!   CALL RealizeAField(exportState, modelGridYYY, numFldsFromBAR, fldsFromBAR, 'dispy', &
!                      TRIM(ADJUSTL(MODEL_NAME)) // " export", rc)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out



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

!    CALL ReadBarData(startTime)

    WRITE (info, *) subname, ' --- initialization phase 2 completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, &
                       line = __LINE__, &
                       file = __FILE__, &
                       rc = dbrc)

  END SUBROUTINE InitializeP2

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   A D V E R T I S E F I E L D S
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
  !>   nFields       The number of fields
  !> @param[inout]
  !>   fieldDefs     An array of FldList_T listing the fields to advertise
  !> @param[inout]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE AdvertiseFields(state, nFields, fieldDefs, rc)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(INOUT) :: state
    INTEGER, INTENT(IN)             :: nFields
    TYPE(FldList_T), INTENT(INOUT)  :: fieldDefs(:)
    INTEGER, INTENT(INOUT)          :: rc

    INTEGER                         :: i
    CHARACTER(LEN = *), PARAMETER   :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':AdvertiseFields)'

    rc = ESMF_SUCCESS

    DO i = 1, nFields
      CALL ESMF_LogWrite('Advertise: '//TRIM(fieldDefs(i) % stdName), ESMF_LOGMSG_INFO, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out

      CALL NUOPC_Advertise(state, &
                           standardName = fieldDefs(i) % stdName, &
                           name = fieldDefs(i) % shortName, &
                           rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END DO

  END SUBROUTINE AdvertiseFields

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   F I E L D S S E T U P
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
  SUBROUTINE FieldsSetup

    IMPLICIT NONE

    INTEGER                        :: rc
    CHARACTER(LEN = *), PARAMETER  :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':FieldsSetup)'

     !==== BPGX
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("east_west_depth_averaged_baroclinic_pressure_gradient")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="eastwest_depth_averaged_baroclinic_pressure_gradient", &
               canonicalUnits="m s-2", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF
     !==== BPGY
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("north_south_depth_averaged_baroclinic_pressure_gradient")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="northsouth_depth_averaged_baroclinic_pressure_gradient", &
               canonicalUnits="m s-2", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF
     !==== SigTS
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("surface_sigmat_density")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="surface_sigmat_density", &
               canonicalUnits="kg m-3", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF
     !==== MLD
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("mixed_layer_depth_ratio")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="mixed_layer_depth_ratio", &
               canonicalUnits="1", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF
     !==== NB
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("buoyancy_frequency_at_the_seabed")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="buoyancy_frequency_at_the_seabed", &
               canonicalUnits="s-1", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF
     !==== NM
     IF (.NOT. NUOPC_FieldDictionaryHasEntry("depth-averaged_buoyancy_frequency")) THEN
       CALL NUOPC_FieldDictionaryAddEntry( &
               standardName="depth-averaged_buoyancy_frequency", &
               canonicalUnits="s-1", rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
        RETURN  ! bail out
      ENDIF

!    !==== KE
!    IF (.NOT. NUOPC_FieldDictionaryHasEntry("depth_averaged_kinetic_energy")) THEN
!      CALL NUOPC_FieldDictionaryAddEntry( &
!              standardName="depth_averaged_kinetic_energy", &
!              canonicalUnits="m2 s-2", rc=rc)
!      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!       RETURN  ! bail out
!     ENDIF
!    !==== CDisp
!    IF (.NOT. NUOPC_FieldDictionaryHasEntry("momentum_dispersion_quadratic_bottom_friction"))_THEN
!      CALL NUOPC_FieldDictionaryAddEntry( &
!              standardName="momentum_dispersion_quadratic_bottom_friction", &
!              canonicalUnits="1", rc=rc)
!      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!       RETURN  ! bail out
!     ENDIF
!    !==== DispX
!    IF (.NOT. NUOPC_FieldDictionaryHasEntry("depth_averaged_x_momentum_dispersion"))_THEN
!      CALL NUOPC_FieldDictionaryAddEntry( &
!              standardName="depth_averaged_x_momentum_dispersion", &
!              canonicalUnits="m s-2", rc=rc)
!      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!       RETURN  ! bail out
!     ENDIF
!    !==== DispY
!    IF (.NOT. NUOPC_FieldDictionaryHasEntry("depth_averaged_y_momentum_dispersion"))_THEN
!      CALL NUOPC_FieldDictionaryAddEntry( &
!              standardName="depth_averaged_y_momentum_dispersion", &
!              canonicalUnits="m s-2", rc=rc)
!      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!       RETURN  ! bail out
!     ENDIF

    !--------- export fields -------------
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "east_west_depth_averaged_baroclinic_pressure_gradient", &
                      shortName = "bpgx")
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "north_south_depth_averaged_baroclinic_pressure_gradient", &
                      shortName = "bpgy")
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "surface_sigmat_density", &
                      shortName = "sigts")
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "mixed_layer_depth_ratio", &
                      shortName = "mld")
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "buoyancy_frequency_at_the_seabed", &
                      shortName = "nb")
    CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
                      stdName = "depth-averaged_buoyancy_frequency", &
                      shortName = "nm")

!   CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
!                     stdName = "depth_averaged_kinetic_energy", &
!                     shortName = "ke")
!   CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
!                     stdName = "momentum_dispersion_quadratic_bottom_friction", &
!                     shortName = "cdisp")
!   CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
!                     stdName = "depth_averaged_x_momentum_dispersion", &
!                     shortName = "dispx")
!   CALL FieldListAdd(num = numFldsFromBAR, fldList = fldsFromBAR, &
!                     stdName = "depth_averaged_y_momentum_dispersion", &
!                     shortName = "dispy")

    WRITE (info, *) subname, ' --- Passed--- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc = rc)

  END SUBROUTINE FieldsSetup

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   F I E L D L I S T A D D
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
  SUBROUTINE FieldListAdd(num, fldList, stdName, data, shortName, unit)

    IMPLICIT NONE

    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    INTEGER, INTENT(INOUT)                             :: num
    TYPE(FldList_T), INTENT(INOUT)                     :: fldList(:)
    CHARACTER(LEN = *), INTENT(IN)                     :: stdName
    REAL(ESMF_KIND_R8), DIMENSION(:), OPTIONAL, target :: data
    CHARACTER(LEN = *), INTENT(IN)                     :: shortName
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL           :: unit

    ! local variables
    INTEGER :: rc
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':FieldListAdd)'

    ! fill in the new entry

    num = num + 1
    IF (num > maxFlds) THEN
      CALL ESMF_LogWrite(TRIM(subname)//": ERROR num gt maxFlds "//TRIM(stdName), &
                         ESMF_LOGMSG_ERROR, line = __LINE__, file = __FILE__, rc = rc)
      RETURN
    END IF

    fldList(num) % stdName = TRIM(ADJUSTL(stdName))
    fldList(num) % shortName = TRIM(ADJUSTL(shortName))

    IF (PRESENT(data)) THEN
      fldList(num) % assoc = .true.
      fldList(num) % farrayPtr => data
    ELSE
      fldList(num) % assoc = .FALSE.
    END IF

    IF (PRESENT(unit)) THEN
      fldList(num) % unit = unit
    END IF

    WRITE (info, *) subname, ' --- Passed--- '

  END SUBROUTINE FieldListAdd

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   R E A L I Z E F I E L D S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Adds a set of fields to an ESMF_State object.
  !>
  !> @details
  !>   Each field is wrapped in an ESMF_Field object.  Memory is either allocated
  !>   by ESMF or an existing data pointer is referenced.
  !>
  !> @param[inout]
  !>   state         The ESMF_State object to add fields to
  !> @param[in]
  !>   theGrid       The ESMF_Grid object on which to define the fields
  !> @param[in]
  !>   nFields       The number of fields
  !> @param[inout]
  !>   fieldDefs     Array of FldList_T indicating the fields to add
  !> @param[in]
  !>   tag           Used to output to the log
  !> @param[out]
  !>   rc            Return code
  !>
  !----------------------------------------------------------------
  SUBROUTINE RealizeFields(state, theGrid, nFields, fieldDefs, tag, rc)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(INOUT) :: state
    TYPE(ESMF_Grid), INTENT(IN)     :: theGrid
    INTEGER, INTENT(IN)             :: nFields
    TYPE(FldList_T), INTENT(INOUT)  :: fieldDefs(:)
    CHARACTER(LEN = *), INTENT(IN)  :: tag
    INTEGER, INTENT(INOUT)          :: rc

    TYPE(ESMF_Field)                :: field
    INTEGER                         :: i
    CHARACTER(LEN = *), PARAMETER   :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':RealizeFields)'

    rc = ESMF_SUCCESS

    DO i = 1, nFields
      field = ESMF_FieldCreate(name = fieldDefs(i) % shortName, grid = theGrid, &
                               typekind = ESMF_TYPEKIND_R8, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out

      IF (NUOPC_IsConnected(state, fieldName = fieldDefs(i) % shortName)) THEN

        CALL NUOPC_Realize(state, field = field, rc = rc)
        IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                               line = __LINE__, &
                               file = __FILE__)) &
          RETURN  ! bail out

        CALL ESMF_LogWrite(subname//tag//" Field "//fieldDefs(i) % stdName//" is connected.", &
                           ESMF_LOGMSG_INFO, &
                           line = __LINE__, &
                           file = __FILE__, &
                           rc = dbrc)
      ELSE
        CALL ESMF_LogWrite(subname//tag//" Field "//fieldDefs(i) % stdName//" is not connected.", &
                           ESMF_LOGMSG_INFO, &
                           line = __LINE__, &
                           file = __FILE__, &
                           rc = dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !IF (associated(fieldDefs(i)%farrayPtr) ) fieldDefs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        CALL ESMF_StateRemove(state, (/fieldDefs(i) % shortName/), rc = rc)
        IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                               line = __LINE__, &
                               file = __FILE__)) &
          RETURN  ! bail out
      END IF
    END DO

    WRITE (info, *) subname, ' --- OUT--- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)

  END SUBROUTINE RealizeFields

!================================================================================

  SUBROUTINE RealizeAField(state, theGrid, nFields, fieldDefs, sName, tag, rc)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(INOUT) :: state
    TYPE(ESMF_Grid), INTENT(IN)     :: theGrid
    INTEGER, INTENT(IN)             :: nFields
    TYPE(FldList_T), INTENT(INOUT)  :: fieldDefs(:)
    CHARACTER(LEN = *), INTENT(IN)  :: sName
    CHARACTER(LEN = *), INTENT(IN)  :: tag
    INTEGER, INTENT(INOUT)          :: rc

    TYPE(ESMF_Field)                :: field
    INTEGER                         :: i, iFLD
    CHARACTER(LEN = *), PARAMETER   :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':RealizeAField)'

    rc = ESMF_SUCCESS

    iFLD = -1
    DO i = 1, nFields
      IF (TRIM(ADJUSTL(ToLowerCase(fieldDefs(i) % shortName))) == &
          TRIM(ADJUSTL(ToLowerCase(sName)))) THEN
        iFLD = i
        EXIT
      END IF
    END DO

    IF (iFLD <= 0) THEN
      WRITE (info, *) subname, ' --- field name ' // TRIM(ADJUSTL((sName))) // ' not in the field list --- '
      CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)
      RETURN
    END IF

    field = ESMF_FieldCreate(name = fieldDefs(iFLD) % shortName, grid = theGrid, &
                             typekind = ESMF_TYPEKIND_R8, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    IF (NUOPC_IsConnected(state, fieldName = fieldDefs(iFLD) % shortName)) THEN

      CALL NUOPC_Realize(state, field = field, rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out

      CALL ESMF_LogWrite(subname//tag//" Field "//fieldDefs(iFLD) % stdName//" is connected.", &
                         ESMF_LOGMSG_INFO, &
                         line = __LINE__, &
                         file = __FILE__, &
                         rc = dbrc)
    ELSE
      CALL ESMF_LogWrite(subname//tag//" Field "//fieldDefs(iFLD) % stdName//" is not connected.", &
                         ESMF_LOGMSG_INFO, &
                         line = __LINE__, &
                         file = __FILE__, &
                         rc = dbrc)
      ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
      !IF (associated(fieldDefs(iFLD)%farrayPtr) ) fieldDefs(iFLD)%farrayPtr = 0.0
      ! remove a not connected Field from State
      CALL ESMF_StateRemove(state, (/fieldDefs(iFLD) % shortName/), rc = rc)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END IF

    WRITE (info, *) subname, ' --- OUT--- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)

  END SUBROUTINE RealizeAField


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
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':ModelAdvance)'
    !tmp vector
!    REAL(ESMF_KIND_R8), POINTER   :: tmp(:)

    !imports

    ! exports
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_bpgx(:, :),  dataPtr_bpgy(:, :)
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_sigts(:, :), dataPtr_mld(:, :)
    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_nb(:, :),    dataPtr_nm(:, :)

    REAL(ESMF_KIND_R8), POINTER   :: dataPtr_ke(:, :), dataPtr_cdisp(:, :), &
                                     dataPtr_dispx(:, :), dataPtr_dispy(:, :)

    TYPE(ESMF_StateItem_Flag)     :: itemType
    TYPE(ESMF_Grid)               :: theGrid
    TYPE(ESMF_Field)              :: lfield
    CHARACTER(LEN = 128)          :: fldName, timeStr

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
                         preString = "------>Advancing " // TRIM(ADJUSTL(MODEL_NAME)) // " from: ", &
                         rc = rc)
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
                        preString = "------------------" //      &
                                    TRIM(ADJUSTL(MODEL_NAME)) // &
                                    "-------------> to: ", rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy = YY, mm = MM, dd = DD, h = H, m = M, s = S, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    PRINT *, TRIM(ADJUSTL(MODEL_NAME)) // " currTime = ", &
             YY, "/", MM, "/", DD, " ", H, ":", M, ":", S
    WRITE (info, *) TRIM(ADJUSTL(MODEL_NAME)) // " currTime = ", &
                    YY, "/", MM, "/", DD, " ", H, ":", M, ":", S
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line = __LINE__, file = __FILE__, rc = rc)

    CALL ESMF_TimeGet(currTime, timeStringISOFrac = timeStr, rc = rc)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    !-----------------------------------------
    !   EXPORT
    !-----------------------------------------
    ! Update the fields from the nearest time in the data NetCDF file
    CALL ReadBarData(currTime)

    !===== Pack and send BPGX
    ALLOCATE(dataPtr_bpgx(NX, NY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'bpgx', fldPtr = dataPtr_bpgx, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_bpgx = BPGX(:, :, 1)

    !===== Pack and send BPGY
    ALLOCATE(dataPtr_bpgy(NX, NYY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'bpgy', fldPtr = dataPtr_bpgy, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_bpgy = BPGY(:, :, 1)

    !===== Pack and send SigTS
    ALLOCATE(dataPtr_sigts(NX, NY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'sigts', fldPtr = dataPtr_sigts, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_sigts = SigTS(:, :, 1)

    !===== Pack and send MLD
    ALLOCATE(dataPtr_mld(NX, NY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'mld', fldPtr = dataPtr_mld, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_mld = MLD(:, :, 1)

    !===== Pack and send NB
    ALLOCATE(dataPtr_nb(NX, NY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'nb', fldPtr = dataPtr_nb, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_nb = NB(:, :, 1)

    !===== Pack and send NM
    ALLOCATE(dataPtr_nm(NX, NY))

    CALL State_GetFldPtr(ST=exportState, fldName = 'nm', fldPtr = dataPtr_nm, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_nm = NM(:, :, 1)

!   !===== Pack and send KE
!   ALLOCATE(dataPtr_ke(NX, NY))

!   CALL State_GetFldPtr(ST=exportState, fldName = 'ke', fldPtr = dataPtr_ke, &
!                        rc = rc, dump = .FALSE., timeStr = timeStr)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out

!   dataPtr_ke = KE(:, :, 1)

!   !===== Pack and send CDisp
!   ALLOCATE(dataPtr_cdisp(NX, NYYY))

!   CALL State_GetFldPtr(ST=exportState, fldName = 'cdisp', fldPtr = dataPtr_cdisp, &
!                        rc = rc, dump = .FALSE., timeStr = timeStr)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out

!   dataPtr_cdisp = CDisp(:, :, 1)

!   !===== Pack and send DispX
!   ALLOCATE(dataPtr_dispx(NX, NYYY))

!   CALL State_GetFldPtr(ST=exportState, fldName = 'dispx', fldPtr = dataPtr_dispx, &
!                        rc = rc, dump = .FALSE., timeStr = timeStr)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out

!   dataPtr_dispx = DispX(:, :, 1)

!   !===== Pack and send DispY
!   ALLOCATE(dataPtr_dispy(NX, NYYY))

!   CALL State_GetFldPtr(ST=exportState, fldName = 'dispy', fldPtr = dataPtr_dispy, &
!                        rc = rc, dump = .FALSE., timeStr = timeStr)
!   IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
!                          line = __LINE__, &
!                          file = __FILE__)) &
!     RETURN  ! bail out

!   dataPtr_dispy = DispY(:, :, 1)


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
  SUBROUTINE State_GetFldPtr(ST, fldName, fldPtr, rc, dump, timeStr)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(IN)                :: ST
    CHARACTER(LEN = *), INTENT(IN)              :: fldName
    REAL(ESMF_KIND_R8), POINTER, INTENT(IN)     :: fldPtr(:)
    INTEGER, INTENT(OUT), OPTIONAL              :: rc
    LOGICAL, INTENT(IN), OPTIONAL               :: dump
    CHARACTER(LEN = *), INTENT(INOUT), OPTIONAL :: timeStr

    ! local variables
    TYPE(ESMF_Field)              :: lfield
    INTEGER                       :: lrc
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':State_GetFldPtr)'

    ! ESMF ref: 21.7.10- Get an item from a State by item name.
    CALL ESMF_StateGet(ST, itemName = TRIM(fldName), field = lfield, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    CALL ESMF_FieldGet(lfield, farrayPtr = fldPtr, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    IF (PRESENT(rc)) rc = lrc

    IF (dump) THEN
      IF (.NOT. PRESENT(timeStr)) timeStr = "_"
      CALL ESMF_FieldWRITE(lfield, &
                           fileName = TRIM(ADJUSTL(ToLowerCase(MODEL_NAME))) // '_field_' // &
                                      TRIM(ADJUSTL(fldName)) // TRIM(ADJUSTL(timeStr)) // '.nc', &
                           rc = rc, overwrite = .TRUE.)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END IF

  END SUBROUTINE State_GetFldPtr

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   S T A T E _ G E T F L D P T R 2 D
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
  !>   fldPtr     A pointer to 2D array
  !> @param[out]
  !>   rc         Return code
  !----------------------------------------------------------------
  SUBROUTINE State_GetFldPtr2D(ST, fldName, fldPtr, rc, dump, timeStr)

    IMPLICIT NONE

    TYPE(ESMF_State), INTENT(IN)                :: ST
    CHARACTER(LEN = *), INTENT(IN)              :: fldName
    REAL(ESMF_KIND_R8), POINTER, INTENT(IN)     :: fldPtr(:, :)
    INTEGER, INTENT(OUT), OPTIONAL              :: rc
    LOGICAL, INTENT(IN), OPTIONAL               :: dump
    CHARACTER(LEN = *), INTENT(INOUT), OPTIONAL :: timeStr

    !===== Local variables
    TYPE(ESMF_Field)              :: lfield
    INTEGER                       :: lrc
    CHARACTER(LEN = *), PARAMETER :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':State_GetFldPtr2D)'

    ! ESMF ref: 21.7.10- Get an item from a State by item name.
    CALL ESMF_StateGet(ST, itemName = TRIM(fldName), field = lfield, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    ! ESMF ref: 26.6.46 - Get a DE-local Fortran array pointer from a Field
    ! Description: get a fortran pointer to DE-local memory allocation within
    ! field, DE-local bounds can be queried at the same time. see section
    ! 26.3.2.
    CALL ESMF_FieldGet(lfield, farrayPtr = fldPtr, localDe = 0, rc = lrc)
    IF (ESMF_LogFoundError(rcToCheck = lrc, msg = ESMF_LOGERR_PASSTHRU, &
                           line = __LINE__, &
                           file = __FILE__)) &
      RETURN  ! bail out

    IF (PRESENT(rc)) rc = lrc

    IF (dump) THEN
      IF (.NOT. PRESENT(timeStr)) timeStr = "_"
      CALL ESMF_FieldWRITE(lfield, &
                           fileName = TRIM(ADJUSTL(ToLowerCase(MODEL_NAME))) // '_field_' // &
                                      TRIM(ADJUSTL(fldName)) // TRIM(ADJUSTL(timeStr)) // '.nc', &
                           rc = rc, overwrite = .TRUE.)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
                             line = __LINE__, &
                             file = __FILE__)) &
        RETURN  ! bail out
    END IF

  END SUBROUTINE State_GetFldPtr2D

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
  !! this simply by calling ModelFinalize.
  !!
  !! @param gcomp the ESMF_GridComp object
  !! @param rc return code
  SUBROUTINE ModelFinalize(gcomp, rc)

    IMPLICIT NONE

    ! input arguments
    TYPE(ESMF_GridComp)  :: gcomp
    INTEGER, INTENT(OUT) :: rc

    ! local variables
    TYPE(ESMF_Clock)     :: clock
    TYPE(ESMF_Time)      :: currTime
    CHARACTER(LEN = *), PARAMETER  :: subname = '(' // TRIM(ADJUSTL(MODEL_NAME)) // ':ModelFinalize)'

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

  END SUBROUTINE ModelFinalize

!================================================================================

END MODULE BarData
