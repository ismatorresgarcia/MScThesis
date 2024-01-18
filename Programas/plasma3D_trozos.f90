!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     plasma class for OfiKinRad input
!     plasma variables (electron density and temperature, populations, rates...)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Eduardo Oliva Gonzalo
!     10-02-2015
!     Last revision 23-09-2015

  module m_plasma3D

    use m_constants
    use m_writeOut
    use fgsl
    
    implicit none

    private

    type, public :: t_plasma3D

       private
       
       class(t_writeOut), pointer :: pWriteOut

       integer :: ntsteps, nlevel, rtsteps, Nz, Nx, Ny, ulevel, llevel
       character(len=MAX_FILE_LENGTH) :: plasmaFile, levelFile, radRateFile, depopUFile, depopLFile
       character(len=MAX_FILE_LENGTH) :: popUFile, popLFile, nueiFile

       real(fgsl_double), dimension(:), pointer :: timePlasma, Te, ne, Zbar, timeRates, nuei
       real(fgsl_double), dimension(:), pointer :: depopU, depopL, NU, NL, totDPU, totDPL, totPU, totPL
       real(fgsl_double), dimension(:), pointer :: totRPU, totRPL
       real(fgsl_double), dimension(:,:), pointer :: levels, radRates, popU, popL
       real(fgsl_double), dimension(:,:,:), pointer :: ne3D, Te3D, Zbar3D, NU3D, NL3D,DPU3D,DPL3D,NZL3D
       real(fgsl_double), dimension(:,:,:), pointer :: PU3D,PL3D, nuei3D, neprof3D
       real(fgsl_double), dimension(:), pointer :: spacex, spacey, spacez                                           
       real(fgsl_double) :: nneutr, tinit, tend, dt0, Lz, Lx, Ly, dx, dy, dz, Aul, DEnergy, zul
       real(fgsl_double) :: totRDPU, totRDPL, plength, pcenter
       integer :: gU, gL, nhgauss

       real(fgsl_double) :: zIR, sigmazIR, IROmega
       logical :: IR_group_vel_flag, pFlag
       
       contains
         
         procedure, public :: setPlasma => allocateVariables
         procedure, public :: initPlasma => initiatePlasma
         procedure, public :: cleanPlasma => deallocateVariables
         procedure, public :: readPlasma => readInputFileOKR
         procedure, public :: readInput => readInputNml
         procedure, public :: getTPlasma => getPointerTPlasma
         procedure, public :: getTe => getPointerTe
         procedure, public :: getNe => getPointerNe
         procedure, public :: getZ => getPointerZ
         procedure, public :: getNU => getPointerNU
         procedure, public :: getNL => getPointerNL
         procedure, public :: getPU => getPointerPU
         procedure, public :: getPL => getPointerPL
         procedure, public :: getDPU => getPointerDPU
         procedure, public :: getDPL => getPointerDPL
         procedure, public :: getNuei => getPointerNuei
         procedure, public :: getAul => getValueAul
         procedure, public :: getZul => getValueZul
         procedure, public :: getDEnergy => getValueDEnergy
         procedure, public :: getDegeneration => getValuegX
         procedure, public :: getSpace => getPointerSpace
         procedure, public :: getLength => getDimLength
         procedure, public :: updatePlasma => fillPlasma
         procedure, public :: outPlasma => outputPlasma
         procedure, public :: initRates => rates
         
      end type t_plasma3D
      
    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     memory allocation
!     does not check if allocated
!     does not return error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine allocateVariables(this)

        implicit none

        class(t_plasma3D), intent(inout) :: this

        integer :: istat

        allocate(this%timePlasma(1:this%ntsteps),STAT=istat)
        this%timePlasma = 0.0

        allocate(this%Te(1:this%ntsteps),STAT=istat)
        this%Te = 0.0
        allocate(this%ne(1:this%ntsteps),STAT=istat)
        this%ne = 0.0
        allocate(this%Zbar(1:this%ntsteps),STAT=istat)
        this%Zbar = 0.0
        allocate(this%nuei(1:this%ntsteps),STAT=istat)
        this%nuei = 0.0
        
        allocate(this%levels(1:this%ntsteps,1:this%nlevel))
        this%levels = 0.0

        allocate(this%radRates(1:this%nlevel,1:this%nlevel))
        this%radRates = 0.0

        allocate(this%timeRates(1:this%rtsteps))
        this%timeRates = 0.0
        allocate(this%depopU(1:this%rtsteps))
        this%depopU = 0.0
        allocate(this%depopL(1:this%rtsteps))
        this%depopL = 0.0

        allocate(this%popU(1:this%rtsteps,1:this%nlevel))
        this%popU = 0.0
        allocate(this%popL(1:this%rtsteps,1:this%nlevel))
        this%popL = 0.0

        allocate(this%totPU(1:this%ntsteps))
        this%totPU = 0.0
        allocate(this%totPL(1:this%ntsteps))
        this%totPU = 0.0
                
        allocate(this%totDPU(1:this%ntsteps),STAT=istat)
        this%totDPU = 0.0
        allocate(this%totDPL(1:this%ntsteps),STAT=istat)
        this%totDPL = 0.0

        allocate(this%totRPU(1:this%ntsteps),STAT=istat)
        this%totRPU = 0.0
        allocate(this%totRPL(1:this%ntsteps),STAT=istat)
        this%totRPL = 0.0

        allocate(this%NZL3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%NZL3D = 0.0
        allocate(this%ne3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%ne3D = 0.0
        allocate(this%Te3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%Te3D = 0.0
        allocate(this%Zbar3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%Zbar3D = 0.0
        allocate(this%DPU3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%DPU3D = 0.0
        allocate(this%DPL3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%DPL3D = 0.0
        allocate(this%PU3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%PU3D = 0.0
        allocate(this%PL3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%PL3D = 0.0
        allocate(this%nuei3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%nuei3D = 0.0
        allocate(this%neprof3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%NZL3D = 0.0

        allocate(this%NU3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%NU3D = 0.0
        allocate(this%NL3D(1:this%Nx,1:this%Ny,1:this%Nz),STAT=istat)
        this%NL3D = 0.0
        
        allocate(this%spacex(1:this%Nx),STAT=istat)
        allocate(this%spacey(1:this%Ny),STAT=istat)
        allocate(this%spacez(1:this%Nz),STAT=istat)

      end subroutine allocateVariables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     deallocation
!     does not check if allocated
!     does not return error message  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine deallocateVariables(this)

        implicit none

        class(t_plasma3D), intent(inout) :: this

        integer :: istat

        deallocate(this%timePlasma,STAT=istat)
        deallocate(this%Te,STAT=istat)
        deallocate(this%ne,STAT=istat)
        deallocate(this%Zbar,STAT=istat)
        deallocate(this%nuei,STAT=istat)
        
        deallocate(this%levels,STAT=istat)
        deallocate(this%radRates,STAT=istat)
        deallocate(this%timeRates,STAT=istat)
        deallocate(this%depopU,STAT=istat)
        deallocate(this%depopL,STAT=istat)
        deallocate(this%popU,STAT=istat)
        deallocate(this%popL,STAT=istat)

        deallocate(this%totPU,STAT=istat)
        deallocate(this%totPL,STAT=istat)

        deallocate(this%totDPU,STAT=istat)
        deallocate(this%totDPL,STAT=istat)

        deallocate(this%totRPU,STAT=istat)
        deallocate(this%totRPL,STAT=istat)

        deallocate(this%NZL3D,STAT=istat)        
        deallocate(this%ne3D,STAT=istat)
        deallocate(this%Te3D,STAT=istat)
        deallocate(this%Zbar3D,STAT=istat)
        deallocate(this%DPU3D,STAT=istat)
        deallocate(this%DPL3D,STAT=istat)
        deallocate(this%PU3D,STAT=istat)
        deallocate(this%PL3D,STAT=istat)
        deallocate(this%NU3D,STAT=istat)
        deallocate(this%NL3D,STAT=istat)
        deallocate(this%nuei3D,STAT=istat)
        deallocate(this%neprof3D,STAT=istat)
        
        deallocate(this%spacex,STAT=istat)
        deallocate(this%spacey,STAT=istat)
        deallocate(this%spacez,STAT=istat)
        
      end subroutine deallocateVariables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to timePlasma
!     pointOut finishes pointing to timePlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerTPlasma(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:), pointer, intent(inout) :: pointOut

        pointOut => this%timePlasma

        end subroutine getPointerTPlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!        this%timeP => this%pPlasma%getFTPlasma()
!     returns a pointer  pointing to timePlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function fPointerTPlasma(this)
        
        implicit none
        
        real(fgsl_double), dimension(:), pointer :: fPointerTPlasma
        class(t_plasma3D), intent(in) :: this

        fPointerTPlasma => this%timePlasma

      end function fPointerTPlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to Te
!     pointOut finishes pointing to Te 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerTe(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%Te3D

        end subroutine getPointerTe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to ne
!     pointOut finishes pointing to ne     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerNe(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%ne3D

        end subroutine getPointerNe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to Zbar
!     pointOut finishes pointing to Zbar     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerZ(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%Zbar3D

      end subroutine getPointerZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to NU3D
!     pointOut finishes pointing to NU3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerNU(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%NU3D

        end subroutine getPointerNU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to NL3D
!     pointOut finishes pointing to NL3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerNL(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%NL3D

      end subroutine getPointerNL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to PU3D
!     pointOut finishes pointing to PU3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerPU(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%PU3D

      end subroutine getPointerPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to PL3D
!     pointOut finishes pointing to PL3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerPL(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%PL3D

      end subroutine getPointerPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to DPU3D
!     pointOut finishes pointing to DPU3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerDPU(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%DPU3D

      end subroutine getPointerDPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to DPL3D
!     pointOut finishes pointing to DPL3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerDPL(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%DPL3D

      end subroutine getPointerDPL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to nuei3D
!     pointOut finishes pointing to nuei3D     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerNuei(this,pointOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:,:,:), pointer, intent(inout) :: pointOut

        pointOut => this%Nuei3D

      end subroutine getPointerNuei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to one of the space vectors
!     pointOut finishes pointing to spacei where i = x,y,z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getPointerSpace(this,pointOut,direction)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), dimension(:), pointer, intent(inout) :: pointOut
        integer, intent(in) :: direction

        if(direction .eq. 1) then
           pointOut => this%spacex
        else if(direction .eq. 2) then
           pointOut => this%spacey
        else if(direction .eq. 3) then
           pointOut => this%spacez
        endif

      end subroutine getPointerSpace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     get a pointer to one of the space vectors
!     pointOut finishes pointing to spacei where i = x,y,z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getDimLength(this,valOut,direction)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), intent(inout) :: valOut
        integer, intent(in) :: direction

        if(direction .eq. 1) then
           valOut = this%Lx
        else if(direction .eq. 2) then
           valOut = this%Ly
        else if(direction .eq. 3) then
           valOut = this%Lz
        endif

      end subroutine getDimLength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getValueAul(this,valOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), intent(out) :: valOut

        valOut = this%Aul

      end subroutine getValueAul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getValueZul(this,valOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), intent(out) :: valOut

        valOut = this%zul

      end subroutine getValueZul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getValueDEnergy(this,valOut)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        real(fgsl_double), intent(out) :: valOut

        valOut = this%DEnergy

      end subroutine getValueDEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getValuegX(this,valOutU,valOutL)
        
        implicit none
        
        class(t_plasma3D), intent(in) :: this
        integer, intent(out) :: valOutU,valOutL

        valOutU = this%gU
        valOutL = this%gL

      end subroutine getValuegX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       reads atomic input files (OfiKinRad/postprocessed)
!       timePlasma, Te, ne, Zbar
!       levels
!       radRates
!       timeRates
!       depopU, depopL, popU, popL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readInputFileOKR(this)

        implicit none

        class(t_plasma3D), intent(inout) :: this

        integer :: funit = 20, k,l
        !character(len=MAX_DUMP_CHAR) :: dumpChar
        real(fgsl_double) :: dumpReal

        open(Unit=funit, File=this%plasmaFile, Status='old')
        do k=1,this%ntsteps
           read(funit,*) this%timePlasma(k), this%Te(k), this%ne(k), this%Zbar(k)
        end do
        close(funit)

        open(Unit=funit, File=this%nueiFile, Status='old')
        do k=1,this%ntsteps
           read(funit,*) this%nuei(k)
        end do
        close(funit)

        open(Unit=funit, File=this%levelFile, Status='old')
        do k=1,this%ntsteps
           ! first column in auxLev is time -> dumped
           read(funit,*) dumpReal, (this%levels(k,l),l=1,this%nlevel)
        end do
        close(funit)

        open(Unit=funit, File=this%radRateFile, Status='old')
        ! first file and first column are indexes
        read(funit,*) dumpReal
        do k=1,this%nlevel
           read(funit,*) dumpReal, (this%radRates(k,l),l=1,this%nlevel)
        end do
        close(funit)
        ! Einstein coefficient of the transition
        this%Aul = this%radRates(this%ulevel,this%llevel)
        ! Dipole matrix element from Einstein A coefficient
        this%zul = sqrt(3*PI*EPS0*HBAR**4*C_SPEED**3*this%Aul/(this%DEnergy*EV2J)**3)
        !write(*,*)"radiative data: ",this%Aul," ",this%zul
        
        open(Unit=funit, File=this%depopUFile, Status='old')
        do k=1,this%rtsteps
           read(funit,*) this%timeRates(k), this%depopU(k)
        end do
        close(funit)

        open(Unit=funit, File=this%depopLFile, Status='old')
        do k=1,this%rtsteps
           !we only read time once, above
           read(funit,*) dumpReal, this%depopL(k)
        end do
        close(funit)

        open(Unit=funit, File=this%popUFile, Status='old')
        do k=1,this%rtsteps
           !time is dumped
           read(funit,*) dumpReal, (this%popU(k,l),l=1,this%nlevel)
        end do
        close(funit)

        open(Unit=funit, File=this%popLFile, Status='old')
        do k=1,this%rtsteps
           !time is dumped
           read(funit,*) dumpReal, (this%popL(k,l),l=1,this%nlevel)
        end do
        close(funit)

      end subroutine readInputFileOKR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     read atomic namelist from input file
!     the file is opened elsewhere (main file)
!     unit is passed as an argument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      subroutine readInputNml(this, unit)
        
        implicit none

        class(t_plasma3D), intent(inout) :: this
        integer,         intent(in)    :: unit

        integer :: time_steps, level_number, rate_time_steps
        real(fgsl_double) :: transition_energy, neutral_density
        character(len=MAX_FILE_LENGTH) :: plasma_file, level_file, rad_file, nuei_file
        character(len=MAX_FILE_LENGTH) :: pop_U_file, pop_L_file, depop_U_file, depop_L_file
        integer :: upper_level, lower_level, upper_level_degeneration, lower_level_degeneration
      
        namelist / atomic / time_steps, level_number, plasma_file, nuei_file, neutral_density, &
             level_file, rad_file, rate_time_steps, depop_U_file, depop_L_file, pop_U_file, &
             pop_L_file, upper_level, lower_level, transition_energy, &
             upper_level_degeneration, lower_level_degeneration

        real(fgsl_double) :: init_time, end_time, init_dt, domain_length, domain_width
        integer :: cell_number_x, cell_number_z, cell_number_y
        integer :: ghost_cell_x, ghost_cell_y, ghost_cell_z
        logical :: adiabatic_approx, ASE_flag, ASE_back_flag

        namelist / common / init_time, end_time, init_dt, domain_length&
             &, cell_number_z, domain_width, cell_number_x, cell_number_y &
             &, ghost_cell_x, ghost_cell_y, ghost_cell_z &
             &, adiabatic_approx, ASE_flag, ASE_back_flag

        logical :: plasma_flag
        real(fgsl_double) :: plasma_length, plasma_center
        integer :: hgauss_index
        
        namelist / plasma / plasma_flag, plasma_length, plasma_center, hgauss_index
        
        logical :: IR_group_vel_flag
        real(fgsl_double) :: tcentIR, tauIR, IRlambda

        namelist / IRlaser / IR_group_vel_flag, tcentIR, tauIR, IRlambda

        read(unit, nml=atomic)

        this%ntsteps = time_steps
        this%nlevel  = level_number
        this%plasmaFile = plasma_file
        this%nneutr = neutral_density
        this%levelFile = level_file
        this%radRateFile = rad_file
        this%rtsteps = rate_time_steps
        this%depopUFile = depop_U_file
        this%depopLFile = depop_L_file
        this%popUFile = pop_U_file
        this%popLFile = pop_L_file
        this%ulevel = upper_level
        this%llevel = lower_level
        this%DEnergy = transition_energy
        this%nueiFile = nuei_file
        this%gU = upper_level_degeneration
        this%gL = lower_level_degeneration
        
        read(unit, nml=common)

        this%tinit = init_time
        this%tend = end_time
        this%dt0 = init_dt

        this%Lz = domain_length
        this%Lx = domain_width
        this%Ly = domain_width
        
        this%Nz = cell_number_z
        this%Nx = cell_number_x
        this%Ny = cell_number_y

        read(unit, nml=plasma)

        this%pFlag = plasma_flag
        this%plength = plasma_length
        this%pcenter = plasma_center
        this%nhgauss = hgauss_index
        
        read(unit, nml=IRlaser)

        this%IR_group_vel_flag = IR_group_vel_flag
        this%zIR = -tcentIR*C_SPEED
        this%sigmazIR = C_SPEED*tauIR/(2*sqrt(2*log(2.0)))
        this%IROmega = 2*PI*C_SPEED/IRlambda
        
      end subroutine readInputNml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fillPlasma(this,ctime,dtstep)

        implicit none

        class(t_plasma3D), intent(inout) :: this
        real(fgsl_double), intent(in) :: ctime, dtstep

        type(fgsl_interp_accel) :: accNe, accTe, accZ, accDPU, accDPL, accPU, accPL, accRPU, accRPL,accNuei
        type(fgsl_spline) :: splineNe, splineTe, splineZ, splineDPU, splineDPL, splinePU, splinePL
        type(fgsl_spline) :: splineRPU, splineRPL, splineNuei
        integer(fgsl_int) :: status
        integer(fgsl_size_t) :: ncells
        integer :: ub1,lb1,m,k,l
        real(fgsl_double) :: itime
        real(fgsl_double) :: radius,r0,rc,rv
        real(fgsl_double) :: sig_r, sig_z, sig_r0, sig_z0, z0_fac, zshift, cen_flag, z0_flag
        real(fgsl_double) :: expion_z0, rlion_z0
        real(fgsl_double), dimension(:), pointer :: expion_cen, rlion_cen !PARAMETROS PARA QUE SEA UN VECTOR LA SIGMOIDE
        real(fgsl_double) :: cen_fac
        real(fgsl_double) :: r_L, r_u, sigma_L, sigma_u, A_ru             !PARÁMETROS PARA LA FUNCION EXPONENCIAL A TROZOS        


        real(fgsl_double) :: sigmaPz

        real(fgsl_double) :: avgNe, avgOmega, avgVgIR
        integer(fgsl_int) :: counter

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!PARTE NUEVA PARA GENERAR LAS SIGMOIDES!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(fgsl_double) :: ki, ri_max, ri_min, z0i, Li !LOS PARÁMETROS DE LA SIGMOIDE
        real(fgsl_double) :: kig, rig_max, rig_min, z0ig, Lig !LOS PARÁMETROS DE LA SIGMOIDE
        real(fgsl_double), dimension(:), pointer :: ri, sig_rig !EL VECTOR DONDE CALCULAMOS LA SIGMOIDE

        integer :: istat
                
        ub1 = ubound(this%ne,1)
        lb1 = lbound(this%ne,1) 

        ncells = ub1-lb1+1

        accNe = fgsl_interp_accel_alloc()
        splineNe = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineNe,this%timePlasma,this%ne)

        accNuei = fgsl_interp_accel_alloc()
        splineNuei = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineNuei,this%timePlasma,this%nuei)

        accTe = fgsl_interp_accel_alloc()
        splineTe = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineTe,this%timePlasma,this%Te)

        accZ = fgsl_interp_accel_alloc()
        splineZ = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineZ,this%timePlasma,this%Zbar)

        accPU = fgsl_interp_accel_alloc()
        splinePU = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splinePU,this%timePlasma,this%totPU)

        accPL = fgsl_interp_accel_alloc()
        splinePL = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splinePL,this%timePlasma,this%totPL)

        accDPU = fgsl_interp_accel_alloc()
        splineDPU = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineDPU,this%timePlasma,this%totDPU)

        accDPL = fgsl_interp_accel_alloc()
        splineDPL = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineDPL,this%timePlasma,this%totDPL)

        accRPU = fgsl_interp_accel_alloc()
        splineRPU = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineRPU,this%timePlasma,this%totRPU)

        accRPL = fgsl_interp_accel_alloc()
        splineRPL = fgsl_spline_alloc(fgsl_interp_akima,ncells)
        status = fgsl_spline_init (splineRPL,this%timePlasma,this%totRPL)

        allocate(ri(1:this%Nz),STAT=istat)
        ri = 0.0
        allocate(sig_rig(1:this%Nz),STAT=istat)
        sig_rig = 0.0
        allocate(expion_cen(1:this%Nz),STAT=istat)
        expion_cen = 0.0
        allocate(rlion_cen(1:this%Nz),STAT=istat)
        rlion_cen = 0.0

        !@EOG sigmoide de la meseta Kr8+
        ri_max = 5e-6;
        ri_min = 5e-6;
        !@EOG unidades inverso de longitud
        ki = 5e6;
        z0i = 3.5e-3;
        Li = ri_max-ri_min !PARÁMETRO AUXILIAR
        !ri(:) = ri_min + Li/(1+exp(ki*(this%spacez(:)-z0i))) !CALCULO DE LA SIGMOIDE

        !@EOG sigmoide de la gaussiana de Kr8+ (sig_r)
        rig_max = 25e-6;
        rig_min = 15e-6;
        !@EOG unidades inverso de longitud
        kig = 6e3;
        z0ig = 6e-3;
        Lig = rig_max-rig_min !PARÁMETRO AUXILIAR
        !sig_rig(:) = rig_min + Lig/(1+exp(kig*(this%spacez(:)-z0ig))) !CALCULO DE LA SIGMOIDE

        !EOG canal (asumo 100 micrometros de ancho del dominio y this%Nx = this%Ny = 100 -> dx = 1 micrometro).
        Rc = 38e-6; !EOG radio del canal
        R0 = 80e-6; !EOG parametro del canal
        Rv = 90e-6; !EOG parametro del canal

        zshift = 0.75; !EOG Posicion relativa del segundo pico de sobreionización (0 inicio del plasma, 0.5 a la mitad, 1.0 al final)
        
        cen_fac = 0.075; !EOG aumento de densidad electronica sobre el valor sin sobreionizacion (cen_fac = 0.0) del segundo pico
        sig_z = 0.25e-3; !Desviación estandar en dirección longitudinal (propagación) de la sobreionización del segundo pico
        sig_r = 15e-6; !Desviación estandar en dirección transversal (radial) de la sobreionización del segundo pico
        !rlion_cen(:) = ri(:); !GUARDAMOS LA SIGMOIDE QUE HEMOS CALCULADO ANTES EN RLION_CEN

        r_L = 5e-6;     !SIGMOIDE DE VALOR CONSTANTE
        r_u = 20e-6;    !RADIO QUE SEPARA LOS DOS TRAMOS DE LA FUNCION A TROZOS
        sigma_L = 15e-6 !DESVIACION ESTANDAR DEL PRIMER TRAMO (radius<=r_u)
        sigma_u = 17e-6 !DESVIACION ESTANDER DEL SEGUNDO TRAMO (radius>r_u)
        
        z0_fac = 0.15;    !EOG aumento de densidad electronica sobre el valor sin sobreionizacion (cen_fac = 0.0) al inicio del plasma
        sig_z0 = 2e-3;    !Desviación estandar en dirección longitudinal (propagación) de la sobreionización del primer pico
        sig_r0 = 15e-6;   !Desviación estandar en dirección transversal (radial) de la sobreionización del primer pico
        rlion_z0 = 40.e-6;
        
        z0_flag = 0.0;
        cen_flag = 0.0;
        
        if(z0_fac .gt. 0.0)then
           z0_flag = 1.0;
        endif

        if(cen_fac .gt. 0.0)then
           cen_flag = 1.0;
        endif

        do k = 1,this%Nx
           do l = 1,this%Ny
           !N2
#if 0
              !EOG forma longitudinal
              if(this%pflag)then
                 sigmaPz = this%plength/(2*(2*log(2.0))**(1/(2*this%nhgauss)))
                 !EOG N2 perfil de densidad de neutros 
                 this%NZL3D(k,l,:) = exp(-0.5*((this%spacez(:)-this%pcenter)/sigmaPz)**(2*this%nhgauss))*exp(-0.5*((k-this%Nx/2.0)/(this%Nx/5))**2-0.5*((l-this%Ny/2.0)/(this%Nx/5))**2)
              else
                 !EOG N2 perfil de densidad de neutros 
                 !this%NZL3D(k,l,:) = exp(-0.5*((k-this%Nx/2.0)/(this%Nx/5))**2-0.5*((l-this%Ny/2.0)/(this%Nx/5))**2)
                 this%NZL3D(k,l,:) = exp(-0.5*((k-this%Nx/2.0)/(this%Nx/2))**2-0.5*((l-this%Ny/2.0)/(this%Nx/2))**2)
              endif
              this%neprof3D(k,l,:) = this%NZL3D(k,l,:)
#endif           
              !EOG canal (asumo 100 micrometros de ancho del dominio).
#if 1
              !EOG fórmula 4
              radius = sqrt((this%spacex(k)-0.5*this%Lx)**2+(this%spacey(l)-0.5*this%Ly)**2)
              if(radius <= Rc)then
                 this%NZL3D(k,l,:) = 1 + (radius/R0)**2
              else
                 this%NZL3D(k,l,:) = (1+(Rc/R0)**2)*((Rv-radius)/(Rv-Rc))
              endif
              
              !EOG con sobreionizacion fórmula 5              
              this%neprof3D(k,l,:) = this%NZL3D(k,l,:)*(1.0 + cen_fac*exp(-0.5*((this%spacez(:)-zshift*this%Lz)/sig_z)**2)* &
                   exp(-0.5*((radius/sig_r)**2)) + z0_fac*exp(-0.5*(this%spacez(:)/sig_z0)**2)*exp(-0.5*(radius/sig_r0)**2))     
                 
              if(radius <= r_u)then      !LA FUNCION EXPONENCIAL A TROZOS
                 expion_cen = exp(-0.5*(max(radius,r_L)/sigma_L)**2)/exp(-0.5*(r_L/sigma_L)**2)
              else
                 A_ru = (exp(-0.5*(max(r_u,r_L)/sigma_L)**2)/exp(-0.5*(r_L/sigma_L)**2))*exp(0.5*(r_u/sigma_u)**2)
                 expion_cen = A_ru*exp(-0.5*(radius/sigma_u)**2)
              endif          

              !expion_cen(:) = exp(-0.5*(max(radius,rlion_cen(:))/sig_r)**2)
              !expion_cen(:) = expion_cen(:)/exp(-0.5*(rlion_cen(:)/sig_r)**2)

              !expion_cen(:) = exp(-0.5*(max(radius,rlion_cen(:))/sig_rig(:))**2)
              !expion_cen(:) = expion_cen(:)/exp(-0.5*(rlion_cen(:)/sig_rig(:))**2)
              
              expion_z0 = exp(-0.5*(max(radius,rlion_z0)/sig_r0)**2)
              expion_z0 = expion_z0/exp(-0.5*(rlion_z0/sig_r0)**2)

              !EOG con sobreionizacion fórmula 6
              this%NZL3D(k,l,:) = this%NZL3D(k,l,:)*(1.0 - cen_flag*exp(-0.5*((this%spacez(:)-zshift*this%Lz)/sig_z)**2)* &
                   exp(-0.5*(radius/sig_r)**2))*expion_cen*(1.0 - z0_flag*exp(-0.5*(this%spacez(:)/sig_z0)**2)*expion_z0) 
#endif
              !HOMOGENEO
#if 0
              this%NZL3D(k,l,:) = 1.0
              this%neprof3D(k,l,:) = this%NZL3D(k,l,:)
#endif
              !exponencial decreciente (plata)
#if 0
              !this%NZL3D(k,l,:) = 1.6333e21*exp(-0.084730*this%spacex(k)*1e6)
              this%neprof3D(k,l,:) = exp(-0.084730*this%spacex(k)*1e6)
              this%NZL3D(k,l,:) = exp(-0.5*((this%spacex(k)-15e-6)/(sqrt(2.0)*5e-6))**2-0.5*((this%spacey(l)-15e-6)/(sqrt(2.0)*5e-6))**2)
#endif
           enddo
        enddo
        
        if(this%IR_group_vel_flag)then

           avgNe = 0.0
           counter = 0
           do m = 1,this%Nz
              if((this%spacez(m) <= this%zIR) .and. (this%spacez(m) >= this%zIR - 2*this%sigmazIR))then
                 avgNe = avgNe + this%ne(m)
                 counter = counter + 1
              endif
           enddo
           avgNe = avgNe/max(counter,1)

           avgOmega = 5.64e4*sqrt(avgNe)
           avgVgIR = C_SPEED*sqrt(1 - (avgOmega/this%IROmega)**2)
           
           do m = 1,this%Nz
              itime = max(ctime - this%spacez(m)/avgVgIR,this%timePlasma(1))
              this%ne3D(:,:,m) = fgsl_spline_eval(splineNe,itime,accNe)*this%neprof3D(:,:,m)
              this%nuei3D(:,:,m) = fgsl_spline_eval(splineNuei,itime,accNuei)*this%neprof3D(:,:,m)
              this%Te3D(:,:,m) = fgsl_spline_eval(splineTe,itime,accTe)
              this%Zbar3D(:,:,m) = fgsl_spline_eval(splineZ,itime,accZ)
              this%PU3D(:,:,m) = (fgsl_spline_eval(splinePU,itime,accPU)*this%ne3D(:,:,m)*1e-18 &
      &     + fgsl_spline_eval(splineRPU,itime,accRPU))*this%nneutr*1e6*this%NZL3D(:,:,m)        
              this%PL3D(:,:,m) = (fgsl_spline_eval(splinePL,itime,accPL)*this%ne3D(:,:,m)*1e-18 &
      &     + fgsl_spline_eval(splineRPL,itime,accRPL))*this%nneutr*1e6*this%NZL3D(:,:,m)                
              this%DPU3D(:,:,m) = fgsl_spline_eval(splineDPU,itime,accDPU)*this%ne3D(:,:,m)*1e-18 &
      &                         + this%totRDPU
              this%DPL3D(:,:,m) = fgsl_spline_eval(splineDPL,itime,accDPL)*this%ne3D(:,:,m)*1e-18 &
      &                         + this%totRDPL
                   
        enddo

        this%zIR = this%zIR + avgVgIR*dtstep

        else
           
        do m = 1,this%Nz
           itime = max(ctime - this%spacez(m)/C_SPEED,this%timePlasma(1))
           this%ne3D(:,:,m) = fgsl_spline_eval(splineNe,itime,accNe)*this%neprof3D(:,:,m)
           this%nuei3D(:,:,m) = fgsl_spline_eval(splineNuei,itime,accNuei)*this%neprof3D(:,:,m)
           this%Te3D(:,:,m) = fgsl_spline_eval(splineTe,itime,accTe)
           this%Zbar3D(:,:,m) = fgsl_spline_eval(splineZ,itime,accZ)
           this%PU3D(:,:,m) = (fgsl_spline_eval(splinePU,itime,accPU)*this%ne3D(:,:,m)*1e-18 &
     &     + fgsl_spline_eval(splineRPU,itime,accRPU))*this%nneutr*1e6*this%NZL3D(:,:,m)        
           this%PL3D(:,:,m) = (fgsl_spline_eval(splinePL,itime,accPL)*this%ne3D(:,:,m)*1e-18 &
     &     + fgsl_spline_eval(splineRPL,itime,accRPL))*this%nneutr*1e6*this%NZL3D(:,:,m)                
           this%DPU3D(:,:,m) = fgsl_spline_eval(splineDPU,itime,accDPU)*this%ne3D(:,:,m)*1e-18 &
     &                         + this%totRDPU
           this%DPL3D(:,:,m) = fgsl_spline_eval(splineDPL,itime,accDPL)*this%ne3D(:,:,m)*1e-18 &
     &                         + this%totRDPL
                   
        enddo

     endif
     
        call fgsl_spline_free(splineNe)
        call fgsl_interp_accel_free(accNe)

        call fgsl_spline_free(splineNuei)
        call fgsl_interp_accel_free(accNuei)

        call fgsl_spline_free(splineTe)
        call fgsl_interp_accel_free(accTe)

        call fgsl_spline_free(splineZ)
        call fgsl_interp_accel_free(accZ)
                        
        call fgsl_spline_free(splinePU)
        call fgsl_interp_accel_free(accPU)

        call fgsl_spline_free(splinePL)
        call fgsl_interp_accel_free(accPL)

        call fgsl_spline_free(splineDPU)
        call fgsl_interp_accel_free(accDPU)

        call fgsl_spline_free(splineDPL)
        call fgsl_interp_accel_free(accDPL)

        call fgsl_spline_free(splineRPU)
        call fgsl_interp_accel_free(accRPU)

        call fgsl_spline_free(splineRPL)
        call fgsl_interp_accel_free(accRPL)

        deallocate(ri,STAT=istat)
        deallocate(expion_cen,STAT=istat)
        deallocate(rlion_cen,STAT=istat)
                
      end subroutine fillPlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine outputPlasma(this)

        implicit none

        class(t_plasma3D), intent(inout) :: this

        call this%pWriteOut%write3Dreal(this%ne3D,1,1.0d0,"electronDensity",15)
        call this%pWriteOut%write3Dreal(this%Te3D,1,1.0d0,"electronTemperature",19)
        !call this%pWriteOut%write3Dreal(this%Zbar3D,1,1.0d0,"Zbar",4)
        call this%pWriteOut%write3Dreal(this%NU3D,1,1.0d0,"NUpper",6)
        call this%pWriteOut%write3Dreal(this%NL3D,1,1.0d0,"NLower",6)
        call this%pWriteOut%write3Dreal(this%NZL3D,1,1.0d0,"neutralDensity",14)
        !call this%pWriteOut%write3Dreal(this%PU3D,1,1.0d0,"PUpper",6)
        !call this%pWriteOut%write3Dreal(this%PL3D,1,1.0d0,"PLower",6)
        !call this%pWriteOut%write3Dreal(this%DPU3D,1,1.0d0,"DUpper",6)
        !call this%pWriteOut%write3Dreal(this%DPL3D,1,1.0d0,"DLower",6)
        call this%pWriteOut%write3Dreal(this%nuei3D,1,1.0d0,"nuei",4)
        
      end subroutine outputPlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initiatePlasma(this, writeOutAddress)

        implicit none

        class(t_plasma3D), intent(inout) :: this
        class(t_writeOut), pointer, intent(in) :: writeOutAddress

        integer :: x,y,z

        this%pWriteOut => writeOutAddress

        this%dx = this%Lx/this%Nx
        this%dy = this%Ly/this%Ny
        this%dz = this%Lz/this%Nz
        
        do x = 1,this%Nx
           this%spacex(x) = x*this%dx - 0.5*this%dx
        enddo
        do y = 1,this%Ny
           this%spacey(y) = y*this%dy - 0.5*this%dy
        enddo
        do z = 1,this%Nz
           this%spacez(z) = z*this%dz - 0.5*this%dz
        enddo
        
      end subroutine initiatePlasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rates(this)

        implicit none

        class(t_plasma3D), intent(inout) :: this

        real :: auxD
        real, dimension(this%nlevel) :: auxP
        integer :: k,l

        !depopulation rates
        this%totRDPU = 0.0
        do k=1,this%ulevel-1
           this%totRDPU = this%totRDPU + this%radRates(this%ulevel,k)
        enddo

        this%totRDPL = 0.0
        do k=1,this%llevel-1
           this%totRDPL = this%totRDPL + this%radRates(this%llevel,k)
        enddo

        !population rates
        l = 1
        do k=1,this%ntsteps
           do
              if(this%timePlasma(k) < this%timeRates(l)) exit
              l = l+1
           enddo
           auxP(1:this%nlevel) = this%popU(l-1,1:this%nlevel) + (this%popU(l,1:this%nlevel) &
     &           - this%popU(l-1,1:this%nlevel))/(this%timeRates(l)-this%timeRates(l-1)) &
     &           * (this%timePlasma(k)-this%timeRates(l-1))

           this%totPU(k) = this%totPU(k) + sum(auxP(1:this%nlevel)*this%levels(k,1:this%nlevel)) &
     &                     - 0.5*auxP(this%ulevel)*this%levels(k,this%ulevel)
           
       this%totRPU(k) = sum(this%radRates(1:this%nlevel,this%ulevel)*this%levels(k,1:this%nlevel)) - &
     &                      this%radRates(this%ulevel,this%ulevel)*this%levels(k,this%ulevel)
           
           auxP(1:this%nlevel) = this%popL(l-1,1:this%nlevel) + (this%popL(l,1:this%nlevel) &
     &           - this%popL(l-1,1:this%nlevel))/(this%timeRates(l)-this%timeRates(l-1)) &
     &           * (this%timePlasma(k)-this%timeRates(l-1)) 

           this%totPL(k) = this%totPL(k) + sum(auxP(1:this%nlevel)*this%levels(k,1:this%nlevel)) &
     &                     - 0.5*auxP(this%llevel)*this%levels(k,this%llevel)

           !EOG RAD-RATE la contribucion del nivel superior al inferior hay que meterla en el solver
       this%totRPL(k) = sum(this%radRates(1:this%nlevel,this%llevel)*this%levels(k,1:this%nlevel)) - &
            &               this%radRates(this%llevel,this%llevel)*this%levels(k,this%llevel) - &
            &               this%radRates(this%ulevel,this%llevel)*this%levels(k,this%ulevel)

           auxD = this%depopU(l-1) + (this%depopU(l) &
     &           - this%depopU(l-1))/(this%timeRates(l)-this%timeRates(l-1)) &
     &           * (this%timePlasma(k)-this%timeRates(l-1))

           this%totDPU(k) = this%totDPU(k) + auxD

           auxD = this%depopL(l-1) + (this%depopL(l) &
     &           - this%depopL(l-1))/(this%timeRates(l)-this%timeRates(l-1)) &
     &           * (this%timePlasma(k)-this%timeRates(l-1))
           
           this%totDPL(k) = this%totDPL(k) + auxD
        enddo
        
      end subroutine rates

    end module m_plasma3D
