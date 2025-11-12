program make_cmd

  use iso_eep_support
  use iso_eep_color

  implicit none

  integer :: ierr
  type(isochrone_set) :: s
  logical, parameter :: do_timing = .false.
  character(len=file_path) :: BC_table_list = 'bc_table.list', cmd_suffix = ''
  logical :: set_fixed_Fe_div_H = .false.
  real(sp) :: extinction_Av, extinction_Rv, Fe_div_H
  integer :: count_rate, time(3)
  ierr=0
  
  if(do_timing) call system_clock(time(1),count_rate)

  call cmd_init(ierr)

  s% cmd_suffix = cmd_suffix

  if(do_timing) call system_clock(time(2),count_rate)

  if(ierr==0) call write_cmds_to_file(s)

  if(do_timing) call system_clock(time(3),count_rate)

  if(do_timing) call report_timing

contains

  subroutine report_timing
    real(dp) :: t(3),c
    t=real(time,kind=dp)
    c=real(count_rate,kind=dp)
    t=t/c
    write(0,*) ' Total execution time: ', t(3)-t(1)
    write(0,*) ' Time to init,read   : ', t(2)-t(1)
    write(0,*) ' Time to write isos  : ', t(3)-t(2)
  end subroutine report_timing

  subroutine cmd_init(ierr)
    integer, intent(out) :: ierr
    integer :: i, j, n_arg
    character(len=file_path) :: phot_string, arg, option, result
    logical ::  spot_flag

    Fe_div_H = 0.0
    extinction_Av = 0.0
    extinction_Rv = 3.1
    spot_flag = .false.
    
    if(command_argument_count()<1)then
       write(*,*) ' make_cmd:   '
       write(*,*) '   usage: ./make_cmd phot_string isochrone_file [Av] [Rv]  '
       write(*,*) '     phot_string = UBVRIplus, etc.                         '
       write(*,*) '     isochrone_file = name of isochrone file to transform  '
       write(*,*) '     OPTIONAL -                                            '
       write(*,*) '     [Av] extinction in V band                             '
       write(*,*) '     [Rv] selective extinction Rv = Av/E(B-V)              '
       ierr=-1
       return
    endif

    call get_command_argument(1,phot_string)
    call get_command_argument(2,s% filename)

    if(cmd_suffix == '') cmd_suffix=phot_string

    do i=1, command_argument_count()
       call get_command_argument(i,arg)
    enddo
    
    !process command arguments
    n_arg = command_argument_count()
    if(n_arg > 2) then
       do i=3,n_arg
          call get_command_argument(i,arg)
          j=index(arg,'=')
          option=arg(1:j-1)
          result=arg(j+1:)
          if(trim(option)=='Av')then
             read(result,*) extinction_Av
          elseif(trim(option)=='Rv')then
             read(result,*) extinction_Rv
          elseif(trim(option)=="SPOTS")then
             if(trim(result)=='true' .or. trim(result)=='TRUE') spot_flag=.true.
          endif
       enddo
    endif
    
    if(ierr==0) call read_isochrone_file(s,ierr)

    if(ierr/=0) write(*,*) 'read_isochrone_file: ierr = ', ierr
    
    s% Av = extinction_Av
    s% Rv = extinction_Rv
    do i=1,s% number_of_isochrones
       s% iso(i)% Av = s% Av
       s% iso(i)% Rv = s% Rv
    enddo

    if(spot_flag) call spotify

    call color_init(phot_string,BC_table_list,set_fixed_Fe_div_H,Fe_div_H,ierr)

  end subroutine cmd_init

  subroutine spotify
    real(dp) :: Teff, R, logg, g, alfa, beta, Tnew, Rnew, gnew, Tspot
    real(dp) :: scaling_factor
    integer :: i, ilogTeff, ilogR, ilogg, j

    write(0,*) 'adding spots!'

    ilogTeff = 0; ilogR = 0; ilogg = 0
    scaling_factor=0.0_dp

    do i = 1, s% iso(1)% ncol
       if(trim(adjustl(s% iso(1)% cols(i)% name))== 'log_Teff') ilogTeff = i
       if(trim(adjustl(s% iso(1)% cols(i)% name))== 'log_R') ilogR = i
       if(trim(adjustl(s% iso(1)% cols(i)% name))== 'log_g') ilogg = i
    enddo

    do j = 1, s% number_of_isochrones
       do i = 1, s% iso(j)% neep
          Teff = pow10(s% iso(j)% data(ilogTeff,i))
          R = pow10(s% iso(j)% data(ilogR,i))
          logg = s% iso(j)% data(ilogg,i)
          g = pow10(logg)
          !calculate spot parameters
          beta = fd(logg)*spot_frac(Teff)
          alfa = 1.0d0 - beta
          Tspot = scaling_factor*Teff
          !caculate adjusted Teff, R, log(g)
          Tnew = pow(alfa*powi(Teff,4) + beta*powi(Tspot,4), 0.25d0)
          Rnew = R * powi(Teff/Tnew,2)
          gnew = g * powi(Tnew/Teff,4)
          !put the results back into the track
          s% iso(j)% data(ilogR,i) = log10(Rnew)
          s% iso(j)% data(ilogTeff,i) = log10(Tnew)
          s% iso(j)% data(ilogg,i) = log10(g)
       enddo
    enddo
    
  end subroutine spotify

 function spot_frac(Teff) result(y)
    real(dp), intent(in) :: Teff
    real(dp) :: y, x
    real(dp), parameter :: A = -1.77514793d-7
    real(dp), parameter :: B = 1.63313609d-3
    real(dp), parameter :: C = -3.35621302d0
    real(dp), parameter :: xmin = 3.0988893333198866d3
    real(dp), parameter :: xmax = 6.1011106351334465d3
    x = Teff + 400.0_dp
    
    if(x > xmin .and. x < xmax) then
       y = A*x*x + B*x + C
    else
       y = 0.0d0
    endif
  end function spot_frac

  function fd(x) result(y)
    real(dp), intent(in) :: x
    real(dp), parameter :: x0 = 4.0d0
    real(dp), parameter :: sigma = 0.25d0 !0.2d0
    real(dp) :: y
    y = 1.0d0 - 1.0d0/(1.0d0 + exp( (x-x0)/sigma))
  end function fd

end program make_cmd
