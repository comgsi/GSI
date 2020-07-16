module gsi_rfv3io_mod
!$$$   module documentation block
!             .      .    .                                       .
! module:     gsi_rfv3io_mod
!   prgmmr:
!
! abstract: IO routines for regional FV3
!
! program history log:
!   2017-03-08  parrish - create new module gsi_rfv3io_mod, starting from
!                           gsi_nemsio_mod as a pattern.
!   2017-10-10  wu      - setup A grid and interpolation coeff in generate_anl_grid
!   2018-02-22  wu      - add subroutines for read/write fv3_ncdf
!   2019-03-13  cTong   - Port CAPS radar DA capability.
!
! subroutines included:
!   sub gsi_rfv3io_get_grid_specs
!   sub read_fv3_files 
!   sub read_fv3_netcdf_guess
!   sub gsi_fv3ncdf2d_read
!   sub gsi_fv3ncdf_read
!   sub gsi_fv3ncdf_readuv
!   sub wrfv3_netcdf
!   sub gsi_fv3ncdf_writeuv
!   sub gsi_fv3ncdf_writeps
!   sub gsi_fv3ncdf_write
!   sub check
!
! variable definitions:
!
! attributes:
!   langauge: f90
!    machine:
!
!$$$ end documentation block

  use kinds, only: r_kind,i_kind
  use gridmod, only: nlon_regional,nlat_regional
  implicit none

!    directory names (hardwired for now)

  character(len=*),parameter :: grid_spec='fv3_grid_spec'            
  character(len=*),parameter :: ak_bk='fv3_akbk'
  character(len=*),parameter :: dynvars='fv3_dynvars'
  character(len=*),parameter :: tracers='fv3_tracer'
  character(len=*),parameter :: sfcdata='fv3_sfcdata'
  integer(i_kind) gfile_grid_spec,gfile_ak_bk,gfile_dynvars
  integer(i_kind) gfile_tracers,gfile_sfcdata
  integer(i_kind) :: gfile
  save gfile
! --- CAPS ---
  character(len=*),parameter :: anadynvars='ana_dynvars'
  character(len=*),parameter :: anatracers='ana_tracer'
! --- CAPS ---


  integer(i_kind) nx,ny,nz
  real(r_kind),allocatable:: grid_lon(:,:),grid_lont(:,:),grid_lat(:,:),grid_latt(:,:)
  real(r_kind),allocatable:: ak(:),bk(:)
  integer(i_kind),allocatable:: ijns2d(:),displss2d(:),ijns(:),displss(:)
  integer(i_kind),allocatable:: ijnz(:),displsz_g(:)

! set default to private
  private
! set subroutines to public
  public :: gsi_rfv3io_get_grid_specs
  public :: read_fv3_files 
  public :: read_fv3_netcdf_guess
  public :: wrfv3_netcdf

  public :: mype_u,mype_v,mype_t,mype_q,mype_p,mype_oz,mype_ql
  public :: mype_qi,mype_qr,mype_qs,mype_qg,mype_qnr,mype_w ! CAPS
  public :: k_slmsk,k_tsea,k_vfrac,k_vtype,k_stype,k_zorl,k_smc,k_stc
  public :: k_snwdph,k_f10m,mype_2d,n2d,k_orog,k_psfc
  public :: ijns,ijns2d,displss,displss2d,ijnz,displsz_g

  integer(i_kind) mype_u,mype_v,mype_t,mype_q,mype_p,mype_oz,mype_ql
  integer(i_kind) mype_qi,mype_qr,mype_qs,mype_qg,mype_qnr,mype_w ! CAPS
  integer(i_kind) k_slmsk,k_tsea,k_vfrac,k_vtype,k_stype,k_zorl,k_smc,k_stc
  integer(i_kind) k_snwdph,k_f10m,mype_2d,n2d,k_orog,k_psfc

contains

subroutine gsi_rfv3io_get_grid_specs(grid_spec,ak_bk,ierr)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    gsi_rfv3io_get_grid_specs
!   pgrmmr: parrish     org: np22                date: 2017-04-03
!
! abstract:  obtain grid dimensions nx,ny and grid definitions
!                grid_x,grid_xt,grid_y,grid_yt,grid_lon,grid_lont,grid_lat,grid_latt
!                nz,ak(nz),bk(nz)
!
! program history log:
!   2017-04-03  parrish - initial documentation
!   2017-10-10  wu - setup A grid and interpolation coeff with generate_anl_grid
!   2018-02-16  wu - read in time info from file coupler.res
!                    read in lat, lon at the center and corner of the grid cell
!                    from file fv3_grid_spec, and vertical grid infor from file fv3_akbk
!                    setup A grid and interpolation/rotation coeff
!   input argument list:
!    grid_spec
!    ak_bk
!    lendian_out
!
!   output argument list:
!    ierr
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
  use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
  use netcdf, only: nf90_inquire_variable
  use mpimod, only: mype
  use mod_fv3_lola, only: generate_anl_grid
  use gridmod,  only:nsig,regional_time,regional_fhr,aeta1_ll,aeta2_ll
  use gridmod,  only:nlon_regional,nlat_regional,eta1_ll,eta2_ll
  use kinds, only: i_kind,r_kind
  use constants, only: half,zero
  use mpimod, only: mpi_comm_world,mpi_itype,mpi_rtype

  implicit none

  character(*),   intent(in   ) :: grid_spec
  character(*),   intent(in   ) :: ak_bk
  integer(i_kind),intent(  out) :: ierr
  integer(i_kind) i,k,ndimensions,iret,nvariables,nattributes,unlimiteddimid
  integer(i_kind) len,gfile_loc
  character(len=128) :: name
  integer(i_kind) myear,mmonth,mday,mhour,mminute,msecond
  real(r_kind),allocatable:: abk_fv3(:)

!!!!! set regional_time
    open(24,file='coupler.res',form='formatted')
    read(24,*)
    read(24,*)
    read(24,*)myear,mmonth,mday,mhour,mminute,msecond
    close(24)
    if(mype==0)  write(6,*)' myear,mmonth,mday,mhour,mminute,msecond=', myear,mmonth,mday,mhour,mminute,msecond
    regional_time(1)=myear
    regional_time(2)=mmonth
    regional_time(3)=mday
    regional_time(4)=mhour
    regional_time(5)=mminute
    regional_time(6)=msecond
    regional_fhr=zero          ! forecast hour set zero for now

!!!!!!!!!!    grid_spec  !!!!!!!!!!!!!!!
    ierr=0
    iret=nf90_open(trim(grid_spec),nf90_nowrite,gfile_grid_spec)
    if(iret/=nf90_noerr) then
       write(6,*)' problem opening ',trim(grid_spec),', Status = ',iret
       ierr=1
       return
    endif

    iret=nf90_inquire(gfile_grid_spec,ndimensions,nvariables,nattributes,unlimiteddimid)
    gfile_loc=gfile_grid_spec
    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       if(trim(name)=='grid_xt') nx=len
       if(trim(name)=='grid_yt') ny=len
    enddo
    nlon_regional=nx
    nlat_regional=ny
    if(mype==0)write(6,*),'nx,ny=',nx,ny

!!!    get nx,ny,grid_lon,grid_lont,grid_lat,grid_latt,nz,ak,bk

    allocate(grid_lat(nx+1,ny+1))
    allocate(grid_lon(nx+1,ny+1))
    allocate(grid_latt(nx,ny))
    allocate(grid_lont(nx,ny))

    do k=ndimensions+1,nvariables
       iret=nf90_inquire_variable(gfile_loc,k,name,len)
       if(trim(name)=='grid_lat') then
          iret=nf90_get_var(gfile_loc,k,grid_lat)
       endif
       if(trim(name)=='grid_lon') then
          iret=nf90_get_var(gfile_loc,k,grid_lon)
       endif
       if(trim(name)=='grid_latt') then
          iret=nf90_get_var(gfile_loc,k,grid_latt)
       endif
       if(trim(name)=='grid_lont') then
          iret=nf90_get_var(gfile_loc,k,grid_lont)
       endif
    enddo

    iret=nf90_close(gfile_loc)

    iret=nf90_open(ak_bk,nf90_nowrite,gfile_loc)
    if(iret/=nf90_noerr) then
       write(6,*)' problem opening ',trim(ak_bk),', Status = ',iret
       ierr=1
       return
    endif
    iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       if(trim(name)=='xaxis_1') nz=len
    enddo
    if(mype==0)write(6,'(" nz=",i5)') nz

    nsig=nz-1

!!!    get ak,bk

    allocate(aeta1_ll(nsig),aeta2_ll(nsig))
    allocate(eta1_ll(nsig+1),eta2_ll(nsig+1))
    allocate(ak(nz),bk(nz),abk_fv3(nz))

    do k=ndimensions+1,nvariables
       iret=nf90_inquire_variable(gfile_loc,k,name,len)
       if(trim(name)=='ak'.or.trim(name)=='AK') then
          iret=nf90_get_var(gfile_loc,k,abk_fv3)
          do i=1,nz
             ak(i)=abk_fv3(nz+1-i)
          enddo
       endif
       if(trim(name)=='bk'.or.trim(name)=='BK') then
          iret=nf90_get_var(gfile_loc,k,abk_fv3)
          do i=1,nz
             bk(i)=abk_fv3(nz+1-i)
          enddo
       endif
    enddo
    iret=nf90_close(gfile_loc)

!!!!! change unit of ak 
    do i=1,nsig+1
       eta1_ll(i)=ak(i)*0.001_r_kind
       eta2_ll(i)=bk(i)
    enddo
    do i=1,nsig
       aeta1_ll(i)=half*(ak(i)+ak(i+1))*0.001_r_kind
       aeta2_ll(i)=half*(bk(i)+bk(i+1))
    enddo
    if(mype==0)then
       do i=1,nz
          write(6,'(" ak,bk(",i3,") = ",2f17.6)') i,ak(i),bk(i)
       enddo
    endif

!!!!!!! setup A grid and interpolation/rotation coeff.
    call generate_anl_grid(nx,ny,grid_lon,grid_lont,grid_lat,grid_latt)

    deallocate (grid_lon,grid_lat,grid_lont,grid_latt)
    deallocate (ak,bk,abk_fv3)

    return
end subroutine gsi_rfv3io_get_grid_specs

subroutine read_fv3_files(mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_fv3_files
!   prgmmr: wu               org: np22                date: 2017-10-10
!
! abstract: read in from fv3 files and figure out available time levels 
!           of background fields starting from read_files as a pattern
!           temporary setup for one first guess file
! program history log:
!
!   input argument list:
!     mype     - pe number
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block

    use kinds, only: r_kind,i_kind
    use mpimod, only: mpi_comm_world,ierror,mpi_rtype,npe
    use guess_grids, only: nfldsig,nfldsfc,ntguessig,ntguessfc,&
         ifilesig,ifilesfc,hrdifsig,hrdifsfc,create_gesfinfo
    use guess_grids, only: hrdifsig_all,hrdifsfc_all
    use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen,nhr_assimilation
    use constants, only: zero,one,r60inv
    use obsmod, only: iadate,time_offset
    use gridmod, only:regional_time
    implicit none

! Declare passed variables
    integer(i_kind),intent(in   ) :: mype

! Declare local parameters
    real(r_kind),parameter:: r0_001=0.001_r_kind

! Declare local variables
    logical(4) fexist
    character(6) filename
    integer(i_kind) in_unit
    integer(i_kind) i,j,iwan,npem1
    integer(i_kind) nhr_half
    integer(i_kind) nminanl,nmings,nming2,ndiff,isecond
    integer(i_kind),dimension(4):: idateg
    integer(i_kind),dimension(5):: idate5
    real(r_kind) hourg,temp,t4dv
    real(r_kind),dimension(202,2):: time_ges

!-----------------------------------------------------------------------------
! Start read_nems_nmmb_files here.
    nhr_half=nhr_assimilation/2
    if(nhr_half*2 < nhr_assimilation) nhr_half=nhr_half+1
    npem1=npe-1

    do i=1,202
       time_ges(i,1) = 999_r_kind
       time_ges(i,2) = 999_r_kind
    end do


! Let a single task query the guess files.
    if(mype==npem1) then

!    Convert analysis time to minutes relative to fixed date
       call w3fs21(iadate,nminanl)
       write(6,*)'READ_netcdf_fv3_FILES:  analysis date,minutes ',iadate,nminanl

!    Check for consistency of times from sigma guess files.
       in_unit=15
       iwan=0
!WWWWWW setup for one first guess file for now
!      do i=0,9 !place holder for FGAT
       i=3

!wwww read in from the external file directly, no internal files sigfxx for FV3
          idate5(1)=  regional_time(1)
          idate5(2)=  regional_time(2)
          idate5(3)=  regional_time(3)
          idate5(4)=  regional_time(4)
          idate5(5)=  regional_time(5)
          isecond  =  regional_time(6)
          hourg    =  zero ! forcast hour

          call w3fs21(idate5,nmings)
          nming2=nmings+60*hourg
          write(6,*)'READ_netcdf_fv3_FILES:  sigma guess file, nming2 ',hourg,idate5,nming2
          t4dv=real((nming2-iwinbgn),r_kind)*r60inv
          if (l4dvar.or.l4densvar) then
             if (t4dv<zero .OR. t4dv>winlen) then
                write(6,*)'ges file not in time range, t4dv=',t4dv
!               cycle ! place holder for FGAT
             endif
          else
             ndiff=nming2-nminanl
!for test with the 3 hr files with FGAT
             if(abs(ndiff) > 60*nhr_half ) then
                write(6,*)'ges file not in time range, ndiff=',ndiff
!               cycle ! place holder for FGAT
             endif
          endif
          iwan=iwan+1
          time_ges(iwan,1) =real((nming2-iwinbgn),r_kind)*r60inv
          time_ges(iwan+100,1)=i+r0_001
!       end do ! i !place holder for FGAT
       time_ges(201,1)=one
       time_ges(202,1)=one
       if(iwan > 1)then
          do i=1,iwan
             do j=i+1,iwan
                if(time_ges(j,1) < time_ges(i,1))then
                   temp=time_ges(i+100,1)
                   time_ges(i+100,1)=time_ges(j+100,1)
                   time_ges(j+100,1)=temp
                   temp=time_ges(i,1)
                   time_ges(i,1)=time_ges(j,1)
                   time_ges(j,1)=temp
                end if
             end do
             if(abs(time_ges(i,1)-time_offset) < r0_001)time_ges(202,1) = i
          end do
       end if
       time_ges(201,1) = iwan+r0_001

!    Check for consistency of times from surface guess files.
       iwan=0
       do i=0,99
          write(filename,200)i
  200     format('sfcf',i2.2)
          inquire(file=filename,exist=fexist)
          if(fexist)then
             idateg(4)=iadate(1); idateg(2)=iadate(2)
             idateg(3)=iadate(3); idateg(1)=iadate(4)
             hourg = zero
             idate5(1)=idateg(4); idate5(2)=idateg(2)
             idate5(3)=idateg(3); idate5(4)=idateg(1); idate5(5)=0
             call w3fs21(idate5,nmings)
             nming2=nmings+60*hourg
             write(6,*)'READ_netcdf_fv3_FILES:  surface guess file, nming2 ',hourg,idateg,nming2
             ndiff=nming2-nminanl
             if(abs(ndiff) > 60*nhr_half ) then
                write(6,*)'ges file not in time range, ndiff=',ndiff
!               cycle ! place holder for FGAT
             endif
             iwan=iwan+1
             time_ges(iwan,2) =real((nming2-iwinbgn),r_kind)*r60inv
             time_ges(iwan+100,2)=i+r0_001
          end if
          if(iwan==1) exit
       end do
       time_ges(201,2)=one
       time_ges(202,2)=one
       if(iwan > 1)then
          do i=1,iwan
             do j=i+1,iwan
                if(time_ges(j,2) < time_ges(i,2))then
                   temp=time_ges(i+100,2)
                   time_ges(i+100,2)=time_ges(j+100,2)
                   time_ges(j+100,2)=temp
                   temp=time_ges(i,2)
                   time_ges(i,2)=time_ges(j,2)
                   time_ges(j,2)=temp
                end if
             end do
             if(abs(time_ges(i,2)-time_offset) < r0_001)time_ges(202,2) = i
          end do
       end if
       time_ges(201,2) = iwan+r0_001
    end if


! Broadcast guess file information to all tasks
    call mpi_bcast(time_ges,404,mpi_rtype,npem1,mpi_comm_world,ierror)

    nfldsig   = nint(time_ges(201,1))
!!nfldsfc   = nint(time_ges(201,2))
    nfldsfc   = nfldsig

! Allocate space for guess information files
    call create_gesfinfo

    do i=1,nfldsig
       ifilesig(i) = -100
       hrdifsig(i) = zero
    end do

    do i=1,nfldsfc
       ifilesfc(i) = -100
       hrdifsfc(i) = zero
    end do

! Load time information for sigma guess field sinfo into output arrays
    ntguessig = nint(time_ges(202,1))
    do i=1,nfldsig
       hrdifsig(i) = time_ges(i,1)
       ifilesig(i) = nint(time_ges(i+100,1))
       hrdifsig_all(i) = hrdifsig(i)
    end do
    if(mype == 0) write(6,*)'READ_netcdf_fv3_FILES:  sigma fcst files used in analysis  :  ',&
         (ifilesig(i),i=1,nfldsig),(hrdifsig(i),i=1,nfldsig),ntguessig


! Load time information for surface guess field info into output arrays
    ntguessfc = nint(time_ges(202,2))
    do i=1,nfldsfc
       hrdifsfc(i) = time_ges(i,2)
       ifilesfc(i) = nint(time_ges(i+100,2))
       hrdifsfc_all(i) = hrdifsfc(i)
    end do

! Below is a temporary fix. The nems_nmmb regional mode does not have a
! surface
! file.  Instead the surface fields are passed through the atmospheric guess
! file.  Without a separate surface file the code above sets ntguessig and 
! nfldsig to zero.  This causes problems later in the code when arrays for
! the surface fields are allocated --> one of the array dimensions is nfldsfc
! and it will be zero.  This portion of the code should be rewritten, but
! until
! it is, the fix below gets around the above mentioned problem.

    ntguessfc = ntguessig
!!nfldsfc   = nfldsig
    do i=1,nfldsfc
       hrdifsfc(i) = hrdifsig(i)
       ifilesfc(i) = ifilesig(i)
       hrdifsfc_all(i) = hrdifsfc(i)
    end do
    if(mype == 0) write(6,*)'READ_nems_nmb_FILES:  surface fcst files used in analysis:  ',&
         (ifilesfc(i),i=1,nfldsfc),(hrdifsfc(i),i=1,nfldsfc),ntguessfc


! End of routine
    return
end subroutine read_fv3_files

subroutine read_fv3_netcdf_guess
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_fv3_netcdf_guess            read fv3 interface file
!   prgmmr: wu               org: np22                date: 2017-07-06
!
! abstract:  read guess for FV3 regional model
! program history log:
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block
    use kinds, only: r_kind,i_kind
    use mpimod, only: npe
    use guess_grids, only: ges_tsen,ges_prsi
    use gridmod, only: lat2,lon2,nsig,ijn,eta1_ll,ijn_s
    use constants, only: one,fv
    use gsi_metguess_mod, only: gsi_metguess_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use mpeu_util, only: die
    use guess_grids, only: ntguessig
! --- CAPS ---
    use constants, only: zero, one_tenth
    use mpimod, only: mype
    use caps_radaruse_mod, only: l_use_log_qx, l_use_log_qx_pval, l_use_log_nt
    use caps_radaruse_mod, only: l_use_dbz_caps
    use caps_radaruse_mod, only: l_cvpnr, cvpnr_pval
! --- CAPS ---

    implicit none

    character(len=24),parameter :: myname = 'read_fv3_netcdf_guess'
    integer(i_kind) k,i,j
    integer(i_kind) it,ier,istatus
    real(r_kind),dimension(:,:),pointer::ges_ps=>NULL()
    real(r_kind),dimension(:,:),pointer::ges_z=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_u=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_v=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_q=>NULL()
!   real(r_kind),dimension(:,:,:),pointer::ges_ql=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_oz=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_tv=>NULL()
! --- CAPS ---
    real(r_kind),dimension(:,:,:),pointer::ges_ql=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_qi=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_qr=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_qs=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_qg=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_qnr=>NULL()
    real(r_kind),dimension(:,:,:),pointer::ges_w=>NULL()

    real(r_kind) :: qr_min, qs_min, qg_min
    real(r_kind) :: qr_thrshd, qs_thrshd, qg_thrshd

    integer,parameter   :: USEZG = 1

! --- CAPS ---


!    setup list for 2D surface fields
    k_f10m =1   !fact10
    k_stype=2   !soil_type
    k_vfrac=3   !veg_frac
    k_vtype=4   !veg_type
    k_zorl =5   !sfc_rough
    k_tsea =6   !sfct ?
    k_snwdph=7  !sno ?
    k_stc  =8   !soil_temp
    k_smc  =9   !soil_moi
    k_slmsk=10  !isli
    k_orog =11  !terrain
    n2d=11

    if (l_use_dbz_caps) then ! CAPS
       if(npe< 14) then
          call die('read_fv3_netcdf_guess','not enough PEs to read in fv3 fields' )
       endif
       mype_u=0
       mype_v=1
       mype_t=2
       mype_p=3
       mype_q=4
       mype_ql=5
       mype_qi=6
       mype_qr=7
       mype_qs=8
       mype_qg=9
       mype_oz=10
       mype_2d=11
       mype_qnr=12
       mype_w=13
    else
       if(npe< 8) then
          call die('read_fv3_netcdf_guess','not enough PEs to read in fv3 fields' )
       endif
       mype_u=0           
       mype_v=1
       mype_t=2
       mype_p=3
       mype_q=4
       mype_ql=5
       mype_oz=6
       mype_2d=7 
    end if
      
    allocate(ijns(npe),ijns2d(npe),ijnz(npe) )
    allocate(displss(npe),displss2d(npe),displsz_g(npe) )

    do i=1,npe
       ijns(i)=ijn_s(i)*nsig
       ijnz(i)=ijn(i)*nsig
       ijns2d(i)=ijn_s(i)*n2d 
    enddo
    displss(1)=0
    displsz_g(1)=0
    displss2d(1)=0
    do i=2,npe
       displss(i)=displss(i-1)+ ijns(i-1)
       displsz_g(i)=displsz_g(i-1)+ ijnz(i-1)
       displss2d(i)=displss2d(i-1)+ ijns2d(i-1)
    enddo

!   do it=1,nfldsig
    it=ntguessig


    ier=0
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ps' ,ges_ps ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'z' , ges_z ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'u' , ges_u ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'v' , ges_v ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'tv' ,ges_tv ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'q'  ,ges_q ,istatus );ier=ier+istatus
!    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql'  ,ges_ql ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'oz'  ,ges_oz ,istatus );ier=ier+istatus
! --- CAPS ---
    if (l_use_dbz_caps) then
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql' ,ges_ql ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qi' ,ges_qi ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qr' ,ges_qr ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qs' ,ges_qs ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qg' ,ges_qg ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qnr',ges_qnr ,istatus );ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'w' , ges_w ,istatus );ier=ier+istatus
    end if
! --- CAPS ---
    if (ier/=0) call die(trim(myname),'cannot get pointers for fv3 met-fields, ier =',ier)

    call gsi_fv3ncdf_readuv(ges_u,ges_v)
    call gsi_fv3ncdf_read(dynvars,'T','t',ges_tsen(1,1,1,it),mype_t)
    call gsi_fv3ncdf_read(dynvars,'DELP','delp',ges_prsi,mype_p)
    ges_prsi(:,:,nsig+1,it)=eta1_ll(nsig+1)
    do i=nsig,1,-1
       ges_prsi(:,:,i,it)=ges_prsi(:,:,i,it)*0.001_r_kind+ges_prsi(:,:,i+1,it)
    enddo
    ges_ps(:,:)=ges_prsi(:,:,1,it)
    call gsi_fv3ncdf_read(tracers,'SPHUM','sphum',ges_q,mype_q)
!    call gsi_fv3ncdf_read(tracers,'LIQ_WAT','liq_wat',ges_ql,mype_ql)
    call gsi_fv3ncdf_read(tracers,'O3MR','o3mr',ges_oz,mype_oz)
! --- CAPS ---
    if (l_use_dbz_caps) then
       call gsi_fv3ncdf_read(dynvars,'W','w',ges_w,mype_w)
       call gsi_fv3ncdf_read(tracers,'LIQ_WAT','liq_wat',ges_ql,mype_ql)
       call gsi_fv3ncdf_read(tracers,'ICE_WAT','ice_wat',ges_qi,mype_qi)
       call gsi_fv3ncdf_read(tracers,'RAINWAT','rainwat',ges_qr,mype_qr)
       call gsi_fv3ncdf_read(tracers,'SNOWWAT','snowwat',ges_qs,mype_qs)
       call gsi_fv3ncdf_read(tracers,'GRAUPEL','graupel',ges_qg,mype_qg)
       call gsi_fv3ncdf_read(tracers,'RAIN_NC','rain_nc',ges_qnr,mype_qnr)

! --- CAPS ---
!!  tsen2tv  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k=1,nsig
        do i=1,lon2
          do j=1,lat2
             ges_tv(j,i,k)=ges_tsen(j,i,k,it)*(one+fv*ges_q(j,i,k))
!            Convert from qx to log(qx) if option is ON 
             if (l_use_log_qx) then ! CAPS
                if ( l_use_log_qx_pval .eq. 0. ) then ! CVlogq
                   if (ges_tsen(j,i,k,it) .gt. 274.15) then
                      qr_min=2.9E-6_r_kind
                      qr_thrshd=qr_min * one_tenth
                      qs_min=0.1E-9_r_kind
                      qs_thrshd=qs_min
                      qg_min=3.1E-7_r_kind
                      qg_thrshd=qg_min * one_tenth
                   else if (ges_tsen(j,i,k,it) .le. 274.15 .and. ges_tsen(j,i,k,it) .ge. 272.15) then
                      qr_min=2.0E-6_r_kind
                      qr_thrshd=qr_min * one_tenth
                      qs_min=1.3E-7_r_kind
                      qs_thrshd=qs_min * one_tenth
                      qg_min=3.1E-7_r_kind
                      qg_thrshd=qg_min * one_tenth
                   else if (ges_tsen(j,i,k,it) .lt. 272.15) then
                      qr_min=0.1E-9_r_kind
                      qr_thrshd=qr_min
                      qs_min=6.3E-6_r_kind
                      qs_thrshd=qs_min * one_tenth
                      qg_min=3.1E-7_r_kind
                      qg_thrshd=qg_min * one_tenth
                   end if

                   if ( ges_qr(j,i,k) .le. qr_thrshd )  ges_qr(j,i,k) = qr_min
                   if ( ges_qs(j,i,k) .le. qs_thrshd )  ges_qs(j,i,k) = qs_min
                   if ( ges_qg(j,i,k) .le. qg_thrshd )  ges_qg(j,i,k) = qg_min

!                  convert hydrometer variables (qr,qs and qg) to log
                   if(mype==0 .AND. i==3 .AND. j==3 .AND. k==3)then
                       write(6,*)'read_fv3_netcdf_guess: ',     &
                           ' reset zero of qr/qs/qg to specified values (~0dbz) ', &
                           'before log transformation. (for dbz assimilation)'
                       write(6,*)'read_fv3_netcdf_guess: convert qr/qs/qg to log.'
                   end if
                   ges_qr(j,i,k) = log(ges_qr(j,i,k))
                   ges_qs(j,i,k) = log(ges_qs(j,i,k))
                   ges_qg(j,i,k) = log(ges_qg(j,i,k))

                else if ( l_use_log_qx_pval .gt. 0. ) then   ! CVpq
!                  convert cloud hydrometer variables (qr,qs and qg) to log
                   if(mype==0 .AND. i==3 .AND. j==3 .AND. k==3)then
                       write(6,*)'read_wrf_mass_guess_netcdf: ',                        &
                                 ' reset zero of qr/qs/qg to specified values(~0dbz) ', &
                                 'before log transformation. (for dbz assimilation)'
                       write(6,*)'read_wrf_mass_guess_netcdf: convert qr/qs/qg to log.'
                   end if
                   !RONG KONG modified below to tranform the ges_qr in
                   !logq space (Gang's version) to q space (arps version)
                   !so that we don't need to make change to the origial background field....
                   if(USEZG==1)then
                      ! chenll
                   !  ges_qr(j,i,k) = log(max(ges_qr(j,i,k),1.0E-6_r_kind))
                   !  ges_qs(j,i,k) = log(max(ges_qs(j,i,k),1.0E-6_r_kind))
                   !  ges_qg(j,i,k) = log(max(ges_qg(j,i,k),1.0E-6_r_kind))
                      ges_qr(j,i,k)=((max(ges_qr(j,i,k),1.0E-6_r_kind))**l_use_log_qx_pval-1)/l_use_log_qx_pval
                      ges_qs(j,i,k)=((max(ges_qs(j,i,k),1.0E-6_r_kind))**l_use_log_qx_pval-1)/l_use_log_qx_pval
                      ges_qg(j,i,k)=((max(ges_qg(j,i,k),1.0E-6_r_kind))**l_use_log_qx_pval-1)/l_use_log_qx_pval

                     !cliu for test 0 hydrometeor mixing ratio background

                    ! ges_qr(j,i,k) = log(1.0E-6_r_kind)
                    ! ges_qs(j,i,k) = log(1.0E-6_r_kind)
                    ! ges_qg(j,i,k) = log(1.0E-6_r_kind)
                   else
                      ges_qr(j,i,k) = ges_qr(j,i,k) + 1.0E-6_r_kind
                      ges_qs(j,i,k) = ges_qs(j,i,k) + 1.0E-6_r_kind
                      ges_qg(j,i,k) = ges_qg(j,i,k) + 1.0E-6_r_kind
                   endif
                   !Rong Kong commented off above!
                end if

             else ! CVq
!               2018/02 comment off the following part of code
!               Following Dr. Xue's suggestion, 
!               temperature is set up based on Dr.s Chengsi Liu and Rong Kong's work in ARPS-3DVar
!               (Match the values in prewgt_reg and wrwrfmassa)

!               if (ges_tsen(j,i,k,it) .gt. 278.15) then
!                  qr_min=1.0E-5_r_kind
!                  qr_thrshd=qr_min * one_tenth
!                  qs_min=0.1E-9_r_kind
!                  qs_thrshd=qs_min
!                  qg_min=1.0E-7_r_kind
!                  qg_thrshd=qg_min * one_tenth
!               else if (ges_tsen(j,i,k,it) .le. 278.15 .and. ges_tsen(j,i,k,it) .ge. 243.15) then
!                  qs_min=1.0E-7_r_kind
!                  qs_thrshd=qs_min * one_tenth
!                  qg_min=6.0E-7_r_kind
!                  qg_thrshd=qg_min * one_tenth
!                  if (ges_tsen(j,i,k,it) .le. 278.15 .and. ges_tsen(j,i,k,it) .ge. 268.15) then
!                     qr_min=5.0E-6_r_kind
!                     qr_thrshd=qr_min * one_tenth
!                  else if (ges_tsen(j,i,k,it) .lt. 268.15 .and. ges_tsen(j,i,k,it) .ge. 243.15) then
!                     qr_min=0.1E-9_r_kind
!                     qr_thrshd=qr_min * one_tenth
!                  end if
!               else if (ges_tsen(j,i,k,it) .lt. 243.15) then
!                  qr_min=0.1E-9_r_kind
!                  qr_thrshd=qr_min
!                  qs_min=9.3E-6_r_kind
!                  qs_thrshd=qs_min * one_tenth
!                  qg_min=7.1E-7_r_kind
!                  qg_thrshd=qg_min * one_tenth
!               end if

!               if ( ges_qr(j,i,k) .le. qr_thrshd )  ges_qr(j,i,k) = qr_min
!               if ( ges_qs(j,i,k) .le. qs_thrshd )  ges_qs(j,i,k) = qs_min
!               if ( ges_qg(j,i,k) .le. qg_thrshd )  ges_qg(j,i,k) = qg_min

                qr_min=zero
                qs_min=zero
                qg_min=zero
                if (mype==0 .AND. i==3 .AND. j==3 .AND. k==3) then
                    write(6,*)'read_wrf_mass_guess: only reset (qr/qs/qg) to 0.0 for negative analysis value. (regular qx)'
                end if
                ges_qr(j,i,k) = max(ges_qr(j,i,k), qr_min)
                ges_qs(j,i,k) = max(ges_qs(j,i,k), qs_min)
                ges_qg(j,i,k) = max(ges_qg(j,i,k), qg_min)
             end if

!            Treatment on qnr
!            log-transform to QNRAIN (g.zhao)
             if ( l_use_log_nt ) then
                ges_qnr(j,i,k) = log(max(ges_qnr(j,i,k),1.0E-2_r_kind))
             end if
             if (l_cvpnr) then
                ges_qnr(j,i,k)=((max(ges_qnr(j,i,k),1.0E-2_r_kind)**cvpnr_pval)-1)/cvpnr_pval !added by cliu
             endif
          enddo
        enddo
      enddo
    else ! in the case of not using l_use_dbz_caps
!!  tsen2tv  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do k=1,nsig
          do j=1,lon2
             do i=1,lat2
                ges_tv(i,j,k)=ges_tsen(i,j,k,it)*(one+fv*ges_q(i,j,k))
             enddo
          enddo
       enddo
    end if

    call gsi_fv3ncdf2d_read(it,ges_z)

end subroutine read_fv3_netcdf_guess

subroutine gsi_fv3ncdf2d_read(it,ges_z)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsi_fv3ncdf2d_read       
!   prgmmr: wu w             org: np22                date: 2017-10-17
!
! abstract: read in 2d fields from fv3_sfcdata file in mype_2d 
!                Scatter the field to each PE 
! program history log:
!   input argument list:
!     it    - time index for 2d fields
!
!   output argument list:
!     ges_z - surface elevation
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block
    use kinds, only: r_kind,i_kind
    use mpimod, only: ierror,mpi_comm_world,npe,mpi_rtype,mype
    use guess_grids, only: fact10,soil_type,veg_frac,veg_type,sfc_rough, &
         sfct,sno,soil_temp,soil_moi,isli
    use gridmod, only: lat2,lon2,itotsub,ijn_s
    use general_commvars_mod, only: ltosi_s,ltosj_s
    use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
    use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
    use netcdf, only: nf90_inquire_variable
    use mod_fv3_lola, only: fv3_h_to_ll,nxa,nya
    use constants, only: grav

    implicit none

    integer(i_kind),intent(in) :: it   
    real(r_kind),intent(in),dimension(:,:),pointer::ges_z
    character(len=128) :: name
    integer(i_kind),allocatable,dimension(:):: dim_id,dim
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:):: a
    real(r_kind),allocatable,dimension(:,:,:):: sfcn2d
    real(r_kind),allocatable,dimension(:,:,:):: sfc
    real(r_kind),allocatable,dimension(:,:):: sfc1
    integer(i_kind) iret,gfile_loc,i,k,len,ndim
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid
    integer(i_kind) kk,n,ns,j,ii,jj,mm1

    mm1=mype+1
    allocate(a(nya,nxa))
    allocate(work(itotsub*n2d))
    allocate( sfcn2d(lat2,lon2,n2d))

 if(mype==mype_2d ) then
    iret=nf90_open(sfcdata,nf90_nowrite,gfile_loc)
    if(iret/=nf90_noerr) then
       write(6,*)' problem opening ',trim(sfcdata),', Status = ',iret
       return
    endif
    iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
    allocate(dim(ndimensions))
    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       dim(k)=len
    enddo
!!!!!!!!!!!! read in 2d variables !!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=ndimensions+1,nvariables
       iret=nf90_inquire_variable(gfile_loc,i,name,len)
       if( trim(name)=='f10m'.or.trim(name)=='F10M' ) then
          k=k_f10m
       else if( trim(name)=='stype'.or.trim(name)=='STYPE' ) then
          k=k_stype
       else if( trim(name)=='vfrac'.or.trim(name)=='VFRAC' ) then
          k=k_vfrac
       else if( trim(name)=='vtype'.or.trim(name)=='VTYPE' ) then
          k=k_vtype
       else if( trim(name)=='zorl'.or.trim(name)=='ZORL' ) then
          k=k_zorl
       else if( trim(name)=='tsea'.or.trim(name)=='TSEA' ) then
          k=k_tsea
       else if( trim(name)=='sheleg'.or.trim(name)=='SHELEG' ) then
          k=k_snwdph
       else if( trim(name)=='stc'.or.trim(name)=='STC' ) then
          k=k_stc 
       else if( trim(name)=='smc'.or.trim(name)=='SMC' ) then
          k=k_smc
       else if( trim(name)=='SLMSK'.or.trim(name)=='slmsk' ) then
          k=k_slmsk
       else
          cycle 
       endif
       iret=nf90_inquire_variable(gfile_loc,i,ndims=ndim)
       if(allocated(dim_id    )) deallocate(dim_id    )
       allocate(dim_id(ndim))
       iret=nf90_inquire_variable(gfile_loc,i,dimids=dim_id)
       if(allocated(sfc       )) deallocate(sfc       )
       allocate(sfc(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
       iret=nf90_get_var(gfile_loc,i,sfc)
       call fv3_h_to_ll(sfc(:,:,1),a,nx,ny,nxa,nya)

       kk=0
       do n=1,npe
          ns=displss2d(n)+(k-1)*ijn_s(n)
          do j=1,ijn_s(n)
             ns=ns+1
             kk=kk+1
             ii=ltosi_s(kk)
             jj=ltosj_s(kk)
             work(ns)=a(ii,jj)
          end do
       end do
    enddo ! i
    iret=nf90_close(gfile_loc)

!!!! read in orog from dynam !!!!!!!!!!!!
    iret=nf90_open(trim(dynvars ),nf90_nowrite,gfile_loc)
    if(iret/=nf90_noerr) then
       write(6,*)' problem opening ',trim(dynvars ),gfile_loc,', Status = ',iret
       return
    endif

    iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
    if(allocated(dim    )) deallocate(dim    )
    allocate(dim(ndimensions))

    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       dim(k)=len
    enddo


    do k=ndimensions+1,nvariables
       iret=nf90_inquire_variable(gfile_loc,k,name,len)
       if(trim(name)=='PHIS'   .or. trim(name)=='phis'  ) then
          iret=nf90_inquire_variable(gfile_loc,k,ndims=ndim)
          if(allocated(dim_id    )) deallocate(dim_id    )
          allocate(dim_id(ndim))
          iret=nf90_inquire_variable(gfile_loc,k,dimids=dim_id)
          allocate(sfc1(dim(dim_id(1)),dim(dim_id(2))) )
          iret=nf90_get_var(gfile_loc,k,sfc1)
          exit
       endif
    enddo     !   k
    iret=nf90_close(gfile_loc)

    k=k_orog
    call fv3_h_to_ll(sfc1,a,nx,ny,nxa,nya)

    kk=0
    do n=1,npe
       ns=displss2d(n)+(k-1)*ijn_s(n)
       do j=1,ijn_s(n)
          ns=ns+1
          kk=kk+1
          ii=ltosi_s(kk)
          jj=ltosj_s(kk)
          work(ns)=a(ii,jj)
       end do
    end do

    deallocate (dim_id,sfc,sfc1,dim)
 endif  ! mype


!!!!!!! scatter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_scatterv(work,ijns2d,displss2d,mpi_rtype,&
      sfcn2d,ijns2d(mm1),mpi_rtype,mype_2d,mpi_comm_world,ierror)

    deallocate ( work )

    fact10(:,:,it)=sfcn2d(:,:,k_f10m)
    soil_type(:,:,it)=sfcn2d(:,:,k_stype)
    veg_frac(:,:,it)=sfcn2d(:,:,k_vfrac)
    veg_type(:,:,it)=sfcn2d(:,:,k_vtype)
    sfc_rough(:,:,it)=sfcn2d(:,:,k_zorl)
    sfct(:,:,it)=sfcn2d(:,:,k_tsea)
    sno(:,:,it)=sfcn2d(:,:,k_snwdph)
    soil_temp(:,:,it)=sfcn2d(:,:,k_stc)
    soil_moi(:,:,it)=sfcn2d(:,:,k_smc)
    ges_z(:,:)=sfcn2d(:,:,k_orog)/grav
    isli(:,:,it)=nint(sfcn2d(:,:,k_slmsk))
    deallocate (sfcn2d,a)
    return
end subroutine gsi_fv3ncdf2d_read

subroutine gsi_fv3ncdf_read(filename,varname,varname2,work_sub,mype_io)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsi_fv3ncdf_read       
!   prgmmr: wu               org: np22                date: 2017-10-10
!
! abstract: read in a field from a netcdf FV3 file in mype_io
!          then scatter the field to each PE 
! program history log:
!
!   input argument list:
!     filename    - file name to read from       
!     varname     - variable name to read in
!     varname2    - variable name to read in
!     mype_io     - pe to read in the field
!
!   output argument list:
!     work_sub    - output sub domain field
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block


    use kinds, only: r_kind,i_kind
    use mpimod, only: ierror,mpi_comm_world,npe,mpi_rtype,mype
    use gridmod, only: lat2,lon2,nsig,nlat,nlon,itotsub,ijn_s
    use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
    use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
    use netcdf, only: nf90_inquire_variable
    use mod_fv3_lola, only: fv3_h_to_ll
    use general_commvars_mod, only: ltosi_s,ltosj_s

    implicit none

    character(*)   ,intent(in   ) :: varname,varname2,filename
    real(r_kind)   ,intent(out  ) :: work_sub(lat2,lon2,nsig) 
    integer(i_kind)   ,intent(in   ) :: mype_io
    character(len=128) :: name
    real(r_kind),allocatable,dimension(:,:,:):: uu
    integer(i_kind),allocatable,dimension(:):: dim_id,dim
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:):: a


    integer(i_kind) n,ns,k,len,ndim
    integer(i_kind) gfile_loc,iret
    integer(i_kind) nz,nzp1,kk,j,mm1,i,ir,ii,jj
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    mm1=mype+1
    allocate (work(itotsub*nsig))

    if(mype==mype_io ) then
       iret=nf90_open(trim(filename),nf90_nowrite,gfile_loc)
       if(iret/=nf90_noerr) then
          write(6,*)' problem opening ',trim(filename),gfile_loc,', Status = ',iret
          return
       endif

       iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
       allocate(dim(ndimensions))
       allocate(a(nlat,nlon))

       do k=1,ndimensions
          iret=nf90_inquire_dimension(gfile_loc,k,name,len)
          dim(k)=len
       enddo


       do k=ndimensions+1,nvariables
          iret=nf90_inquire_variable(gfile_loc,k,name,len)
          if(trim(name)==varname .or. trim(name)==varname2) then
             iret=nf90_inquire_variable(gfile_loc,k,ndims=ndim)
             if(allocated(dim_id    )) deallocate(dim_id    )
             allocate(dim_id(ndim))
             iret=nf90_inquire_variable(gfile_loc,k,dimids=dim_id)
             if(allocated(uu        )) deallocate(uu        )
             allocate(uu(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
             iret=nf90_get_var(gfile_loc,k,uu)
             exit
          endif
       enddo     !   k
       nz=nsig
       nzp1=nz+1
       do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(uu(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work(ns)=a(ii,jj)
             end do
          end do
       enddo ! i

       iret=nf90_close(gfile_loc)
       deallocate (uu,a,dim,dim_id)


    endif !mype

    call mpi_scatterv(work,ijns,displss,mpi_rtype,&
       work_sub,ijns(mm1),mpi_rtype,mype_io,mpi_comm_world,ierror)

    deallocate (work)
    return
end subroutine gsi_fv3ncdf_read

subroutine gsi_fv3ncdf_readuv(ges_u,ges_v)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsi_fv3ncdf_readuv
!   prgmmr: wu w             org: np22                date: 2017-11-22
!
! abstract: read in a field from a netcdf FV3 file in mype_u,mype_v
!           then scatter the field to each PE 
! program history log:
!
!   input argument list:
!
!   output argument list:
!     ges_u       - output sub domain u field
!     ges_v       - output sub domain v field
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block
    use kinds, only: r_kind,i_kind
    use mpimod, only: ierror,mpi_comm_world,npe,mpi_rtype,mype
    use gridmod, only: lat2,lon2,nsig,itotsub,ijn_s
    use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
    use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
    use netcdf, only: nf90_inquire_variable
    use mod_fv3_lola, only: fv3_h_to_ll,nya,nxa,fv3uv2earth
    use general_commvars_mod, only: ltosi_s,ltosj_s

    implicit none

    real(r_kind)   ,intent(out  ) :: ges_u(lat2,lon2,nsig) 
    real(r_kind)   ,intent(out  ) :: ges_v(lat2,lon2,nsig) 
    character(len=128) :: name
    real(r_kind),allocatable,dimension(:,:,:):: uu,temp1
    integer(i_kind),allocatable,dimension(:):: dim_id,dim
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:):: a
    real(r_kind),allocatable,dimension(:,:):: u,v

    integer(i_kind) n,ns,k,len,ndim
    integer(i_kind) gfile_loc,iret
    integer(i_kind) nz,nzp1,kk,j,mm1,i,ir,ii,jj
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    allocate (work(itotsub*nsig))
    mm1=mype+1
    if(mype==mype_u .or. mype==mype_v) then
       iret=nf90_open(dynvars,nf90_nowrite,gfile_loc)
       if(iret/=nf90_noerr) then
          write(6,*)' problem opening ',trim(dynvars),', Status = ',iret
          return
       endif

       iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)

       allocate(dim(ndimensions))
       allocate(a(nya,nxa))

       do k=1,ndimensions
          iret=nf90_inquire_dimension(gfile_loc,k,name,len)
          dim(k)=len
       enddo

       allocate(u(dim(1),dim(4)))
       allocate(v(dim(1),dim(4)))

       do k=ndimensions+1,nvariables
          iret=nf90_inquire_variable(gfile_loc,k,name,len)
          if(trim(name)=='u'.or.trim(name)=='U' .or.  &
             trim(name)=='v'.or.trim(name)=='V' ) then 
             iret=nf90_inquire_variable(gfile_loc,k,ndims=ndim)
             if(allocated(dim_id    )) deallocate(dim_id    )
             allocate(dim_id(ndim))
             iret=nf90_inquire_variable(gfile_loc,k,dimids=dim_id)
!    NOTE:   dimension of variables on native fv3 grid.  
!            u and v have an extra row in one of the dimensions
             if(allocated(uu)) deallocate(uu)
             allocate(uu(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
             iret=nf90_get_var(gfile_loc,k,uu)
             if(trim(name)=='u'.or.trim(name)=='U') then
                allocate(temp1(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
                temp1=uu
             else if(trim(name)=='v'.or.trim(name)=='V') then
                exit
             endif
          endif
       enddo     !   k
! transfor to earth u/v, interpolate to analysis grid, reverse vertical order
       nz=nsig
       nzp1=nz+1
       do i=1,nz
          ir=nzp1-i 
          call fv3uv2earth(temp1(:,:,i),uu(:,:,i),nx,ny,u,v)
          if(mype==mype_u)then
             call fv3_h_to_ll(u,a,nx,ny,nxa,nya)
          else
             call fv3_h_to_ll(v,a,nx,ny,nxa,nya)
          endif
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work(ns)=a(ii,jj)
             end do
          end do
       enddo ! i

       deallocate(temp1,a)
       deallocate (dim,dim_id,uu,v,u)

       iret=nf90_close(gfile_loc)
    endif ! mype

!!  scatter to ges_u,ges_v !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_scatterv(work,ijns,displss,mpi_rtype,&
         ges_u,ijns(mm1),mpi_rtype,mype_u,mpi_comm_world,ierror)
    call mpi_scatterv(work,ijns,displss,mpi_rtype,&
         ges_v,ijns(mm1),mpi_rtype,mype_v,mpi_comm_world,ierror)
    deallocate(work)
end subroutine gsi_fv3ncdf_readuv

subroutine wrfv3_netcdf
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    wrfv3_netcdf           write out FV3 analysis increments
!   prgmmr: wu               org: np22                date: 2017-10-23
!
! abstract:  write FV3 analysis  in netcdf format
!
! program history log:
!   2019-04-18  cTong - import CAPS radar DA capabilities 
!   2019-11-22  cTong - modify "add_saved" to properly output analyses
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    use kinds, only: r_kind,i_kind
    use guess_grids, only: ntguessig,ges_tsen
    use gsi_metguess_mod, only: gsi_metguess_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use mpeu_util, only: die
! --- CAPS ---
    use gridmod, only: nsig,lon2,lat2, pt_ll, aeta1_ll
    use mpimod, only: mype
    use constants, only: zero,one_tenth,one,r0_01
    use constants, only: r10, r100, rd
    use caps_radaruse_mod, only: l_use_log_qx, l_use_log_qx_pval, cld_nt_updt, l_use_log_nt
    use caps_radaruse_mod, only: l_use_dbz_caps
    use caps_radaruse_mod, only: init_mm_qnr
    use caps_radaruse_mod, only: l_cvpnr, cvpnr_pval
    use gridmod, only: pt_ll,aeta1_ll

! --- CAPS ---
    implicit none

! Declare local constants
    logical add_saved
! variables for cloud info
    integer(i_kind) ier,istatus,it
    real(r_kind),pointer,dimension(:,:  ):: ges_ps  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_u   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_v   =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_q   =>NULL()
! --- CAPS ----
    integer(i_kind) i,j,k
! variables for log q conversion
    real(r_kind) :: qr_min, qs_min, qg_min
    real(r_kind) :: qr_thrshd, qs_thrshd, qg_thrshd
    real(r_kind) :: qr_tmp, qs_tmp, qg_tmp
    real(r_kind) :: qnr_tmp
    real(r_kind),pointer,dimension(:,:,:):: ges_ql  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qi  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qr  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qs  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qg  =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_qnr =>NULL()
    real(r_kind),pointer,dimension(:,:,:):: ges_w   =>NULL()
! temporary arrays for dbz-related variables
    real(r_kind) :: P1D,T1D,Q1D,RHO,QR1D
    integer,parameter   :: USEZG = 1
    real(r_kind),parameter:: D608=0.608_r_kind ! CAPS
    real(r_kind),dimension(lat2,lon2,nsig):: tmparr_qr, tmparr_qs, tmparr_qg, tmparr_qnr
! --- CAPS ----

    it=ntguessig
    ier=0
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ps' ,ges_ps ,istatus );ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'u' , ges_u ,istatus);ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'v' , ges_v ,istatus);ier=ier+istatus
    call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'q'  ,ges_q ,istatus);ier=ier+istatus
! --- CAPS ---
    if (l_use_dbz_caps) then
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql' ,ges_ql,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qi' ,ges_qi,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qr' ,ges_qr,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qs' ,ges_qs,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qg' ,ges_qg,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qnr',ges_qnr,istatus);ier=ier+istatus
       call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'w' , ges_w ,istatus);ier=ier+istatus
    end if
! --- CAPS ---
    if (ier/=0) call die('get ges','cannot get pointers for fv3 met-fields, ier =',ier)

! --- CAPS ---
    if (l_use_dbz_caps) then
! Initialized temp arrays with ges. Will be recalculated later if log(q) is used
       tmparr_qr =ges_qr
       tmparr_qs =ges_qs
       tmparr_qg =ges_qg
       tmparr_qnr=ges_qnr

       do k=1,nsig
         do i=1,lon2
           do j=1,lat2

!             initialize hydrometeors as zero
              qr_tmp=zero
              qs_tmp=zero
              qg_tmp=zero

              if ( l_use_log_qx ) then 
                 if ( l_use_log_qx_pval .eq. 0. ) then ! CVlogq

                    if (ges_tsen(j,i,k,it) .gt. 274.15) then
                       qr_min=2.9E-6_r_kind
                       qr_thrshd=qr_min * one_tenth
                       qs_min=0.1E-9_r_kind
                       qs_thrshd=qs_min
                       qg_min=3.1E-7_r_kind
                       qg_thrshd=qg_min * one_tenth
                    else if (ges_tsen(j,i,k,it) .le. 274.15 .and. ges_tsen(j,i,k,it) .ge. 272.15) then
                       qr_min=2.0E-6_r_kind
                       qr_thrshd=qr_min * one_tenth
                       qs_min=1.3E-7_r_kind
                       qs_thrshd=qs_min * one_tenth
                       qg_min=3.1E-7_r_kind
                       qg_thrshd=qg_min * one_tenth
                    else if (ges_tsen(j,i,k,it) .lt. 272.15) then
                       qr_min=0.1E-9_r_kind
                       qr_thrshd=qr_min
                       qs_min=6.3E-6_r_kind
                       qs_thrshd=qs_min * one_tenth
                       qg_min=3.1E-7_r_kind
                       qg_thrshd=qg_min * one_tenth
                    end if

!                   convert log(qr/qs/qg) back to qr/qs/qg
                    if (mype==0 .AND. i==3 .AND. j==3 .AND. k==3) then
                       write(6,*)'wrfv3_netcdf: convert log(qr/qs/qg) back to qr/qs/qg.'
                       write(6,*)'wrfv3_netcdf: then reset (qr/qs/qg) to 0.0 for some cases.'
                    end if

                    qr_tmp=exp(ges_qr(j,i,k))
                    qs_tmp=exp(ges_qs(j,i,k))
                    qg_tmp=exp(ges_qg(j,i,k))

!                   if no update or very tiny value of qr/qs/qg, re-set/clear it off to zero
                    if ( abs(qr_tmp - qr_min) .lt. (qr_min*r0_01) ) then
                       qr_tmp=zero
                    else if (qr_tmp .lt. qr_thrshd) then
                       qr_tmp=zero
                    end if

                    if ( abs(qs_tmp - qs_min) .lt. (qs_min*r0_01) ) then
                       qs_tmp=zero
                    else if (qs_tmp .lt. qs_thrshd) then
                       qs_tmp=zero
                    end if

                    if ( abs(qg_tmp - qg_min) .lt. (qg_min*r0_01) ) then
                       qg_tmp=zero
                    else if (qg_tmp .lt. qg_thrshd) then
                       qg_tmp=zero
                    end if

                    IF (mype==0 .AND. qr_tmp<0) WRITE(*,*)'CCT: qr_tmp,qr_min,qr_thrshd=',qr_tmp,qr_min,qr_thrshd

                 else if ( l_use_log_qx_pval .gt. 0. ) then   ! CVpq
!                   convert log(qr/qs/qg) back to qr/qs/qg
                    if (mype==0 .AND. i==3 .AND. j==3 .AND. k==3) then
                       write(6,*)'wrwrfmass_netcdf: convert log(qr/qs/qg) back to qr/qs/qg.'
                       write(6,*)'wrwrfmass_netcdf: then reset (qr/qs/qg) to 0.0 for some cases.'
                    end if

                    if(USEZG==1)then
                      !qr_tmp=max(exp(ges_qr(j,i,k))-1.0E-6_r_kind,0.0)
                      !qs_tmp=max(exp(ges_qs(j,i,k))-1.0E-6_r_kind,0.0)
                      !qg_tmp=max(exp(ges_qg(j,i,k))-1.0E-6_r_kind,0.0)
                      !chenll
                       qr_tmp=max((l_use_log_qx_pval*ges_qr(j,i,k)+1)**(1/l_use_log_qx_pval)-1.0E-6_r_kind,0.0)
                       qs_tmp=max((l_use_log_qx_pval*ges_qs(j,i,k)+1)**(1/l_use_log_qx_pval)-1.0E-6_r_kind,0.0)
                       qg_tmp=max((l_use_log_qx_pval*ges_qg(j,i,k)+1)**(1/l_use_log_qx_pval)-1.0E-6_r_kind,0.0)
                    else
                       qr_tmp=max(ges_qr(j,i,k)-1.0E-6_r_kind,0.0)
                       qs_tmp=max(ges_qs(j,i,k)-1.0E-6_r_kind,0.0)
                       qg_tmp=max(ges_qg(j,i,k)-1.0E-6_r_kind,0.0)
                    endif

                    !Rong Kong add a upper limit to hydrometeors to disable overshooting
                    if ( 1 == 2 ) then  !threshold test 1
                       qr_tmp=min(qr_tmp,2E-2_r_kind)
                       qs_tmp=min(qs_tmp,2E-2_r_kind)
                       qg_tmp=min(qg_tmp,3E-2_r_kind)
                    else
                       qr_tmp=min(qr_tmp,1E-2_r_kind)
                       qs_tmp=min(qs_tmp,1.0E-2_r_kind)
                       qg_tmp=min(qg_tmp,2E-2_r_kind)
                    endif
                    !Rong Kong added above
                 end if
              
              else   ! For CVq

!                temperature match the setup in read_wrf_mass_guess and
!                prewgt_reg (based on Drs. Chengsi Liu and Rong's work in ARPS-3DVar)
!                if (ges_tsen(j,i,k,it) .gt. 278.15) then
!                   qr_min=1.0E-5_r_kind
!                   qr_thrshd=qr_min * one_tenth
!                   qs_min=0.1E-9_r_kind
!                   qs_thrshd=qs_min
!                   qg_min=1.0E-7_r_kind
!                   qg_thrshd=qg_min * one_tenth
!                else if (ges_tsen(j,i,k,it) .le. 278.15 .and. ges_tsen(j,i,k,it) .ge. 243.15) then
!                   qs_min=1.0E-7_r_kind
!                   qs_thrshd=qs_min * one_tenth
!                   qg_min=6.0E-7_r_kind
!                   qg_thrshd=qg_min * one_tenth
!                   if (ges_tsen(j,i,k,it) .le. 278.15 .and. ges_tsen(j,i,k,it) .ge. 268.15) then
!                      qr_min=5.0E-6_r_kind
!                      qr_thrshd=qr_min * one_tenth
!                   else if (ges_tsen(j,i,k,it) .lt. 268.15 .and. ges_tsen(j,i,k,it) .ge. 243.15) then
!                      qr_min=0.1E-9_r_kind
!                      qr_thrshd=qr_min * one_tenth
!                   end if
!                else if (ges_tsen(j,i,k,it) .lt. 243.15) then
!                   qr_min=0.1E-9_r_kind
!                   qr_thrshd=qr_min
!                   qs_min=9.3E-6_r_kind
!                   qs_thrshd=qs_min * one_tenth
!                   qg_min=7.1E-7_r_kind
!                   qg_thrshd=qg_min * one_tenth
!                end if

                 qr_min=zero
                 qs_min=zero
                 qg_min=zero

                 if (mype==0 .AND. i==3 .AND. j==3 .AND. k==3) then
                     write(6,*)'wrwrfmass_netcdf: only reset (qr/qs/qg) to 0.0 for negative analysis value. (regular qx)'
                 end if

                 qr_tmp=ges_qr(j,i,k)-1.0E-8_r_kind
                 qs_tmp=ges_qs(j,i,k)-1.0E-8_r_kind
                 qg_tmp=ges_qg(j,i,k)-1.0E-8_r_kind


                 qr_tmp=max(qr_tmp,qr_min)
                 qs_tmp=max(qs_tmp,qs_min)
                 qg_tmp=max(qg_tmp,qg_min)

              end if         ! log_qx

              tmparr_qr(j,i,k)=qr_tmp
              tmparr_qs(j,i,k)=qs_tmp
              tmparr_qg(j,i,k)=qg_tmp

!             =========================================================================
!             !
!             re-initialisation of QNRAIN with analyzed qr and N0r(which is single-moment parameter)
!             equation is used in subroutine init_MM of initlib3d.f90 in arps pacakge
              qnr_tmp = zero
              if ( cld_nt_updt .eq. 2 ) then
                 T1D=ges_tsen(j,i,k,it)                                 ! sensible temperature (K)
                 P1D=r100*(aeta1_ll(k)*(r10*ges_ps(j,i)-pt_ll)+pt_ll)   ! pressure hPa --> Pa
                 Q1D=ges_q(j,i,k)/(one-ges_q(j,i,k))                    ! mixing ratio
                 RHO=P1D/(rd*T1D*(one+D608*Q1D))                        ! air density in kg m^-3
                 QR1D=qr_tmp
                 CALL init_mm_qnr(RHO,QR1D,qnr_tmp)
                 qnr_tmp = max(qnr_tmp, zero)

                 if (mype==0 .AND. i==25 .AND. j==25 .AND. k==15) then
                      write(6,*)'wrwrfmass_netcdf: re-init qnrain (number concentration) with analysis qrain and N0_rain. at (25,25,15) & 
                                     qnr,rho,qr,t,qv,pres=',qnr_tmp, RHO, QR1D, T1D, Q1D, P1D
                 end if
              else

!                log_transform to QNRAIN 
                 if (l_use_log_nt) then
                    if(USEZG==1)then
                       qnr_tmp=max(exp(ges_qnr(j,i,k))-1.0E-2_r_kind,0.0)
                    else
                       qnr_tmp=max(ges_qnr(j,i,k)-1.0E-2_r_kind,0.0)
                    endif
                    ! if (qnr_tmp <= one ) qnr_tmp = zero
                 elseif (l_cvpnr) then
                    qnr_tmp=max((cvpnr_pval*ges_qnr(j,i,k)+1)**(1/cvpnr_pval)-1.0E-2_r_kind,0.0)
                 else
                    qnr_tmp=ges_qnr(j,i,k)
                 end if

              end if
              tmparr_qnr(j,i,k)=qnr_tmp

           end do
         end do
       end do

!   write out ! CAPS uses ana+dynvars/tracers
       add_saved=.true.
       call gsi_fv3ncdf_write(anadynvars,'T',ges_tsen(1,1,1,it),mype_t,add_saved)
       call gsi_fv3ncdf_writeuv(ges_u,ges_v,mype_v,add_saved)
       call gsi_fv3ncdf_writeps(anadynvars,'delp',ges_ps,mype_p,add_saved)
       call gsi_fv3ncdf_write(anadynvars,'W',ges_w,mype_w,add_saved)
       call gsi_fv3ncdf_write(anatracers,'sphum',  ges_q ,mype_q ,add_saved)
       add_saved=.true.
       call gsi_fv3ncdf_write(anatracers,'liq_wat',ges_ql,mype_ql,add_saved)
       call gsi_fv3ncdf_write(anatracers,'ice_wat',ges_qi,mype_qi,add_saved)
       call gsi_fv3ncdf_write(anatracers,'rainwat',tmparr_qr,mype_qr,add_saved)
       call gsi_fv3ncdf_write(anatracers,'snowwat',tmparr_qs,mype_qs,add_saved)
       call gsi_fv3ncdf_write(anatracers,'graupel',tmparr_qg,mype_qg,add_saved)
       call gsi_fv3ncdf_write(anatracers,'rain_nc',tmparr_qnr,mype_qnr,add_saved)
    else

       add_saved=.true.

!   write out
       call gsi_fv3ncdf_write(dynvars,'T',ges_tsen(1,1,1,it),mype_t,add_saved)
       call gsi_fv3ncdf_write(tracers,'sphum',ges_q   ,mype_q,add_saved)
       call gsi_fv3ncdf_writeuv(ges_u,ges_v,mype_v,add_saved)
       call gsi_fv3ncdf_writeps(dynvars,'delp',ges_ps,mype_p,add_saved)

    end if

    
end subroutine wrfv3_netcdf

subroutine gsi_fv3ncdf_writeuv(varu,varv,mype_io,add_saved)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    gsi_nemsio_writeuv
!   pgrmmr: wu
!
! abstract: gather u/v fields to mype_io, put u/v in FV3 model defined directions & orders
!           then write out
!
! program history log:
!
!   input argument list:
!    varu,varv
!    add_saved - true: add analysis increments to readin guess then write out
!              - false: write out total analysis fields
!    mype_io
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    use mpimod, only:  mpi_rtype,mpi_comm_world,ierror,npe,mype
    use gridmod, only: lat2,lon2,nlon,nlat,lat1,lon1,nsig, &
                       ijn,displs_g,itotsub,iglobal, &
                       nlon_regional,nlat_regional
    use mod_fv3_lola, only: fv3_ll_to_h,fv3_h_to_ll, &
                            fv3uv2earth,earthuv2fv3
    use general_commvars_mod, only: ltosi,ltosj
    use netcdf, only: nf90_open,nf90_close,nf90_noerr
    use netcdf, only: nf90_write,nf90_inq_varid
    use netcdf, only: nf90_put_var,nf90_get_var
    use caps_radaruse_mod, only: l_use_dbz_caps ! CAPS

    implicit none

    real(r_kind)   ,intent(in   ) :: varu(lat2,lon2,nsig)
    real(r_kind)   ,intent(in   ) :: varv(lat2,lon2,nsig)
    integer(i_kind),intent(in   ) :: mype_io
    logical        ,intent(in   ) :: add_saved

    integer(i_kind) :: ugrd_VarId,gfile_loc,vgrd_VarId
    integer(i_kind) i,j,mm1,n,k,ns,kr,m
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:,:):: work_sub,work_au,work_av
    real(r_kind),allocatable,dimension(:,:,:):: work_bu,work_bv
    real(r_kind),allocatable,dimension(:,:):: u,v,workau2,workav2
    real(r_kind),allocatable,dimension(:,:):: workbu2,workbv2

    mm1=mype+1

    allocate(    work(max(iglobal,itotsub)*nsig),work_sub(lat1,lon1,nsig))
!!!!!! gather analysis u !! revers k !!!!!!!!!!!
    do k=1,nsig
       kr=nsig+1-k
       do i=1,lon1
          do j=1,lat1
             work_sub(j,i,kr)=varu(j+1,i+1,k)
          end do
       end do
    enddo
    call mpi_gatherv(work_sub,ijnz(mm1),mpi_rtype, &
          work,ijnz,displsz_g,mpi_rtype,mype_io,mpi_comm_world,ierror)

    if(mype==mype_io) then
       allocate( work_au(nlat,nlon,nsig),work_av(nlat,nlon,nsig))
       ns=0
       do m=1,npe
          do k=1,nsig
             do n=displs_g(m)+1,displs_g(m)+ijn(m) 
             ns=ns+1
             work_au(ltosi(n),ltosj(n),k)=work(ns)
             end do
          enddo
       enddo
    endif  ! mype

!!!!!! gather analysis v !! reverse k !!!!!!!!!!!!!!!!!!
    do k=1,nsig
       kr=nsig+1-k
       do i=1,lon1
          do j=1,lat1
             work_sub(j,i,kr)=varv(j+1,i+1,k)
          end do
       end do
    enddo
    call mpi_gatherv(work_sub,ijnz(mm1),mpi_rtype, &
          work,ijnz,displsz_g,mpi_rtype,mype_io,mpi_comm_world,ierror)

    if(mype==mype_io) then
       ns=0
       do m=1,npe
          do k=1,nsig
             do n=displs_g(m)+1,displs_g(m)+ijn(m) 
            ns=ns+1
             work_av(ltosi(n),ltosj(n),k)=work(ns)
             end do
          enddo
       enddo
       deallocate(work,work_sub)
       allocate( u(nlon_regional,nlat_regional+1))
       allocate( v(nlon_regional+1,nlat_regional))
       allocate( work_bu(nlon_regional,nlat_regional+1,nsig))
       allocate( work_bv(nlon_regional+1,nlat_regional,nsig))
       
       if ( l_use_dbz_caps )then ! CAPS
          call check( nf90_open(trim(anadynvars ),nf90_write,gfile_loc) )
       else
          call check( nf90_open(trim(dynvars ),nf90_write,gfile_loc) )
       end if 
       call check( nf90_inq_varid(gfile_loc,'u',ugrd_VarId) )
       call check( nf90_inq_varid(gfile_loc,'v',vgrd_VarId) )

       if(add_saved)then
          allocate( workau2(nlat,nlon),workav2(nlat,nlon))
          allocate( workbu2(nlon_regional,nlat_regional+1))
          allocate( workbv2(nlon_regional+1,nlat_regional))
!!!!!!!!  readin work_b !!!!!!!!!!!!!!!!
          call check( nf90_get_var(gfile_loc,ugrd_VarId,work_bu) )
          call check( nf90_get_var(gfile_loc,vgrd_VarId,work_bv) )
          do k=1,nsig
             call fv3uv2earth(work_bu(1,1,k),work_bv(1,1,k),nlon_regional,nlat_regional,u,v)
             call fv3_h_to_ll(u,workau2,nlon_regional,nlat_regional,nlon,nlat)
             call fv3_h_to_ll(v,workav2,nlon_regional,nlat_regional,nlon,nlat)
!!!!!!!! find analysis_inc:  work_a !!!!!!!!!!!!!!!!
             work_au(:,:,k)=work_au(:,:,k)-workau2(:,:)
             work_av(:,:,k)=work_av(:,:,k)-workav2(:,:)
             call fv3_ll_to_h(work_au(:,:,k),u,nlon,nlat,nlon_regional,nlat_regional,.true.)
             call fv3_ll_to_h(work_av(:,:,k),v,nlon,nlat,nlon_regional,nlat_regional,.true.)
             call earthuv2fv3(u,v,nlon_regional,nlat_regional,workbu2,workbv2)
!!!!!!!!  add analysis_inc to readin work_b !!!!!!!!!!!!!!!!
             work_bu(:,:,k)=work_bu(:,:,k)+workbu2(:,:)
             work_bv(:,:,k)=work_bv(:,:,k)+workbv2(:,:)
          enddo
          deallocate(workau2,workbu2,workav2,workbv2)
       else
          do k=1,nsig
             call fv3_ll_to_h(work_au(:,:,k),u,nlon,nlat,nlon_regional,nlat_regional,.true.)
             call fv3_ll_to_h(work_av(:,:,k),v,nlon,nlat,nlon_regional,nlat_regional,.true.)
             call earthuv2fv3(u,v,nlon_regional,nlat_regional,work_bu(:,:,k),work_bv(:,:,k))
          enddo
       endif

       deallocate(work_au,work_av,u,v)
       if ( l_use_dbz_caps )then ! CAPS
          print *,'write out u/v to ',trim(anadynvars )
       else
          print *,'write out u/v to ',trim(dynvars )
       end if 
       call check( nf90_put_var(gfile_loc,ugrd_VarId,work_bu) )
       call check( nf90_put_var(gfile_loc,vgrd_VarId,work_bv) )
       call check( nf90_close(gfile_loc) )
       deallocate(work_bu,work_bv)
    end if !mype_io

    if(allocated(work))deallocate(work)
    if(allocated(work_sub))deallocate(work_sub)

end subroutine gsi_fv3ncdf_writeuv

subroutine gsi_fv3ncdf_writeps(filename,varname,var,mype_io,add_saved)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    gsi_nemsio_writeps
!   pgrmmr: wu
!
! abstract: write out analyzed "delp" to fv_core.res.nest02.tile7.nc
!
! program history log:
!
!   input argument list:
!    varu,varv
!    add_saved
!    mype     - mpi task id
!    mype_io
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    use mpimod, only: mpi_rtype,mpi_comm_world,ierror,mype
    use gridmod, only: lat2,lon2,nlon,nlat,lat1,lon1,nsig
    use gridmod, only: ijn,displs_g,itotsub,iglobal
    use gridmod,  only: nlon_regional,nlat_regional,eta1_ll,eta2_ll
    use mod_fv3_lola, only: fv3_ll_to_h,fv3_h_to_ll
    use general_commvars_mod, only: ltosi,ltosj
    use netcdf, only: nf90_open,nf90_close
    use netcdf, only: nf90_write,nf90_inq_varid
    use netcdf, only: nf90_put_var,nf90_get_var
    implicit none

    real(r_kind)   ,intent(in   ) :: var(lat2,lon2)
    integer(i_kind),intent(in   ) :: mype_io
    logical        ,intent(in   ) :: add_saved
    character(*)   ,intent(in   ) :: varname,filename

    integer(i_kind) :: VarId,gfile_loc
    integer(i_kind) i,j,mm1,k,kr,kp
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:):: work_sub,work_a
    real(r_kind),allocatable,dimension(:,:,:):: work_b,work_bi
    real(r_kind),allocatable,dimension(:,:):: workb2,worka2


    mm1=mype+1
    allocate(    work(max(iglobal,itotsub)),work_sub(lat1,lon1) )
    do i=1,lon1
       do j=1,lat1
          work_sub(j,i)=var(j+1,i+1)
       end do
    end do
    call mpi_gatherv(work_sub,ijn(mm1),mpi_rtype, &
          work,ijn,displs_g,mpi_rtype,mype_io,mpi_comm_world,ierror)

    if(mype==mype_io) then
       allocate( work_a(nlat,nlon))
       do i=1,iglobal
          work_a(ltosi(i),ltosj(i))=work(i)
       end do
       allocate( work_bi(nlon_regional,nlat_regional,nsig+1))
       allocate( work_b(nlon_regional,nlat_regional,nsig))
       call check( nf90_open(trim(filename),nf90_write,gfile_loc) )
       call check( nf90_inq_varid(gfile_loc,trim(varname),VarId) )

       if(add_saved)then
          allocate( workb2(nlon_regional,nlat_regional))
          allocate( worka2(nlat,nlon))
!!!!!!!! read in guess delp  !!!!!!!!!!!!!!
          call check( nf90_get_var(gfile_loc,VarId,work_b) )
          work_bi(:,:,1)=eta1_ll(nsig+1)
          do i=2,nsig+1
             work_bi(:,:,i)=work_b(:,:,i-1)*0.001_r_kind+work_bi(:,:,i-1)
          enddo
          call fv3_h_to_ll(work_bi(:,:,nsig+1),worka2,nlon_regional,nlat_regional,nlon,nlat)
!!!!!!! analysis_inc Psfc: work_a
          work_a(:,:)=work_a(:,:)-worka2(:,:)
          call fv3_ll_to_h(work_a,workb2,nlon,nlat,nlon_regional,nlat_regional,.true.)
          do k=1,nsig+1
             kr=nsig+2-k
!!!!!!! ges_prsi+hydrostatic analysis_inc !!!!!!!!!!!!!!!!
             work_bi(:,:,k)=work_bi(:,:,k)+eta2_ll(kr)*workb2(:,:)
          enddo
       else
          allocate( workb2(nlon_regional,nlat_regional)) ! added by cTong ! CAPS
          call fv3_ll_to_h(work_a,workb2,nlon,nlat,nlon_regional,nlat_regional,.true.)
          do k=1,nsig+1
             kr=nsig+2-k
!!!!!!! Psfc_ges+hydrostatic analysis_inc !!!!!!!!!!!!!!!!
             work_bi(:,:,k)=eta1_ll(kr)+eta2_ll(kr)*workb2(:,:)
          enddo
       endif
!          delp     
       do k=nsig,1,-1
          kp=k+1
          work_b(:,:,k)=(work_bi(:,:,kp)-work_bi(:,:,k))*1000._r_kind
       enddo

       print *,'write out ',trim(varname),' to ',trim(filename)

       call check( nf90_put_var(gfile_loc,VarId,work_b) )
       call check( nf90_close(gfile_loc) )

       if(add_saved) deallocate(worka2)
       deallocate(workb2)
       deallocate(work_b,work_a,work_bi)

    end if !mype_io

    deallocate(work,work_sub)
end subroutine gsi_fv3ncdf_writeps

subroutine gsi_fv3ncdf_write(filename,varname,var,mype_io,add_saved)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    gsi_nemsio_write
!   pgrmmr: wu
!
! abstract:
!
! program history log:
!
!   input argument list:
!    varu,varv
!    add_saved
!    mype     - mpi task id
!    mype_io
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    use mpimod, only: mpi_rtype,mpi_comm_world,ierror,npe,mype
    use gridmod, only: lat2,lon2,nlon,nlat,lat1,lon1,nsig
    use gridmod, only: ijn,displs_g,itotsub,iglobal
    use mod_fv3_lola, only: fv3_ll_to_h
    use mod_fv3_lola, only: fv3_h_to_ll
    use general_commvars_mod, only: ltosi,ltosj
    use netcdf, only: nf90_open,nf90_close
    use netcdf, only: nf90_write,nf90_inq_varid
    use netcdf, only: nf90_put_var,nf90_get_var
    implicit none

    real(r_kind)   ,intent(in   ) :: var(lat2,lon2,nsig)
    integer(i_kind),intent(in   ) :: mype_io
    logical        ,intent(in   ) :: add_saved
    character(*)   ,intent(in   ) :: varname,filename

    integer(i_kind) :: VarId,gfile_loc
    integer(i_kind) i,j,mm1,k,kr,ns,n,m
    real(r_kind),allocatable,dimension(:):: work
    real(r_kind),allocatable,dimension(:,:,:):: work_sub,work_a
    real(r_kind),allocatable,dimension(:,:,:):: work_b
    real(r_kind),allocatable,dimension(:,:):: workb2,worka2
    real(r_kind) :: clip


    mm1=mype+1

    allocate(    work(max(iglobal,itotsub)*nsig),work_sub(lat1,lon1,nsig))
!!!!!!!! reverse z !!!!!!!!!!!!!!
    do k=1,nsig
       kr=nsig+1-k
       do i=1,lon1
          do j=1,lat1
             work_sub(j,i,kr)=var(j+1,i+1,k)
          end do
       end do
    enddo
    call mpi_gatherv(work_sub,ijnz(mm1),mpi_rtype, &
          work,ijnz,displsz_g,mpi_rtype,mype_io,mpi_comm_world,ierror)

    if(mype==mype_io) then
       allocate( work_a(nlat,nlon,nsig))
       ns=0
       do m=1,npe
          do k=1,nsig
             do n=displs_g(m)+1,displs_g(m)+ijn(m) 
                ns=ns+1
                work_a(ltosi(n),ltosj(n),k)=work(ns)
             end do
          enddo
       enddo

       allocate( work_b(nlon_regional,nlat_regional,nsig))

       call check( nf90_open(trim(filename),nf90_write,gfile_loc) )
       call check( nf90_inq_varid(gfile_loc,trim(varname),VarId) )


       if(add_saved)then
          allocate( workb2(nlon_regional,nlat_regional))
          allocate( worka2(nlat,nlon))
          call check( nf90_get_var(gfile_loc,VarId,work_b) )

          do k=1,nsig
             call fv3_h_to_ll(work_b(:,:,k),worka2,nlon_regional,nlat_regional,nlon,nlat)
!!!!!!!! analysis_inc:  work_a !!!!!!!!!!!!!!!!
             work_a(:,:,k)=work_a(:,:,k)-worka2(:,:)
             call fv3_ll_to_h(work_a(1,1,k),workb2,nlon,nlat,nlon_regional,nlat_regional,.true.)
             work_b(:,:,k)=work_b(:,:,k)+workb2(:,:)
          enddo
          deallocate(worka2,workb2)
       
          ! --- treat negative hydrometers as above zero (Jun from CAPS)
          ! This routine is adopted from 'cliptracers' at 'gridio_wrf.f90' at EnKF
          if ( varname .eq. 'sphum'   .or. varname .eq. 'rainwat' .or. &
               varname .eq. 'snowwat' .or. varname .eq. 'graupel' .or. varname .eq. 'rain_nc'  )  then
             clip = tiny(work_b(1,1,1))
             where ( work_b < clip) work_b = clip
          end if
       else
          do k=1,nsig
             call fv3_ll_to_h(work_a(1,1,k),work_b(1,1,k),nlon,nlat,nlon_regional,nlat_regional,.true.)
          enddo
       endif

       print *,'write out ',trim(varname),' to ',trim(filename)

       call check( nf90_put_var(gfile_loc,VarId,work_b) )
       call check( nf90_close(gfile_loc) )
       deallocate(work_b,work_a)
    end if !mype_io

    deallocate(work,work_sub)

end subroutine gsi_fv3ncdf_write
subroutine check(status)
    use kinds, only: i_kind
    use netcdf, only: nf90_noerr,nf90_strerror
    integer(i_kind), intent ( in) :: status

    if(status /= nf90_noerr) then
       print *,'ncdf error ', trim(nf90_strerror(status))
       stop  
    end if
end subroutine check


end module gsi_rfv3io_mod
