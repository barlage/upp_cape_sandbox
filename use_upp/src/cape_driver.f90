program cape_driver

use netcdf
use upp_physics, only: calcape
use ctlblk_mod , only: ista_2l,iend_2u,ista,iend,jsta_2l,jend_2u,jsta,jend,im,me,spval,lm
use vrbls2d    , only: ieql, teql
use vrbls3d    , only: t, q, pmid, zint
use masks      , only: lmh
use lookup_mod , only: ptbl , ttbl, rdq, rdth, rdp, rdthe, pl, thl, qs0, sqs, sthe, the0, plq, &
                       ttblq,rdpq,rdtheq,stheq,the0q

implicit none

integer, parameter :: num_atmo_levels = 127
integer, parameter :: cape_type = 2           ! only 2 works for now
integer, parameter :: parcel_level_type2 = 0  ! parcel level relative to surface (0 for surface) for cape_type 2
real   , parameter :: dpbnd_type1 = 10000.0   ! set dpbnd for type 1 cape
integer, parameter :: num_increments = 10     ! n will produce -n:n vector, length 2n+1 (0 for sounding only)
real   , parameter :: t_range = 10.0          ! one sided range of temperature
real   , parameter :: q_range = 5.0           ! one sided rang of spec humid [g/kg]
logical, parameter :: create_netcdf = .true.  ! create a netcdf file, recommended for num_increments > 0

real, dimension(num_atmo_levels)   :: sounding_t, sounding_q, sounding_p, sounding_dz, rtmp
real, dimension(num_atmo_levels+1) :: sounding_z

real                 :: DPBND      ! Depth over which one searches for most unstable parcel.
real,    allocatable :: P1D  (:,:) ! Array of pressure of parcels to lift.
real,    allocatable :: T1D  (:,:) ! Array of temperature of parcels to lift.
real,    allocatable :: Q1D  (:,:) ! Array of specific humidity of parcels to lift.
integer, allocatable :: L1D  (:,:) ! ==not used== Array of model level of parcels to lift.
real,    allocatable :: CAPE (:,:) ! Convective available potential energy (J/kg).
real,    allocatable :: CINS (:,:) ! Convective inhibition (J/kg).
real,    allocatable :: PPARC(:,:) ! Pressure level of parcel lifted when one searches over a particular depth to compute CAPE/CIN.
real,    allocatable :: ZEQL (:,:) ! Height of equilibrium level
real,    allocatable :: THUND(:,:) ! Thunder parameter?

real :: PT

integer :: level, total_bins, irange, jrange, type2_parcel

integer :: ncid, dimid, varid, status   ! netcdf identifiers
integer :: dim_id_i, dim_id_j           ! netcdf dimension identifiers

total_bins = 2*num_increments + 1

lm      = num_atmo_levels
ista_2l = 1
iend_2u = total_bins
ista    = 1
iend    = total_bins
jsta_2l = 1
jend_2u = total_bins
jsta    = 1
jend    = total_bins
im      = huge(1)   ! used?
me      = huge(1)   ! used?
spval   = 9.9e10

allocate(P1D  (ista:iend,jsta:jend))
allocate(T1D  (ista:iend,jsta:jend))
allocate(Q1D  (ista:iend,jsta:jend))
allocate(L1D  (ista:iend,jsta:jend))
allocate(CAPE (ista:iend,jsta:jend))
allocate(CINS (ista:iend,jsta:jend))
allocate(PPARC(ista:iend,jsta:jend))
allocate(ZEQL (ista:iend,jsta:jend))
allocate(THUND(ista:iend,jsta:jend))

allocate(LMH  (ista:iend,jsta:jend)) ! MASKS
allocate(IEQL (ista:iend,jsta:jend)) ! VRBLS2D
allocate(TEQL (ista:iend,jsta:jend)) ! VRBLS2D
allocate(T    (ista:iend,jsta:jend,num_atmo_levels))   ! VRBLS3D
allocate(Q    (ista:iend,jsta:jend,num_atmo_levels))   ! VRBLS3D
allocate(PMID (ista:iend,jsta:jend,num_atmo_levels))   ! VRBLS3D
allocate(ZINT (ista:iend,jsta:jend,num_atmo_levels+1)) ! VRBLS3D

!open(10,file='sounding.txt',form='formatted')
open(10,file='sounding_reformat',form='formatted')

sounding_z(num_atmo_levels+1) = 0.0

do level = num_atmo_levels, 1, -1

!  read(10,*) sounding_p(level), sounding_q(level), sounding_t(level), rtmp(level), rtmp(level), sounding_dz(level)
  read(10,*) sounding_p(level), sounding_q(level), sounding_t(level), sounding_dz(level)
  sounding_z(level) = sounding_z(level+1) - sounding_dz(level)
  
end do
sounding_p = 100.0*sounding_p  ! convert to Pa

type2_parcel = num_atmo_levels - parcel_level_type2

do irange = ista, iend   ! temperature
do jrange = jsta, jend   ! specific humidity

  T    (irange,jrange,:) = sounding_t
  Q    (irange,jrange,:) = sounding_q
  PMID (irange,jrange,:) = sounding_p
  ZINT (irange,jrange,:) = sounding_z

 if(num_increments > 0) then
  T(irange,jrange,type2_parcel) = sounding_t(type2_parcel) + &
           t_range*(irange - num_increments-1)/num_increments
  Q(irange,jrange,type2_parcel) = sounding_q(type2_parcel) + &
           q_range/1000.0*(jrange - num_increments-1)/num_increments
 end if

end do
end do

DPBND = dpbnd_type1
P1D(ista:iend,jsta:jend)   = PMID(ista:iend,jsta:jend,type2_parcel)
T1D(ista:iend,jsta:jend)   = T(ista:iend,jsta:jend,type2_parcel)
Q1D(ista:iend,jsta:jend)   = Q(ista:iend,jsta:jend,type2_parcel)
L1D   = num_atmo_levels
LMH   = num_atmo_levels

PT  = PMID(ista,jsta,1)  ! sounding top pressure
THL = 210.0             ! from INITPOST.f
PLQ = 70000.0           ! from INITPOST.f

call table(PTBL,TTBL,PT,RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)
call tableq(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)

call calcape(cape_type,DPBND,P1D,T1D,Q1D,L1D,CAPE,CINS,PPARC,ZEQL,THUND)

if(create_netcdf) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create the output filename and netcdf file (overwrite old)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  status = nf90_create("cape_output.nc", NF90_CLOBBER, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "idim"   , total_bins , dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "jdim"   , total_bins , dim_id_j)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_def_var(ncid, "cape", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "cape")
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "cin", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "cin")
      if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_def_var(ncid, "temperature", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "temperature at lowest level")
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "specific_humidity", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "specific_humidity at lowest level")
      if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_enddef(ncid)

  status = nf90_inq_varid(ncid, "cape", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , cape)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "cin", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , cins)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "temperature", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , T(ista:iend,jsta:jend,num_atmo_levels))
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "specific_humidity", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , 1000.0*Q(ista:iend,jsta:jend,num_atmo_levels))
    if (status /= nf90_noerr) call handle_err(status)

 status = nf90_close(ncid)
 
else

 print*, 'CAPE: ', cape
 print*, 'CIN:  ', cins

end if

end program

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err
