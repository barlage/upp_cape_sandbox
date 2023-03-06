program cape_driver

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

real, dimension(num_atmo_levels)   :: sounding_t, sounding_q, sounding_p, sounding_dz, rtmp
real, dimension(num_atmo_levels+1) :: sounding_z

real                 :: DPBND      ! Depth over which one searches for most unstable parcel.
real,    allocatable :: P1D  (:,:) ! Array of pressure of parcels to lift.
real,    allocatable :: T1D  (:,:) ! Array of temperature of parcels to lift.
real,    allocatable :: Q1D  (:,:) ! Array of specific humidity of parcels to lift.
integer, allocatable :: L1D  (:,:) ! (not used) Array of model level of parcels to lift.
real,    allocatable :: CAPE (:,:) ! Convective available potential energy (J/kg).
real,    allocatable :: CINS (:,:) ! Convective inhibition (J/kg).
real,    allocatable :: PPARC(:,:) ! Pressure level of parcel lifted when one searches over a particular depth to compute CAPE/CIN.
real,    allocatable :: ZEQL (:,:) ! Height of equilibrium level
real,    allocatable :: THUND(:,:) ! Thurder parameter?

real :: PT

integer :: level

lm      = num_atmo_levels
ista_2l = 1
iend_2u = 1
ista    = 1
iend    = 1
jsta_2l = 1
jend_2u = 1
jsta    = 1
jend    = 1
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

open(10,file='sounding_reformat',form='formatted')

sounding_z(num_atmo_levels+1) = 0.0

do level = num_atmo_levels, 1, -1

  read(10,*) sounding_p(level), sounding_q(level), sounding_t(level), sounding_dz(level)
  sounding_z(level) = sounding_z(level+1) - sounding_dz(level)
  
end do
sounding_p = 100.0*sounding_p  ! convert to Pa

T    (ista,jsta,:) = sounding_t
Q    (ista,jsta,:) = sounding_q
PMID (ista,jsta,:) = sounding_p
ZINT (ista,jsta,:) = sounding_z

DPBND = huge(1.0)
P1D   = sounding_p(num_atmo_levels)
T1D   = sounding_t(num_atmo_levels)
Q1D   = sounding_q(num_atmo_levels)
L1D   = num_atmo_levels
LMH   = num_atmo_levels

PT = PMID(ista,jsta,1)  ! sounding top pressure
THL = 210.0             ! from INITPOST.f
PLQ = 70000.0           ! from INITPOST.f

call table(PTBL,TTBL,PT,RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)
call tableq(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)

call calcape(cape_type,DPBND,P1D,T1D,Q1D,L1D,CAPE,CINS,PPARC,ZEQL,THUND)

print*, "cape: ", cape
print*, "cin:  ", cins

end
