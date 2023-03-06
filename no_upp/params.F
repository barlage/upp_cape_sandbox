module params_mod

implicit none

! CONSTANTS.
    real, parameter :: D00=0.0
    real, parameter :: H1M12=1.E-12
    real, parameter :: H99999=99999.
    real, parameter :: H10E5=1.E5
    real, parameter :: CAPA=0.28589641E0
    real, parameter :: CP=1004.6
    real, parameter :: ELIVW=2.72E6
    real, parameter :: ELOCP=ELIVW/CP
    real, parameter :: FMW=18.015
    real, parameter :: FMD=28.964
    real, parameter :: EPS=FMW/FMD
    real, parameter :: ONEPS=1.-EPS
    real, parameter :: G=9.81

end module params_mod
