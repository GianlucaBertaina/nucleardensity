!     This module contains mathematical and physical constants
      MODULE constants
      implicit none
      save
!     conversion from a.u. [codice in a.u.]
       real(8)   , parameter  :: FROMauTOev      = 27.2114d0
       real(8)   , parameter  :: FROMevTOau      = 1.d0/FROMauTOev
       real(8)   , parameter  :: FROMauTOcm      = 219474.625d0
       real(8)   , parameter  :: FROMcmTOau      = 1.d0/FROMauTOcm
       real(8)   , parameter  :: FROMauTOang     = 0.5291772108d0
       real(8)   , parameter  :: FROMangTOau     = 1/0.5291772108d0
       real(8)   , parameter  :: FROMangTOau_vel = 1.d0/21.87676d0
       real(8)   , parameter  :: pi              = dacos(-1.d0)
       real(8)   , parameter  :: twopi           = 2.d0 * pi
       real(8)   , parameter  :: zero_real       = 0.d0
       complex(8), parameter  :: aye             = (0.d0,1.d0)


! DA AGGIUSTARE PERCHE' POCO CHIARE: (DA CHIEDERE)
       real(8)   , parameter  :: toautime   = 4.1341d4
       real(8)   , parameter  :: toauhess   = 1.066D-6
       real(8)   , parameter  :: toaueng    = 3.8088d-4



! old constants: 
!      toev = 27.2114d0
!      tocm = 219474.625d0
!      toau = 1.d0/0.529177d0
!      toauvel = 1.d0/21.87676d0
!      toaueng = 3.8088d-4




      END MODULE constants
