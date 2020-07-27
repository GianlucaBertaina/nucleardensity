!    Semiclassical Nuclear Densities. Full notice in LICENSE file. See also README.md
!    Copyright (C) 2020 Chiara Donatella Aieta, Marco Micciarelli, Gianluca Bertaina, Michele Ceotto


!     This module contains all global units for I/O
      MODULE io_units_def
      implicit none
      save

       integer,   parameter   :: unit_input            = 30
       integer,   parameter   :: unit_output           = 31
       integer,   parameter   :: unit_eq_geo           = 32
       integer,   parameter   :: unit_omega            = 33
       integer,   parameter   :: unit_cnorm            = 34
       integer,   parameter   :: unit_wfn              = 35
       integer,   parameter   :: unit_trajMC           = 360
       integer,   parameter   :: unit_cube             = 37
       integer,   parameter   :: unit_geo_exp          = 38
       integer,   parameter   :: unit_geo_exp_h        = 39 
       integer,   parameter   :: unit_bonds_in         = 10000
       integer,   parameter   :: unit_bonds_out        = 10001
       integer,   parameter   :: unit_dihedrals_in     = 20000
       integer,   parameter   :: unit_dihedrals_out    = 20001
       integer,   parameter   :: unit_angles_in        = 30000
       integer,   parameter   :: unit_angles_out       = 30001

      END MODULE io_units_def


      MODULE mpimod
        implicit none
        include "mpif.h"
        integer my_rank, num_procs
      END MODULE
