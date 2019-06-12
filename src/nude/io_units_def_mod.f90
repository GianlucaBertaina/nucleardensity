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


      END MODULE io_units_def


      MODULE mpimod
        implicit none
        include 'mpif.h'
        integer my_rank, num_procs
      END MODULE
