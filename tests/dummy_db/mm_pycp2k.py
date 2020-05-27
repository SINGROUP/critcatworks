from pycp2k import CP2K
from ase.lattice.cubic import Diamond

#====================== Create the structure with ASE ==========================

#================= Define and setup the calculator object ======================
calc = CP2K()
calc.working_directory = "."
calc.project_name = "cu_mm_bulk"
calc.mpi_n_processes = 2

#==================== Define shortcuts for easy access =========================
CP2K_INPUT = calc.CP2K_INPUT
GLOBAL = CP2K_INPUT.GLOBAL
FORCE_EVAL = CP2K_INPUT.FORCE_EVAL_add()  # Repeatable items have to be first created
MM = FORCE_EVAL.MM
FORCEFIELD = MM.FORCEFIELD
SUBSYS = FORCE_EVAL.SUBSYS
SPLINE = FORCEFIELD.SPLINE_add()
POISSON = MM.POISSON
EWALD = POISSON.EWALD
#======================= Write the simulation input ============================
GLOBAL.Run_type = "ENERGY"
GLOBAL.Print_level = "LOW"
FORCE_EVAL.Method = "FIST"
SPLINE.Emax_spline = 10000
FORCEFIELD.Do_nonbonded = True
CHCu = FORCEFIELD.CHARGE_add()  # Section_parameters can be provided as argument.
CHCu.Atom = "Cu"
CHCu.Charge = 0.0
CHPt = FORCEFIELD.CHARGE_add()  # Section_parameters can be provided as argument.
CHPt.Atom = "Pt"
CHPt.Charge = 0.0
LJCu = FORCEFIELD.NONBONDED.LENNARD_JONES_add()  # Section_parameters can be provided as argument.
# print(LJCu)
LJCu.Atoms = "Cu Cu"
LJCu.Epsilon = "[K_e] {}".format(164.56)
LJCu.Sigma = "[angstrom] {}".format(3.601)
LJCu.Rcut = "[angstrom] {}".format(25.0)

LJPt = FORCEFIELD.NONBONDED.LENNARD_JONES_add()
LJPt.Atoms = "Pt Pt"
LJPt.Epsilon = "[K_e] {}".format(164.56)
LJPt.Sigma = "[angstrom] {}".format(3.601)
LJPt.Rcut = "[angstrom] {}".format(25.0)

SUBSYS.CELL.Abc = "[angstrom] 20 20 20"
SUBSYS.CELL.Periodic = "NONE"

SUBSYS.COORD.Default_keyword = [
    ["Cu", 0.0, 0.0, 0.0],
    ["Cu", 4.0, 0.0, 0.0],
]
SUBSYS.COORD.Unit = "angstrom"

POISSON.Periodic = "None"
EWALD.Ewald_type = "None"


"""
&GLOBAL                  ! section to select the kind of calculation
   RUN_TYPE ENERGY       ! select type of calculation. In this case: ENERGY (=Single point calculation)
&END GLOBAL
&FORCE_EVAL              ! section with parameters and system description
  METHOD FIST            ! Molecular Mechanics method
  &MM                    ! specification of MM parameters
    &FORCEFIELD          ! parameters needed to describe the potential
    &SPLINE
    EMAX_SPLINE 10000    ! numeric parameter to ensure calculation stability. Should not be changed
    &END
        &NONBONDED       ! parameters for the non bonded interactions
          &LENNARD-JONES ! Lennard-Jones parameters
          atoms Cu Cu
          EPSILON    [K_e] 164.56
          SIGMA [angstrom]   3.601
          RCUT  [angstrom]  25.0
        &END LENNARD-JONES
          &LENNARD-JONES ! Lennard-Jones parameters
          atoms Pt Pt
          EPSILON    [K_e] 164.56
          SIGMA [angstrom]   3.601
          RCUT  [angstrom]  25.0
        &END LENNARD-JONES
      &END NONBONDED
      &CHARGE
        ATOM Cu
        CHARGE 0.0
      &END CHARGE
      &CHARGE
        ATOM Pt
        CHARGE 0.0
      &END CHARGE
    &END FORCEFIELD
    &POISSON              ! solver for non periodic calculations
     PERIODIC NONE
      &EWALD
        EWALD_TYPE none
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS                 ! system description
    &CELL
     ABC [angstrom] 20 20 20
     PERIODIC NONE
    &END CELL
    &COORD
      UNIT angstrom
      Cu  0 0 0
      Cu  4 0 0
    &END COORD
   &END SUBSYS
&END FORCE_EVAL
"""
#============ Run the simulation or just write the input file ================
calc.write_input_file()
calc.run()
