
WORKFLOW='inversion'    # inversion, migration, modeling
SOLVER='specfem3d'      # specfem2d, specfem3d
SYSTEM='slurm'          # serial, pbs, slurm
OPTIMIZE='base'      # base
PREPROCESS='legacy'    # base
POSTPROCESS='base'   # base

MISFIT='Waveform'
MATERIALS='Acoustic'
DENSITY='Constant'
PRECOND=None


# WORKFLOW
BEGIN=1                 # first iteration
END=30                  # last iteration
NREC=729                # number of receivers
NSRC=5                  # number of sources


# PREPROCESSING
READER='su_specfem3d'   # data file format
CHANNELS='z'            # data channels
NORMALIZE=0             # normalize
BANDPASS=0              # bandpass
MUTE=0                  # mute direct arrival
FREQLO=0.               # low frequency corner
FREQHI=0.               # high frequency corner
MUTECONST=0.            # mute constant
MUTESLOPE=0.            # mute slope


# POSTPROCESSING
SMOOTH=0.               # smoothing radius
SCALE=1.                # scaling factor


# OPTIMIZATION
STEPMAX=10              # maximum trial steps
STEPINIT=0.25           # step length safeguard
STEPFACTOR=0.75


# SOLVER
NT=800                  # number of time steps
DT=0.0025               # time step
F0=0.0                  # dominant frequency


# SYSTEM
NTASK=NSRC              # number of tasks
NPROC=16                # number of processers
WALLTIME=100            # walltime

