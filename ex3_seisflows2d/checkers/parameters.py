
WORKFLOW='inversion'    # inversion, migration
SOLVER='specfem2d'      # specfem2d, specfem3d
SYSTEM='serial'         # serial, pbs, slurm
OPTIMIZE='base'         # base, newton
PREPROCESS='legacy'     # base
POSTPROCESS='base'      # base

MISFIT='Waveform'
MATERIALS='LegacyAcoustic'
DENSITY='Constant'


# WORKFLOW
BEGIN=1                 # first iteration
END=5                   # last iteration
NREC=132                # number of receivers
NSRC=25                 # number of sources
SAVEGRADIENT=1          # save gradient how often


# PREPROCESSING
READER='su_specfem2d'   # data file format
CHANNELS='y'            # data channels
NORMALIZE=0             # normalize
NORMALIZE_ALL=0         # normalize
BANDPASS=0              # bandpass
MUTE=0                  # mute direct arrival
FREQLO=0.               # low frequency corner
FREQHI=0.               # high frequency corner
MUTECONST=0.            # mute constant
MUTESLOPE=0.            # mute slope


# POSTPROCESSING
SMOOTH=20.              # smoothing radius
SCALE=6.0e6             # scaling factor


# OPTIMIZATION
PRECOND=None            # preconditioner type
STEPMAX=10              # maximum trial steps
STEPTHRESH=0.1          # step length safeguard


# SOLVER
NT=4800                 # number of time steps
DT=0.06                 # time step
F0=0.084                # dominant frequency


# SYSTEM
NTASK=1                 # must satisfy 1 <= NTASK <= NSRC
NPROC=1                 # processors per task

