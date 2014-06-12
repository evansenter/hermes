BINDIR            = ~/bin
LIBDIR            = ~/lib
GCC               = gcc
GXX               = g++
LIB               = ../../lib
SHARED_HEADER     = ../../h
HEADER            = .
CODE              = .
GCC_VERSION       = $(shell expr `$(GCC) -dumpversion`)
GCC_MAJ_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
GCC_MIN_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
GCC_PATCH_VERSION = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION   = $(shell expr $(GCC_MAJ_VERSION) \+ $(GCC_MIN_VERSION) \+ $(GCC_PATCH_VERSION))
GCC_GTEQ_4.6.0    = $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0    = $(shell expr $(GCC_NUM_VERSION) \>= 40900)
GXX_VERSION       = $(shell expr `$(GXX) -dumpversion`)
GXX_MAJ_VERSION   = $(shell expr `echo $(GXX_VERSION) | cut -d . -f 1` \* 10000)
GXX_MIN_VERSION   = $(shell expr `echo $(GXX_VERSION) | cut -d . -f 2` \* 100)
GXX_PATCH_VERSION = $(shell expr `echo $(GXX_VERSION) | cut -d . -f 3`)
GXX_NUM_VERSION   = $(shell expr $(GXX_MAJ_VERSION) \+ $(GXX_MIN_VERSION) \+ $(GXX_PATCH_VERSION))
GXX_GTEQ_4.6.0    = $(shell expr $(GXX_NUM_VERSION) \>= 40600)
GXX_GTEQ_4.9.0    = $(shell expr $(GXX_NUM_VERSION) \>= 40900)
