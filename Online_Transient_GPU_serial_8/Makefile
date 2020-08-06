########

#Compiler and Linker
CC          := g++

#The Target Binary Program
TARGET      := trans

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := source
INCDIR	    := include
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := source
OUTDIR      := output

#extensions, no need to modify
SRCEXT      := cpp
SRCEXTCU      := cu
DEPEXT      := d
OBJEXT      := o
CUDA_ROOT = /usr/local/cuda

#ARCH = sm_30

NVCC = $(CUDA_ROOT)/bin/nvcc


#NVCCFLAGS = -O3 $(INCLUDES) -arch $(ARCH)

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

#$(NVCC) -c source/trans_cuda.cu -o obj/trans_cuda,o

#LIBINCLUDE += -I/Users/peng/graphlu_mac/graphlu_test/util
#LIBINCLUDE += -I/Users/peng/graphlu_mac/graphlu_test/include
#LIBINCLUDE += -I$(shell pwd)/$(INCDIR)

#LIBINCLUDE += -I/home/firehawk/Dropbox/summer_intern/transient_analysis/for_Xiang/graphlu_code/graphlu_test/util
#LIBINCLUDE += -I/home/firehawk/Dropbox/summer_intern/transient_analysis/for_Xiang/graphlu_code/graphlu_test/include
#LIBINCLUDE += -I/home/firehawk/boost_1_53_0
LIBINCLUDE += -I/home/bowen/zantai/Transient/graphlu_code/graphlu_test/util
LIBINCLUDE += -I/home/bowen/zantai/Transient/graphlu_code/graphlu_test/include
#LIBINCLUDE += -I/home/bowen/zantai/Transient/boost_1_73_0
LIBINCLUDE += -I/usr/local/boost
LIBINCLUDE += -I/usr/local/cuda/include
#LIBINCLUDE += -I/home/xiangzhang/graphlu_code/graphlu_test/util
#LIBINCLUDE += -I/home/xiangzhang/graphlu_code/graphlu_test/include
#LIBINCLUDE += -I/home/xiangzhang/boost_1_53_0
LIBINCLUDE += -I$(shell pwd)/$(INCDIR)

#LDFLAGS += -L/Users/peng/graphlu_mac/graphlu_test/util/
#LDFLAGS += -L/home/firehawk/Dropbox/summer_intern/transient_analysis/for_Xiang/graphlu_code/graphlu_test/util/
#LDFLAGS += -L/home/firehawk/boost_1_53_0
#LDFLAGS += -L/home/bowen/zantai/Transient/graphlu_code/graphlu_test/util/
#LDFLAGS += -L/home/bowen/zantai/Transient/boost_1_73_0
LDFLAGS += -L/usr/local/boost
LDFLAGS += -L/home/xiangzhang/graphlu_code/graphlu_test/util/
#LDFLAGS += -L/home/xiangzhang/boost_1_53_0
LDFLAGS += -lgraphlu_util -lm -pthread
#LDFLAGS += -np 4
#LDFLAGS += -L/Users/peng/graphlu_mac/graphlu_test/lib/
#LDFLAGS += -L/home/firehawk/Dropbox/summer_intern/transient_analysis/for_Xiang/graphlu_code/graphlu_test/lib/
LDFLAGS += -L/home/bowen/zantai/Transient/graphlu_code/graphlu_test/lib/
LDFLAGS += -L/usr/local/cuda/lib64/
#LDFLAGS += -L/home/xiangzhang/graphlu_code/graphlu_test/lib/
LDFLAGS += -lgraphlu -lcudart

#Flags, Libraries and Includes
CFLAGS      := -Wall -pedantic -std=c++11 -pthread
NVCCCFLAGS      := -arch=sm_35 -std=c++11 
LIB         := $(LDFLAGS)
INC         := $(LIBINCLUDE) -I/usr/local/include
INCDEP      := $(LIBINCLUDE)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#SOURCESCU     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXTCU))
#OBJECTSCU     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCESCU:.$(SRCEXTCU)=.$(OBJEXT)))

#Default Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(OUTDIR)

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR) $(TARGETDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR) $(BUILDDIR) $(OUTDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) obj/transient.o obj/ode_solver.o obj/generatorModels.o obj/apply_perturbation.o $^ $(LIB)
#Compile
#$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXTCU)
#	@mkdir -p $(dir $@)
#	$(NVCC) $(NVCCCFLAGS)  $(INC) -c -o $@ $<
#	@$(NVCC) $(NVCCCFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXTCU) > $(BUILDDIR)/$*.$(DEPEXT)
#	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
#	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
#	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
#	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp


#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(NVCC) $(INC) $(NVCCCFLAGS) -c source/transient.cu -o obj/transient.o
	#$(NVCC) $(INC) $(NVCCCFLAGS) -c source/trans_cuda.cu -o obj/trans_cuda.o
	$(NVCC) $(INC) $(NVCCCFLAGS) -c source/ode_solver.cu -o obj/ode_solver.o
	$(NVCC) $(INC) $(NVCCCFLAGS) -c source/generatorModels.cu -o obj/generatorModels.o
	$(NVCC) $(INC) $(NVCCCFLAGS) -c source/apply_perturbation.cu -o obj/apply_perturbation.o
	$(CC) $(CFLAGS)  $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all remake clean cleaner resources
