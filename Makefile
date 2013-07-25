ifndef DES_HOME
DES_HOME=${HOME}/DES_HOME
endif

build:
	@echo "PSM: Ready to install to ${DES_HOME}"

install: 
	@echo "PSM: Installing to ${DES_HOME}"
	-mkdir -p ${DES_HOME}
	-rsync -vCaq bin/setup.bash ${DES_HOME}/bin/PSM_setup.bash
	-rsync -vCaq bin/setup.csh ${DES_HOME}/bin/PSM_setup.csh
	-rsync -vCaq bin/runPSM ${DES_HOME}/bin/runPSM
	-mkdir -p ${DES_HOME}/share/java
	-rsync -vCaq lib/* ${DES_HOME}/share/java/PSM

installunitrunning:
	@echo "PSM: Installing to ${INSTALL_DIR}"
	-mkdir -p ${INSTALL_DIR}/pyPSM
	-rsync -vCaq pyPSM/python ${INSTALL_DIR}/pyPSM
	

clean:
	@echo "PSM: Nothing to clean"
