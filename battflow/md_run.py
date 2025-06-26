import gromacs
import gromacs.run

import os
import shutil

from pathlib import Path


def md_simulation_run(BASE_DIR, config, work_path, md_em_path, md_eq_path, md_prod_path):
    
    class MDrunnerMPI(gromacs.run.MDrunner):
        """Custom MDrunner that uses mpiexec to launch gmx_mpi mdrun."""
        
        # mdrun binary (can be a string or a tuple of candidates)
        mdrun = config["md_run_env"]["mdrun"]
        
        mpiexec = config["md_run_env"]["mpiexec"]
    
    #execution of the MD simulation workflow
    print(md_em_path)
    print(md_eq_path)
    print(md_prod_path)
    
    #enter the correct folder
    os.chdir(md_em_path)

    #prepare files for energy minimisation
    shutil.copy(BASE_DIR / "resources" / "md_run_files" / "minim.mdp", md_em_path)
    gromacs.grompp(f= "minim.mdp", c="electrolyte.gro", p="topol_gaff2.top", o="em.tpr", maxwarn = 99)

    #run energy minimisation
    mdrun_mpi = MDrunnerMPI(v=True, deffnm="em")
    rc = mdrun_mpi.run(ncores=config["md_run_env"]["ncores"])

    #enter the correct folder
    os.chdir(md_eq_path)

    #prepare files for calibration run
    shutil.copy(BASE_DIR / "resources" / "md_run_files" / "eq.mdp", md_eq_path)
    shutil.copy(md_em_path / "em.gro", md_eq_path)

    gromacs.grompp(f= "eq.mdp", c="em.gro", r="em.gro", p="topol_gaff2.top", o="eq.tpr", maxwarn = 99)

    #run calibration
    mdrun_mpi = MDrunnerMPI(v=True, deffnm="eq")
    rc = mdrun_mpi.run(ncores=config["md_run_env"]["ncores"])

    #enter the correct folder
    os.chdir(md_prod_path)

    #prepare files for production run
    shutil.copy(BASE_DIR / "resources" / "md_run_files" / "prod.mdp", md_prod_path)
    shutil.copy(md_eq_path / "eq.gro", md_prod_path)
    shutil.copy(md_eq_path / "eq.cpt", md_prod_path)

    gromacs.grompp(f= "prod.mdp", c="eq.gro", t="eq.cpt", p="topol_gaff2.top", o="prod_0_1.tpr", maxwarn = 99)

    # run the production run
    mdrun_mpi = MDrunnerMPI(v=True, deffnm="prod_0_1")

    rc = mdrun_mpi.run(ncores=config["md_run_env"]["ncores"])   

    #Safely return to the work folder
    os.chdir(work_path)
