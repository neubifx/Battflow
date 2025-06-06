from pathlib import Path
import shutil


def prepare_topol(doc_id, mols, ans, cats, ions, n_mols_box, md_em_path, md_eq_path, md_prod_path, config):  
    
    # Store content of .itp files and their names
    path_files = [] #full path of files
    stored_files = [] #content of each file
    filesname = []
    for files in md_em_path.iterdir():
            
        if files.suffix == ".itp":
            path_files.append(files)
            filesname.append(files.name)
            
            with open(files, "r") as f:
                lines = f.readlines()
                stored_files.append(lines) # list with readlines lists for each .itp


    #start loop to store copied lines from .itp files

    copied_lines = []
    for _, itp_file in enumerate(stored_files):
        copied_lines.append(f";{filesname[_]}")
        
        for i, line in enumerate(itp_file):
            
            if line.strip() == "[ atomtypes ]":
                           
                for j in range(i + 1, len(itp_file)):
                    copied_line = itp_file[j].strip()
                    
                    copied_lines.append(copied_line)
                    #print(copied_line)
                    if itp_file[j].strip() == "":
                        break #there is only one [ atomtypes ] so break this loop is enough

    #new loop to delete lines after atom types from .itp files
    new_itp = []
    skip = False
    for _, itp_file in enumerate(stored_files):       
        new_itp_lines = []
        
        for i, line in enumerate(itp_file):

            if line.strip() == "[ atomtypes ]":
                skip = True
                
                continue

            if skip:
                
                if line.strip() == "":
                    skip = False
                    
                continue
                
            new_itp_lines.append(line)    
        
        new_itp.append(new_itp_lines)

    #copy new .itp to each file
    for i, file in enumerate(path_files):
    
        for line in new_itp[i]:
        
            with open (file, "w") as f:
                f.writelines(line + "\n" if not line.endswith("\n") else line for line in new_itp[i]) #pretty print
            
    
    
    # add new lines at the end of copied_lines calling the .itp files
    itp_to_top = []
    for name in filesname:
        
        itp_lines = [f";Include {name} topology\n", f'#include "{name}"\n', "\n"]
        copied_lines.extend(itp_lines)
    
    
    #concatenate lists with components names to use now
    component_names  = mols + ans + cats + ions
    
    #define main topol file inside md_path
    main_topol = md_em_path / Path(config["md_simulations"]["topol_main"]).name
    
    with open(main_topol, "r") as f:
        lines = f.readlines()
    
    # build new format of topol file
    new_lines = []
    
    for i, line in enumerate(lines):  
        new_lines.append(line)
        if line.strip() == "[ atomtypes ]":        
            new_lines.extend(copied_lines)
            
        elif line.strip() == "[ system ]":        
            new_lines.append(";doc_id")
            new_lines.append(f"{doc_id}")
    
        elif line.strip() == "[ molecules ]":  
            new_lines.append(";Compounds  nmols")
            
            for i in range(len(component_names)):
                molecules = f"{component_names[i]} {n_mols_box[i]}"
                new_lines.append(molecules)
    
    #copy back to main .topol file
    with open(main_topol, "w") as f:
        f.writelines(line + "\n" if not line.endswith("\n") else line for line in new_lines) #pretty print


    #start copying everything to other folders
    for files in md_em_path.iterdir():
            
        if files.suffix == ".itp" or files.suffix == ".top":
            shutil.copy(files, md_eq_path)
            shutil.copy(files, md_prod_path)