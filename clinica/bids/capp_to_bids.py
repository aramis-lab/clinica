"""
Covert the CAPP dataset into the BIDS specification.
Requires fsl.

ToDo:
- Integrate a list that indicates, in case of rescan, what is the best version to consider
- Provide support for the metadata
- Improve the mergeDTI function merging the .bvec and .bval files with python instead that with bash
- Insert a timeline inside the warning file
- Fix the code style

@author: Sabrina Fontanella
"""

from os import path
from glob import glob
import os
from shutil import copy
import logging
import nibabel as nib
from missing_mods_tracker import MissingModsTracker

def remove_rescan(list_path):
    """
    Remove all the folders containing the keyword 'rescan' from
    a given list of folders.

    Args:
        list_path (str): the list of the files to analize.

    Returns:
        The list of non rescanned folders.

    """

    #Extract all the folders without the substring 'rescan'
    noResc_lst = [s for s in list_path if 'rescan' not in s]
    if len(noResc_lst)!=len(list_path):
        for r_file in list(set(list_path)-set(noResc_lst)):
            logging.warning('Rescan found '+r_file+' Ignored.')
    return noResc_lst


def merge_DTI(folder_input, folder_output, name):
    """
    Merge all the DTI files of a given subject.

    For the merging only DTI folders containing all .nii.gz, .bval and .bvec are considered,
    otherwise the folder is ignored.

    Args:
        folder_input (str) : the folder containing all the DTI to merge.
        folder_output (str) : the folder where store the merged file.
        name (str): name to give to the merged file.
    """
    img = []
    bval = []
    bvec = []
    dti_list = remove_rescan(glob(path.join(folder_input, '*DTI*')))
    incomp_folders = []
    nr_dti = len(dti_list)
    if nr_dti ==0:
        print 'No DTI found for '+folder_input
        return -1
    else:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)
        for folder in dti_list:
            if len(glob(path.join(folder,'*.bval'))) !=0 and len(glob(path.join(folder,'*.bvec'))) !=0:
                img.append(glob(path.join(folder,'*.nii*'))[0])
                bval.append(glob(path.join(folder,'*.bval'))[0])
                bvec.append(glob(path.join(folder,'*.bvec'))[0])
            else:
                print '.bvec or .bval not found for DTI folder '+folder+ ' Skipped'
                incomp_folders.append(folder)

        # if it has been found at least a DTI folder complete with bvec, bval and nii.gz
        if len(img)>0:
            #merge all the .nii.gz file with fslmerge
            os.system('fslmerge -t '+path.join(folder_output,name+'.nii.gz')+' '+" ".join(img))
            #merge all the .bval files
            os.system('paste '+" ".join(bval) +' > '+ path.join(folder_output,name+'.bval'))
            #merge all the .bvec files
            os.system('paste ' + " ".join(bvec) + ' > ' + path.join(folder_output, name + '.bvec'))

        if len(incomp_folders)>0:
            return incomp_folders
        else:
            return None


def choose_correction(dir, to_consider, mod):
    """
    Decides what is the best file type to choose depending on a priority list.

    The selection criteria is based on the input list toConsider that contains
    the files to consider in order of priority.

    Args:
        dir: the directory containing the files.
        to_consider: the list of files to consider in order of priority.
        mod: the modality

    Returns:
        The best file (according to the list to_consider) for the given modality or:
        -1 if if the modality is missing.
        0 if none of the desided corrections is available..

    """
    # extract all the files availabe for a certain modality
    correction_list = remove_rescan(glob(path.join(dir, '*' + mod + '*')))
    if len(correction_list)==0:
        return -1
    if len(correction_list) == 1:
        return correction_list[0].split(os.sep)[-1]
    else:
        for i in range(0, len(to_consider)):
            if any(to_consider[i] in c for c in correction_list):
                return to_consider[i]
        return 0


def convert(source_dir, dest_dir):
    """
    :param source_dir: directory of the input dataset
    :param dest_dir: output directory
    """
    capp_spath = []
    capp_ids = []
    bids_ids = []
    t1_selected = 'none'
    t1_to_consider = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    sessions = ["M00", "M18"]
    anat = ['T1', 'T2']
    anat_modalities = ['T1', 'T2', 'FLAIR']
    # mod_to_ignore = ['other', 'MAP', 'Survey', 'GRE', 'DTI', 'TSE']
    mod_to_consider = ['T1', 'T2']
    map_problems = ["13001PBA20150623M18B0MAPph_S016.nii.gz", "13002PRJ20150922M18B0MAPph_S014.nii.gz","11001PGM20130704M00B0MAPph_S010.bval","07002PPP20150116M18B0MAPph_S009.nii.gz", "07003PGM20141217M18B0MAPph_S010.nii.gz" ]
    bids_suff = {
        'T1': '_T1w',
        'T2': '_T2w',
        'Flair': '_FLAIR',
        'MapPh': '_phasediff',
        'Map'  : '_magnitude',
        'fMRI' : '_bold',
        'dwi'  : '_dwi'
    }

    mmt = MissingModsTracker()
    os.mkdir(dest_dir)
    readme = open(path.join(dest_dir,'README.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'logfile.log'), format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %I:%M' )
    participants = open(path.join(dest_dir,'participants.tsv'), 'w')
    cities_folder = glob(path.join(source_dir,'*','*','*'))

    #create the lists of bids_ids extracting the list of the subjects from the dataset
    for cf in cities_folder:
        for sub_path in glob(path.join(cf, "*")):
            sub_id = sub_path.split(os.sep)[-1]
            # the substring 'xxx' inside a folder subject name indicates an acquisition to ignore
            if 'xxx' in sub_id:
                logging.warning('Anomalous subject folder: '+sub_id+'. Ignored.\n')
            else:
                capp_spath.append(sub_path)
                capp_ids.append(sub_id)
                bids_ids.append('sub-CAPP'+sub_id)
                os.mkdir(path.join(dest_dir, bids_ids[-1]))

    participants.write("------------------------------\n")
    participants.write("participant_id   BIDS_id\n")
    participants.write("------------------------------\n")
    for c_ids, b_ids in zip(capp_ids, bids_ids):
        participants.writelines(c_ids+"     "+b_ids+'\n')
    print "*******************************"
    print "CAPP to BIDS converter"
    print "*******************************"
    print "Number of subjects: ", len(capp_spath)
    #for each subject extract the list of files and convert them into BIDS specification
    for sub_path in capp_spath:
        print "Converting:", sub_path
        logging.info("Converting:"+sub_path)
        #for each subsession
        for ses in sessions:
            if os.path.exists(path.join(sub_path,ses)):
                #extracting the index of the subject
                subj_bindex = capp_spath.index(sub_path)
                os.mkdir(path.join(dest_dir, bids_ids[subj_bindex], 'ses-'+ses))
                ses_dir_bids = path.join(dest_dir, bids_ids[subj_bindex], 'ses-'+ses)
                bids_file_name = bids_ids[subj_bindex]+'_ses-'+ses
                session_path = path.join(sub_path,ses)
                type_folders = glob(path.join(sub_path,ses,"NIFTI","*"))
                curr_dir = path.join(sub_path,ses, 'NIFTI')

                # convert the fieldmap data
                map = remove_rescan(glob(path.join(sub_path, ses, "NIFTI", "*MAP_*", '*.nii.gz')))
                mapPh = remove_rescan(glob(path.join(sub_path, ses, "NIFTI", "*MAPph_*", '*.nii.gz')))
                if len(map)==0:
                    logging.warning('Missing magnitude image(s) for '+curr_dir)
                    mmt.add_missing_mod('Map')
                if len(mapPh)==0:
                    logging.warning('Missing phase image(s) for '+curr_dir)
                    mmt.add_missing_mod('MapPh')
                #If the information regarding the Fieldmap data are complete
                if len(map)>0 and len(mapPh)>0:
                    map_ph_name = mapPh[0].split(os.sep)[-1]
                    map_name = map[0].split(os.sep)[-1]
                    #toSolve: there are some files that produce an error when loaded with Nibabel!
                    if (map_ph_name not in map_problems) and (map_name not in map_problems):
                        #open the files with Nibabel
                        map_nib = nib.load(map[0])
                        mapPh_nib = nib.load(mapPh[0])
                        dim_map = (map_nib.header['dim'])[4]
                        dim_mapPh = (mapPh_nib.header['dim'])[4]
                        os.mkdir(path.join(ses_dir_bids, "fmap"))
                        #Case 1: one phase difference image and at least one magnitude image
                        if dim_mapPh ==1 and dim_map >0:
                            copy(mapPh[0], path.join(ses_dir_bids, "fmap", bids_file_name +bids_suff['MapPh']+'.nii.gz'))
                            os.system('fslsplit ' + map[0]+ ' '+path.join(ses_dir_bids, "fmap/"+bids_file_name+bids_suff['Map']))
                        #Case 2: two phase images and two magnitude images'
                        elif dim_mapPh==2 and dim_map==2:
                            os.system('fslsplit ' + mapPh[0] + ' ' + path.join(ses_dir_bids, "fmap/" + bids_file_name +bids_suff['MapPh']))
                            os.system('fslsplit ' + map[0] + ' ' + path.join(ses_dir_bids,"fmap/" + bids_file_name +bids_suff['Map']))

                #decide which T1 will be considered for the conversion
                t1_selected = choose_correction(curr_dir, t1_to_consider, 'T1')
                if t1_selected!=-1 and t1_selected!=0:
                    os.mkdir(path.join(ses_dir_bids, "anat"))
                    data_files = glob(path.join(curr_dir, '*' + t1_selected + '*', '*.nii*'))[0]
                    copy(data_files, path.join(ses_dir_bids,'anat',bids_file_name+bids_suff['T1']+'.nii.gz'))
                else:
                    if t1_selected == -1:
                        logging.info('No T1 found for ' + curr_dir)
                        mmt.add_missing_mod('T1')
                    else:
                        logging.warning('None of the desiderd T1 corrections is available for ' +   curr_dir)



                #extract and convert the T2 FLAIR modality if is available
                t2Fl_capp_lst = remove_rescan(glob(path.join(curr_dir, '*' + '*T2FLAIR*')))
                if len(t2Fl_capp_lst)>0:
                    if not os.path.exists(path.join(path.join(ses_dir_bids, "anat"))):
                        os.mkdir(path.join(path.join(ses_dir_bids, "anat")))
                    t2Fl_capp_file = glob(path.join(t2Fl_capp_lst[0],'*.nii.gz*'))
                    copy(t2Fl_capp_file[0], path.join(ses_dir_bids, 'anat', bids_file_name + bids_suff['Flair'] + '.nii.gz'))
                else:
                    logging.info('No FLAIR found for '+curr_dir)
                    mmt.add_missing_mod('FLAIR')

                #combine and convert all valid DTI folders for the same subject
                out = merge_DTI(curr_dir, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                if out != None: #missing or incomplete DTI
                    if out == -1:
                        logging.info('No DTI found for ' + curr_dir)
                        mmt.add_missing_mod('DTI')
                    else:
                        for e in out:
                            logging.warning('.bvec or .bval not found for DTI folder ' + e + ' Skipped')


                #extract and convert the fMRI
                fmriC_lst = remove_rescan(glob(path.join(curr_dir, '*' + '*fMRI*')))
                if len(fmriC_lst) > 0:
                    os.mkdir(path.join(path.join(ses_dir_bids, "func")))
                    fmriC_file = glob(path.join(fmriC_lst[0], '*.nii*'))[0]
                    copy(fmriC_file, path.join(ses_dir_bids, 'func', bids_file_name +'_task-rest'+bids_suff['fMRI']+'.nii.gz'))
                else:
                    logging.info('Non fMRI found for '+curr_dir)
                    mmt.add_missing_mod('fMRI')

        logging.info("Conversion for the subject terminated.\n")

    #printing the statistics about missing modalities into a file
    readme.write('Number of subjects converted: '+str(len(capp_spath))+'\n')
    readme.write('********************************\n')
    readme.write('Number of missing files per modalities (more details inside the logfile):\n')
    mmt.remove_mods_unused()
    missing_list = mmt.get_missing_list()
    for key in missing_list:
        readme.write(key+':'+str(missing_list[key])+"\n")
    readme.write('********************************\n')
    readme.close()