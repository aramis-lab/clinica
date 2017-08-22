#participant_tsv='/Volumes/aramis-projects/CLINICA/CLINICA_datasets/BIDS/AIBL_BIDS/participants.tsv'

#path_to_bids='/Volumes/aramis-projects/CLINICA/CLINICA_datasets/BIDS/AIBL_BIDS/'

def find_mci_converted(participant_tsv,path_to_bids):
    #this method finds the mci converted.
    #participant_tsv=tsv file in BIDS directory"
    participants = pd.io.parsers.read_csv(os.path.join(participant_tsv), sep='\t')
    rid=participants.alternative_id_1
    path=participants.participant_id
    mci=[] #total of MCI
    mci_converted_m18=[] #mci who converted at m18
    mci_converted_m36=[]#mci who converted at m36
    mci_converted_m54=[]#mci who converted at 54
    mci_only_bl=[] #patients where we have only the M00 visit so we cannot say if they converted
    mci_non_converted=[]
    for i in path:
        #iterates for each patient
        session_file=os.path.join(path_to_bids,i,i+ '_sessions.tsv')
        if os.path.exists(session_file):
            session_file_read=pd.io.parsers.read_csv(session_file, sep='\t')
            index = session_file_read.session_id[session_file_read.session_id == 'M00'].index.tolist()
            #read the diagnosis in the sessions file
            diagnosis = session_file_read.diagnosis[index]
            if diagnosis[0]=='MCI':
                mci.append(i)
                if len(session_file_read)>1:
                    for j in xrange(len(session_file_read)):
                        if 'M18' in session_file_read['session_id'][j]:
                            indexm18=session_file_read.session_id[session_file_read.session_id == 'M18'].index.tolist()
                            diagnosism18 = session_file_read.diagnosis[indexm18]
                            if diagnosism18[1] == 'AD':
                                mci_converted_m18.append(i)
                        if 'M36' in session_file_read['session_id'][j] and i not in mci_converted_m18:
                            indexm36 = session_file_read.session_id[session_file_read.session_id == 'M36'].index.tolist()
                            diagnosism36 = session_file_read.diagnosis[indexm36]
                            if diagnosism36[indexm36[0]] == 'AD':
                                mci_converted_m36.append(i)
                        if 'M54' in session_file_read['session_id'][j] and i not in mci_converted_m18 and i not in mci_converted_m36:
                            indexm54 = session_file_read.session_id[session_file_read.session_id == 'M54'].index.tolist()
                            diagnosism54 = session_file_read.diagnosis[indexm54]
                            if diagnosism54[indexm54[0]] == 'AD':
                                mci_converted_m54.append(i)
                if len(session_file_read)==1:
                    mci_only_bl.append(i) #patient where we have only the baseline so we cannot say if they convert or not
    mci_non_converted=len(mci_converted_m54)+len(mci_converted_m36)+len(mci_only_bl) - len(mci)
    print len(mci)
    print len(mci_converted_m18)
    print len(mci_converted_m36)
    print len(mci_converted_m54)
    print len(mci_only_bl)
    print len(mci_non_converted)

#aibl_dartel='/Users/simona.bottani/AIBL/aibl_dartel.tsv' PATIENTS WHERE

def find_parameters_statistics(participant_tsv,path_to_bids,aibl_dartel,output_table):
    #output_table=where to save the final tsv
    #this function create a tsv file where for each participant of aibl_dartel it's reported the age, diagnosis and mmscore at the baseline

    import re
    import pandas
    participants = pd.io.parsers.read_csv(os.path.join(participant_tsv), sep='\t')
    dartel=pd.io.parsers.read_csv(os.path.join(aibl_dartel), sep='\t')
    rid_dartel=dartel.participant_id
    rid=participants.alternative_id_1
    path=participants.participant_id
    #sub_without_exam_date=[]
    dict=[]
    all_gender=[]
    all_diagnosis=[]
    all_date=[]
    all_mmscore=[]
    #for i in path:
    for i in rid_dartel:
        #for each patient used in the dartel for the paper
        session_file=os.path.join(path_to_bids,i,i+ '_sessions.tsv')
        if os.path.exists(session_file):
            session_file_read=pd.io.parsers.read_csv(session_file, sep='\t')
        if i!='sub-AIBL1503':
            date_bl=session_file_read.examination_date[0]
            #if date_bl!=-4:
            year_bl=re.search('[0-9].*/[0-9].*/([0-9].*)', str(date_bl)).group(1)  # string from image directory
        else:
            year_bl=2014

        dob = participants.loc[(participants["participant_id"] == i), 'date_of_birth']
        gender=participants.loc[(participants["participant_id"] == i), 'sex']
        index = session_file_read.session_id[session_file_read.session_id == 'M00'].index.tolist()
        diagnosis = session_file_read.diagnosis[index][0]
        mmscore=session_file_read.MMS[index][0]
        all_mmscore.append(mmscore) #mmscore at M00
        all_diagnosis.append(diagnosis) #diagnosis at M00
        age=int(int(year_bl)-dob) #age at M00
        all_date.append(age)
        all_gender.append(gender) #gender


    all_sex=[]
    for j in xrange(len(all_gender)):
        #sex=all_gender[j][j]
        sex=all_gender[j]
        all_sex.append(sex)

    dict = pandas.DataFrame({'subjects': rid_dartel,
                             'sex': all_sex,
                             'age': all_date,
                             'diagnosis': all_diagnosis,
                             'mmscore': all_mmscore
                             })
    dict.to_csv(os.path.join(output_table), sep='\t', index=False, encoding='utf8')
    #where to save the table#




import pandas as pd
import numpy as np

# AIBL
table_tsv = '/Users/simona.bottani/AIBL/3-table.tsv'
def statistics(table_tsv):
    #all the computaton to create the table for the paper
    participants = pd.io.parsers.read_csv(table_tsv, sep='\t')
    participant_id = list(participants.subjects)
    sex = np.asarray(list(participants.sex))
    age = np.asarray(list(participants.age))
    mmse = np.asarray(list(participants.mmscore))
    diagnosis = np.asarray(list(participants.diagnosis))

    for i in range(len(sex)):
        print (str(participant_id[i]) + ' ' + str(sex[i]) + ' ' + str(age[i]) + ' ' + str(mmse[i]))

# General statistics

    N_AD = sum(diagnosis == 'AD')
    mask_non_AD = diagnosis != 'AD'

    N_CN = sum(diagnosis == 'CN')
    mask_non_CN = diagnosis != 'CN'

    N_MCI = sum(diagnosis == 'MCI')
    mask_non_MCI = diagnosis != 'MCI'


# Number of subjects
    print('N AD = ' + str(N_AD) + ', N CN = ' + str(N_CN)+', N MCI = ' + str(N_MCI))
    if N_AD + N_CN  + N_MCI!= len(sex):
        raise(Exception('AD + CN + MCI != total number of subjects'))

# Age
    print('Age AD = ' + str(np.mean(np.ma.masked_array(age, mask=mask_non_AD))) + ' +/- ' + str(np.std(np.ma.masked_array(age, mask_non_AD))))
    print('Age CN = ' + str(np.mean(np.ma.masked_array(age, mask=mask_non_CN))) + ' +/- ' + str(np.std(np.ma.masked_array(age, mask_non_CN))))
    print('Age MCI = ' + str(np.mean(np.ma.masked_array(age, mask=mask_non_MCI))) + ' +/- ' + str(np.std(np.ma.masked_array(age, mask_non_MCI))))
    print('Mean age = ' + str(np.mean(age)) + ' +/- ' + str(np.std(age)))

# Male / Female
    female_and_AD = np.array(diagnosis == 'AD').astype(int) * np.array(sex == 'F').astype(int)
    female_and_CN = np.array(diagnosis == 'CN').astype(int) * np.array(sex == 'F').astype(int)
    female_and_MCI = np.array(diagnosis == 'MCI').astype(int) * np.array(sex == 'F').astype(int)
    male_and_AD = np.array(diagnosis == 'AD').astype(int) * np.array(sex == 'M').astype(int)
    male_and_CN = np.array(diagnosis == 'CN').astype(int) * np.array(sex == 'M').astype(int)
    male_and_MCI = np.array(diagnosis == 'MCI').astype(int) * np.array(sex == 'M').astype(int)

    print(' Female and AD = ' + str(sum(female_and_AD)))
    print(' Female and CN = ' + str(sum(female_and_CN)))
    print(' Female and MCI = ' + str(sum(female_and_MCI)))
    print(' male and AD = ' + str(sum(male_and_AD)))
    print(' male and CN = ' + str(sum(male_and_CN)))
    print(' male and MCI = ' + str(sum(male_and_MCI)))

# MMSE
    nan_mmse = mmse != mmse

    AD_mmse = np.mean(np.ma.masked_array(mmse, mask_non_AD))
    AD_mmse_std = np.std(np.ma.masked_array(mmse, mask_non_AD))
    CN_mmse = np.mean(np.ma.masked_array(mmse, np.asarray(mask_non_CN.astype(int)) + np.asarray(nan_mmse.astype(int))))
    CN_mmse_std = np.std(np.ma.masked_array(mmse, np.asarray(mask_non_CN.astype(int)) + np.asarray(nan_mmse.astype(int))))
    MCI_mmse = np.mean(np.ma.masked_array(mmse, np.asarray(mask_non_MCI.astype(int)) + np.asarray(nan_mmse.astype(int))))
    MCI_mmse_std = np.std(np.ma.masked_array(mmse, np.asarray(mask_non_MCI.astype(int)) + np.asarray(nan_mmse.astype(int))))

    print('AD mmse = ' + str(AD_mmse) + ' +/- ' + str(AD_mmse_std))
    print('CN mmse = ' + str(CN_mmse) + ' +/- ' + str(CN_mmse_std))
    print('MCI mmse = ' + str(MCI_mmse) + ' +/- ' + str(MCI_mmse_std))








