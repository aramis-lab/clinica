from os import walk
from os.path import join
from pandas.io import parsers
from clinica.pipeline.machine_learning.voxel_based_svm import linear_svm_binary_classification


def get_all_subj_dx(diagnose_file, data_dir, dx_filter=None):

    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        filenames.sort()
        subjects = [x for x in filenames if x.endswith('.nii')]
        break

    labels = []
    new_subjects = []

    patient_list = parsers.read_csv(diagnose_file, sep=',')

    for i in range(len(subjects)):

        s = subjects[i]
        s_id = s[s.find('c1') + 3:].split('-')[0]
        res = patient_list[patient_list.S_ID == int(s_id)]
        if res.shape[0] < 1:
            continue
        patient = res.iloc[0]

        if dx_filter is not None:
            if patient.DX in dx_filter:
                new_subjects.append(subjects[i])
                labels.append(patient.DX)
        else:
            new_subjects.append(subjects[i])
            labels.append(patient.DX)

    return new_subjects, labels


input_dir = '/aramis/dataARAMIS/users/alexandre.morin/Outputs/segment_spm12/smooth_8'
diag_file = '/aramis/dataARAMIS/users/jorge.samper/alex.morin/DX.csv'
dx_filter = ['Depr', 'EOAD', 'LOAD', 'DCB', 'DCL', 'DFT', 'APPl', 'APPs']

filenames, diagnose_list = get_all_subj_dx(diag_file, input_dir, dx_filter)
image_list = [join(input_dir, f) for f in filenames]
output_dir = '/aramis/dataARAMIS/users/jorge.samper/test'

print 'Images: '
print len(image_list)
print 'Diagnosis'
print len(diagnose_list)

linear_svm_binary_classification(image_list, diagnose_list, output_dir, balanced=True, outer_folds=10, inner_folds=10, n_threads=40, save_subject_classification=True)

