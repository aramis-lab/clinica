def reading_image(directory):
    import nibabel as nib
    reading_input_image = nib.load(directory)
    image = reading_input_image.get_data()
    return image


def read_csv_label(atlas_id):
    import os
    import pandas
    directory_atlas=os.path.join(os.path.expandvars('$CLINICA_HOME'), 'clinica', 'resources', 'atlases_SPM')
    read_csv_atlas=os.path.join(directory_atlas,atlas_id+'.csv')
    label_list=[]

    if num_columns_csv(read_csv_atlas) > 5 :
        read_ROIname_voxel = pandas.read_csv(read_csv_atlas, sep=';', usecols=[2,5])# da cambiare
        read_labels=read_ROIname_voxel[read_ROIname_voxel.Voxel>0]
    elif ((num_columns_csv(read_csv_atlas)<5) and (num_columns_csv(read_csv_atlas)>2)):
        read_labels = pandas.read_csv(read_csv_atlas, sep=';', usecols=[2])
    else:
        read_labels = pandas.read_csv(read_csv_atlas, sep=';', usecols=[1])


    for j in range(len(read_labels)):
        label_list.append(read_labels.iloc[j,0])
    return label_list

def num_columns_csv(filename):
    import csv
    with open(filename, 'rU') as f1:

        csvlines = csv.reader(f1, delimiter=';')
        for lineNum, line in enumerate(csvlines):
            if lineNum == 0:
                headings = ';'.join(line)
    return len(line)

def statistics(atlas_id,image):
    import os
    import numpy as np
    import pandas
    directory_atlas=os.path.join(os.path.expandvars('$CLINICA_HOME'), 'clinica', 'resources', 'atlases_SPM')

    directory_atlas_image=os.path.join(directory_atlas, atlas_id+'.nii')

    atlas= reading_image(directory_atlas_image)

    labels = list(set(atlas.ravel()))

    stats_scalar=np.zeros((len(labels),2))
    label_list=read_csv_label(directory_atlas, atlas_id)
    #label_list=np.array(label_list,dtype="|S100")
    for index, n in enumerate(labels):

        atlas_label_index = np.array(np.where(atlas == n))
        stats_scalar[index,0]=index
        labeled_voxel = image[atlas_label_index[0, :], atlas_label_index[1, :], atlas_label_index[2, :]]
        average_voxel = labeled_voxel.mean()
        stats_scalar[index,1] = average_voxel

    data = pandas.DataFrame({'index': stats_scalar[:,0],
                    'label_name': label_list,
                    'mean_scalar': stats_scalar[:,1]
                   })


    return data


