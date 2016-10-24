

def get_tissue_tuples(tissue_map, tissue_classes, dartel_tissues, save_warped_unmodulated, save_warped_modulated):
    '''
    Method to obtain the list of tuples, one for each tissue class, with the following fields:
         - tissue probability map (4D), 1-based index to frame
         - number of gaussians
         - which maps to save [Native, DARTEL] - a tuple of two boolean values
         - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values

    :param tissue_map: Path to tissue maps
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :return: List of tuples according to NewSegment input por tissues
    '''

    tissues = []

    for i in range(1, 7):
        n_gaussians = 2

        if i == 4 or i == 5:
            n_gaussians = i - 1

        native_space = False
        dartel_input = False
        warped_unmodulated = False
        warped_modulated = False

        if i in tissue_classes:
            native_space = True
            if save_warped_unmodulated:
                warped_unmodulated = True
            if save_warped_modulated:
                warped_modulated = True

        if i in dartel_tissues:
            dartel_input = True

        tissues.append(((tissue_map, i), n_gaussians, (native_space, dartel_input), (warped_unmodulated, warped_modulated)))

    return tissues


def get_class_images(class_images,index_list):
    """
    utility method to extract class images
    from a multi session class_images set:
    class nb : tissue type
    1 : Grey Matter
    2 : White Matter
    3 : CerebroSpinal Fluid
    4 : Skull
    5 : Out-of-brain soft tissue
    6 : Head surrounding

    :param class_images: image set from which to extract images.
    :param index_list: index list of the classes to extract.
    :return: extracted images in a list of list (without empty lists).

    Example
    -------
    >>> class_n_images = get_class_images(class_images,[1,2])
    """

    # Declare class images list
    class_n_images = {}
    for idx in index_list:
        class_n_images[idx]=[]

    for session in class_images:
        for idx in index_list:
            class_n_images[idx].extend(session[idx-1])

    result = []
    for idx in index_list:
        if class_n_images[idx]:
            result.append(class_n_images[idx])

    return result


def group_nested_images_by_subject(class_images):
    return [[s for tissue in subject for s in tissue] for subject in class_images]


def group_images_list_by_subject(images_list, tissue_classes):
    if len(tissue_classes) == 0:
        return []

    n_tissue_classes = len(tissue_classes)
    n_subjects = len(images_list)/n_tissue_classes
    grouped_images = []

    for subj in range(n_subjects):
        subj_images = []
        for tissue in range(n_tissue_classes):
            subj_images.append(images_list[(tissue * n_subjects) + subj])
        grouped_images.append(subj_images)

    return grouped_images

